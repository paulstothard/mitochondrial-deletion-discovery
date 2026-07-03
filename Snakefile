import csv
import os
from copy import deepcopy
from pathlib import Path

import yaml

os.environ.setdefault("MPLCONFIGDIR", str(Path(".snakemake") / "mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)


def deep_update(base, override):
    merged = deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_update(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


DEFAULTS = load_yaml("config/defaults.yaml")
CFG = deep_update(DEFAULTS, config)
DATASET_CONFIG = workflow.configfiles[-1] if workflow.configfiles else ""
DATASET = CFG["dataset"]["name"]
SPECIES = CFG["dataset"]["species"]
OUTDIR = f'{CFG["project"]["output_dir"]}/{DATASET}'
DELIVERABLES_DIR = f"{OUTDIR}/{DATASET}_deliverables"
WORKDIR = CFG["project"]["work_dir"]
REF = CFG["references"][SPECIES]
MT_LENGTH = int(REF["mt_length"])
MT_NAMES = ",".join(REF["mt_contig_names"])
RESOLVED_SAMPLES = f"metadata/generated/{DATASET}.samples.tsv"
RESOLVED_CONFIG = f"{OUTDIR}/config/resolved_config.yaml"
START_FROM = str(
    config.get("workflow_start_from")
    or CFG.get("workflow", {}).get("start_from")
    or "raw"
).strip().lower().replace("-", "_")
if START_FROM not in {"raw", "trimmed"}:
    raise ValueError("workflow.start_from/workflow_start_from must be 'raw' or 'trimmed'")

FIRST_PASS_SELECTION = str(
    config.get("first_pass_read_selection")
    or CFG.get("mapping", {}).get("first_pass_read_selection")
    or CFG.get("mt_realign", {}).get("input_strategy")
    or "mt_evidence_reads"
).strip().lower().replace("-", "_")
if FIRST_PASS_SELECTION not in {"mt_evidence_reads", "nuclear_unmapped_reads", "whole_genome_mt_best"}:
    raise ValueError("mapping.first_pass_read_selection must be 'mt_evidence_reads', 'nuclear_unmapped_reads', or 'whole_genome_mt_best'")

FIRST_PASS_ALIGNER = str(CFG.get("mapping", {}).get("first_pass_aligner", "star")).strip().lower()
if FIRST_PASS_SELECTION == "nuclear_unmapped_reads" and FIRST_PASS_ALIGNER not in {"star", "minimap2"}:
    raise ValueError("nuclear_unmapped_reads currently supports mapping.first_pass_aligner 'star' or 'minimap2'")
if FIRST_PASS_SELECTION == "whole_genome_mt_best" and FIRST_PASS_ALIGNER not in {"star", "minimap2"}:
    raise ValueError("whole_genome_mt_best currently supports mapping.first_pass_aligner 'star' or 'minimap2'")

# In trimmed mode, trimmed FASTQs are external inputs. Moving the producer-rule
# outputs out of the normal result paths prevents Snakemake from backtracking
# into downloads/trimming when raw FASTQs are absent.
RAW_RULE_DIR = "fastq" if START_FROM == "raw" else ".disabled/raw_fastq"
TRIM_RULE_DIR = "trimmed" if START_FROM == "raw" else ".disabled/trimmed"
READ_INPUT_RULE_DIR = "qc" if START_FROM == "raw" else ".disabled/qc_read_input"
FASTQC_RULE_DIR = "qc" if START_FROM == "raw" else ".disabled/qc_fastqc"
TRIM_QC_RULE_DIR = "qc" if START_FROM == "raw" else ".disabled/qc_trim"
CLASSIFY_MT_RULE_DIR = (
    "mt_reads"
    if FIRST_PASS_SELECTION == "mt_evidence_reads" and FIRST_PASS_ALIGNER == "star"
    else ".disabled/mt_reads_full_genome_evidence"
)
NUCLEAR_UNMAPPED_STAR_RULE_DIR = (
    "mt_reads" if FIRST_PASS_SELECTION == "nuclear_unmapped_reads" and FIRST_PASS_ALIGNER == "star" else ".disabled/mt_reads_nuclear_unmapped_star"
)
NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR = (
    "mt_reads" if FIRST_PASS_SELECTION == "nuclear_unmapped_reads" and FIRST_PASS_ALIGNER == "minimap2" else ".disabled/mt_reads_nuclear_unmapped_minimap2"
)
WHOLE_GENOME_MT_MINIMAP2_RULE_DIR = (
    "mt_reads" if FIRST_PASS_SELECTION == "whole_genome_mt_best" and FIRST_PASS_ALIGNER == "minimap2" else ".disabled/mt_reads_whole_genome_mt_minimap2"
)
WHOLE_GENOME_MT_STAR_RULE_DIR = (
    "mt_reads" if FIRST_PASS_SELECTION == "whole_genome_mt_best" and FIRST_PASS_ALIGNER == "star" else ".disabled/mt_reads_whole_genome_mt_star"
)


def parse_rotation_start(value):
    if isinstance(value, str) and value.strip().lower() == "half":
        return MT_LENGTH // 2 + 1
    return int(value)


MT_ROTATIONS = CFG.get("mt_realign", {}).get("rotations") or [
    {"name": "normal", "start": 1},
    {"name": "half", "start": MT_LENGTH // 2 + 1},
]
ROTATION_STARTS = {str(item["name"]): parse_rotation_start(item["start"]) for item in MT_ROTATIONS}
ROTATION_NAMES = list(ROTATION_STARTS)


def star_option_string(section, key="star_chimeric_options"):
    options = CFG.get(section, {}).get(key, {}) or {}
    parts = []
    for name, value in options.items():
        if value is False or value is None:
            continue
        if value is True:
            parts.append(f"--{name}")
        else:
            parts.append(f"--{name} {value}")
    return " ".join(parts)


def first_pass_aligner_is(name):
    return FIRST_PASS_SELECTION == "nuclear_unmapped_reads" and FIRST_PASS_ALIGNER == name


def read_sample_ids(sample_tsv):
    with open(sample_tsv, "r", encoding="utf-8", newline="") as handle:
        return [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]


def sample_ids(wildcards):
    resolved = checkpoints.resolve_samples.get().output.samples
    return read_sample_ids(resolved)


def sample_outputs(pattern):
    return lambda wildcards: expand(pattern, sample=sample_ids(wildcards))


def rotated_sample_outputs(pattern):
    return lambda wildcards: expand(pattern, sample=sample_ids(wildcards), rotation=ROTATION_NAMES)


def known_sequence_searches_configured():
    return bool(CFG.get("analysis", {}).get("known_sequence_searches", []) or [])


def known_sequence_r1_inputs(wildcards):
    if not known_sequence_searches_configured():
        return []
    return expand(f"{OUTDIR}/mt_reads/{{sample}}.mt_evidence.fastq.gz", sample=sample_ids(wildcards))


def known_sequence_r2_inputs(wildcards):
    return []


def known_sequence_count_inputs(wildcards):
    if not known_sequence_searches_configured():
        return []
    return expand(f"{OUTDIR}/mt_reads/{{sample}}.mt_read_summary.json", sample=sample_ids(wildcards))


def bool_flag(section, key, flag):
    return flag if CFG.get(section, {}).get(key, False) else ""


rule all:
    input:
        f"{OUTDIR}/junctions/junction_clusters.tsv",
        f"{OUTDIR}/annotations/mt_features.tsv",
        f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        f"{OUTDIR}/analysis/breakpoint_reference_support.tsv",
        f"{OUTDIR}/matrices/exact_deletion_raw_counts.tsv",
        f"{OUTDIR}/matrices/exact_deletion_support_per_million_mt_reads.tsv",
        f"{OUTDIR}/matrices/affected_feature_raw_counts.tsv",
        f"{OUTDIR}/matrices/affected_feature_support_per_million_mt_reads.tsv",
        f"{OUTDIR}/analysis/deletion_burden.tsv",
        f"{OUTDIR}/analysis/exact_deletion_comparison.tsv",
        f"{OUTDIR}/analysis/affected_feature_comparison.tsv",
        f"{OUTDIR}/analysis/feature_impact_class_comparison.tsv",
        f"{OUTDIR}/analysis/deletion_size_distribution_tests.tsv",
        f"{OUTDIR}/analysis/deletion_size_bin_summary.tsv",
        f"{OUTDIR}/analysis/factorial_model_summary.tsv",
        f"{OUTDIR}/analysis/deletion_metadata_associations.tsv",
        f"{OUTDIR}/analysis/per_gene_affected_burden.tsv",
        f"{OUTDIR}/analysis/qc_summary.tsv",
        f"{OUTDIR}/plots/deletion_burden_by_sample.pdf",
        f"{OUTDIR}/plots/unique_exact_deletions_by_sample.pdf",
        f"{OUTDIR}/plots/deletion_burden_factorial_interaction.pdf",
        f"{OUTDIR}/plots/unique_exact_deletions_factorial_interaction.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_unweighted.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_support_weighted.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_support_weighted_log_y.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_small.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_medium.pdf",
        f"{OUTDIR}/plots/deletion_size_distribution_large.pdf",
        f"{OUTDIR}/plots/deletion_rainfall_left_breakpoint.pdf",
        f"{OUTDIR}/plots/deletion_rainfall_right_breakpoint.pdf",
        f"{OUTDIR}/plots/deletion_rainfall_midpoint.pdf",
        f"{OUTDIR}/plots/breakpoint_pair_support_map.pdf",
        f"{OUTDIR}/plots/affected_feature_support.pdf",
        f"{OUTDIR}/plots/affected_feature_counts.pdf",
        f"{OUTDIR}/plots/affected_feature_proportions.pdf",
        f"{OUTDIR}/plots/feature_impact_classes.pdf",
        f"{OUTDIR}/plots/per_gene_affected_burden.pdf",
        f"{OUTDIR}/plots/exact_deletion_recurrence.pdf",
        f"{OUTDIR}/plots/exact_deletion_pca.pdf",
        f"{OUTDIR}/plots/exact_deletion_bray_curtis_mds.pdf",
        f"{OUTDIR}/plots/affected_feature_pca.pdf",
        f"{OUTDIR}/plots/affected_feature_bray_curtis_mds.pdf",
        f"{OUTDIR}/analysis/known_sequence_search_summary.tsv",
        f"{OUTDIR}/analysis/known_sequence_search_hits.tsv",
        f"{OUTDIR}/.report/index.html",
        f"{OUTDIR}/.report/read_lists/manifest.tsv",
        f"{DELIVERABLES_DIR}/DELIVERABLES_COMPLETE.txt",


checkpoint resolve_samples:
    input:
        defaults="config/defaults.yaml",
        dataset_config=DATASET_CONFIG,
    output:
        samples=RESOLVED_SAMPLES,
        config=RESOLVED_CONFIG,
        run_table=f"metadata/generated/{DATASET}.sra_run_table.csv",
    params:
        bioproject=lambda wildcards: f'--bioproject {CFG["dataset"]["bioproject"]}' if CFG["dataset"].get("bioproject") else "",
        sra_run_table=lambda wildcards: f'--sra-run-table {CFG["dataset"]["sra_run_table"]}' if CFG["dataset"].get("sra_run_table") else "",
        sample_table=lambda wildcards: f'--sample-table {CFG["dataset"]["sample_table"]}' if CFG["dataset"].get("sample_table") else "",
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/resolve_samples.py --defaults {input.defaults} --dataset-config {input.dataset_config} "
        "--dataset {DATASET} --species {SPECIES} {params.bioproject} {params.sra_run_table} {params.sample_table} "
        "--out-run-table {output.run_table} --out-samples {output.samples} --out-config {output.config}"


rule download_genome:
    output:
        fasta=f"{WORKDIR}/references/{SPECIES}/genome.fa",
    params:
        source=REF.get("genome_path") or REF.get("genome_url"),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/fetch_reference.py --source {params.source} --output {output.fasta}"


rule download_annotation:
    output:
        gtf=f"{WORKDIR}/references/{SPECIES}/annotation.gtf",
    params:
        source=REF.get("annotation_path") or REF.get("annotation_url"),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/fetch_reference.py --source {params.source} --output {output.gtf}"


rule extract_mt_reference:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/genome.fa",
    output:
        fasta=f"{WORKDIR}/references/{SPECIES}/mt.fa",
        json=f"{WORKDIR}/references/{SPECIES}/mt_reference.json",
    params:
        names=MT_NAMES,
        length=MT_LENGTH,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/extract_mt_reference.py --genome {input.fasta} --mt-contig-names {params.names} "
        "--expected-length {params.length} --out-fasta {output.fasta} --out-json {output.json}"


rule extract_nuclear_reference:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/genome.fa",
    output:
        fasta=f"{WORKDIR}/references/{SPECIES}/nuclear.fa",
        json=f"{WORKDIR}/references/{SPECIES}/nuclear_reference.json",
    params:
        names=MT_NAMES,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/extract_nuclear_reference.py --genome {input.fasta} --mt-contig-names {params.names} "
        "--out-fasta {output.fasta} --out-json {output.json}"


rule filter_nuclear_annotation:
    input:
        gtf=f"{WORKDIR}/references/{SPECIES}/annotation.gtf",
    output:
        gtf=f"{WORKDIR}/references/{SPECIES}/nuclear.annotation.gtf",
    params:
        names=MT_NAMES,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/filter_nuclear_gtf.py --gtf {input.gtf} --mt-contig-names {params.names} --out-gtf {output.gtf}"


rule extract_mt_features:
    input:
        gtf=f"{WORKDIR}/references/{SPECIES}/annotation.gtf",
    output:
        features=f"{OUTDIR}/annotations/mt_features.tsv",
    params:
        names=MT_NAMES,
        mt_length=MT_LENGTH,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/extract_mt_features.py --gtf {input.gtf} --mt-contig-names {params.names} "
        "--mt-length {params.mt_length} --output {output.features}"


rule make_rotated_mt_reference:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/mt.fa",
    output:
        fasta=f"{WORKDIR}/references/{SPECIES}/mt.{{rotation}}.fa",
        metadata=f"{WORKDIR}/references/{SPECIES}/mt.{{rotation}}.json",
    params:
        start=lambda wildcards: ROTATION_STARTS[wildcards.rotation],
        name=lambda wildcards: wildcards.rotation,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/make_rotated_mt_reference.py --input {input.fasta} --start {params.start} "
        "--name {params.name} --output {output.fasta} --metadata {output.metadata}"


rule index_full_genome:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/genome.fa",
        gtf=f"{WORKDIR}/references/{SPECIES}/annotation.gtf",
    output:
        directory(f"{WORKDIR}/indexes/{SPECIES}/star_full"),
    params:
        sjdb=48,
        sa=CFG["mapping"]["star_genome_sa_index_nbases"],
    threads:
        CFG["mapping"]["star_threads"]
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p {output} && STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb} --genomeSAindexNbases {params.sa}"


rule index_nuclear_star:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/nuclear.fa",
        gtf=f"{WORKDIR}/references/{SPECIES}/nuclear.annotation.gtf",
    output:
        directory(f"{WORKDIR}/indexes/{SPECIES}/star_nuclear"),
    params:
        sjdb=48,
        sa=CFG["mapping"]["star_genome_sa_index_nbases"],
    threads:
        CFG["mapping"]["star_threads"]
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p {output} && STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb} --genomeSAindexNbases {params.sa}"


rule index_nuclear_minimap2:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/nuclear.fa",
    output:
        mmi=f"{WORKDIR}/indexes/{SPECIES}/minimap2_nuclear.mmi",
    threads:
        2
    params:
        extra=CFG["mapping"].get("first_pass_minimap2_index_extra", ""),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.mmi}) && minimap2 {params.extra} -d {output.mmi} {input.fasta}"


rule index_full_minimap2:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/genome.fa",
    output:
        mmi=f"{WORKDIR}/indexes/{SPECIES}/minimap2_full.mmi",
    threads:
        2
    params:
        extra=CFG["mapping"].get("first_pass_minimap2_index_extra", ""),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.mmi}) && minimap2 {params.extra} -d {output.mmi} {input.fasta}"


rule index_rotated_mt:
    input:
        fasta=f"{WORKDIR}/references/{SPECIES}/mt.{{rotation}}.fa",
    output:
        mmi=f"{WORKDIR}/indexes/{SPECIES}/minimap2_mt_{{rotation}}.mmi",
    threads:
        2
    params:
        extra=CFG["mt_realign"].get("minimap2_index_extra", ""),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.mmi}) && minimap2 {params.extra} -d {output.mmi} {input.fasta}"


rule prepare_reads:
    input:
        samples=RESOLVED_SAMPLES,
    output:
        r1=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R2.fastq.gz",
        summary=f"{OUTDIR}/{READ_INPUT_RULE_DIR}/{{sample}}/read_input.json",
    log:
        log=f"{OUTDIR}/logs/downloads/{{sample}}.log",
    params:
        sample=lambda wildcards: wildcards.sample,
        method=CFG.get("downloads", {}).get("method", "fasterq_dump"),
        fasterq_threads=CFG.get("downloads", {}).get("fasterq_threads", 2),
        prefetch=CFG.get("downloads", {}).get("prefetch", False),
        dataset_config=DATASET_CONFIG,
    resources:
        download=1
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/prepare_reads.py --sample {params.sample} --sample-table {input.samples} "
        "--defaults config/defaults.yaml --config {params.dataset_config} "
        "--method {params.method} --fasterq-threads {params.fasterq_threads} "
        "--prefetch {params.prefetch} --out-r1 {output.r1} --out-r2 {output.r2} "
        "--summary {output.summary} --log {log.log}"


rule fastqc_raw:
    input:
        r1=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R2.fastq.gz",
    output:
        done=f"{OUTDIR}/{FASTQC_RULE_DIR}/{{sample}}/fastqc_raw.done",
    params:
        outdir=lambda wildcards: f"{OUTDIR}/qc/{wildcards.sample}/fastqc_raw",
        run_fastqc=CFG.get("qc", {}).get("run_fastqc", True),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "if [ '{params.run_fastqc}' = 'False' ]; then "
        "printf 'FastQC disabled by configuration\\n' > {params.outdir}/fastqc_disabled.txt; "
        "elif python scripts/fastq_gz_has_records.py {input.r2}; then fastqc -o {params.outdir} {input.r1} {input.r2}; "
        "else fastqc -o {params.outdir} {input.r1}; fi && "
        "touch {output.done}"


rule trim_reads:
    input:
        r1=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{RAW_RULE_DIR}/{{sample}}_R2.fastq.gz",
        fastqc=f"{OUTDIR}/{FASTQC_RULE_DIR}/{{sample}}/fastqc_raw.done",
    output:
        r1=f"{OUTDIR}/{TRIM_RULE_DIR}/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{TRIM_RULE_DIR}/{{sample}}_R2.fastq.gz",
        json=f"{OUTDIR}/{TRIM_QC_RULE_DIR}/{{sample}}/fastp.json",
        html=f"{OUTDIR}/{TRIM_QC_RULE_DIR}/{{sample}}/fastp.html",
        decision=f"{OUTDIR}/{TRIM_QC_RULE_DIR}/{{sample}}/qc_decision.json",
        counts=f"{OUTDIR}/{TRIM_QC_RULE_DIR}/{{sample}}/fragment_counts.tsv",
    params:
        min_len=CFG["qc"]["minimum_length_after_trimming"],
        extra=CFG["qc"].get("fastp_extra", ""),
        sample=lambda wildcards: wildcards.sample,
        skip="--skip" if not CFG.get("qc", {}).get("trim_reads", True) else "",
    threads:
        CFG["qc"].get("fastp_threads", 2)
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/run_fastp.py --sample {params.sample} --in-r1 {input.r1} --in-r2 {input.r2} "
        "--out-r1 {output.r1} --out-r2 {output.r2} --json {output.json} --html {output.html} "
        "--decision {output.decision} --counts {output.counts} --min-length {params.min_len} "
        "--threads {threads} --extra '{params.extra}' {params.skip}"


rule known_sequence_search:
    input:
        samples=lambda wildcards: checkpoints.resolve_samples.get().output.samples,
        r1=known_sequence_r1_inputs,
        r2=known_sequence_r2_inputs,
        counts=known_sequence_count_inputs,
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
    output:
        hits=f"{OUTDIR}/analysis/known_sequence_search_hits.tsv",
        summary=f"{OUTDIR}/analysis/known_sequence_search_summary.tsv",
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/search_known_sequences.py --defaults config/defaults.yaml --config {input.config} "
        "--samples {input.samples} --r1-files {input.r1} --r2-files {input.r2} "
        "--mt-summaries {input.counts} "
        "--hits {output.hits} --summary {output.summary}"


rule align_full_genome:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/star_full",
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        bam=f"{OUTDIR}/alignments/full/{{sample}}.bam",
        bai=f"{OUTDIR}/alignments/full/{{sample}}.bam.bai",
        log=f"{OUTDIR}/alignments/full/{{sample}}.Log.final.out",
        junction=f"{OUTDIR}/alignments/full/{{sample}}.Chimeric.out.junction",
    params:
        prefix=lambda wildcards: f"{OUTDIR}/alignments/full/.star_tmp/{wildcards.sample}.",
        star_options=star_option_string("mapping"),
    threads:
        CFG["mapping"]["star_threads"]
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.bam}) $(dirname {params.prefix}) && "
        "READ2=$(if python scripts/fastq_gz_has_records.py {input.r2}; then printf ' {input.r2}'; fi) && "
        "STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.r1} $READ2 "
        "--readFilesCommand 'gzip -cd' --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate "
        "{params.star_options} --twopassMode Basic && "
        "mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} && "
        "cp {params.prefix}Log.final.out {output.log} && "
        "if [ -f {params.prefix}Chimeric.out.junction ]; then mv {params.prefix}Chimeric.out.junction {output.junction}; else touch {output.junction}; fi && "
        "samtools index {output.bam}"


rule classify_mt_reads:
    input:
        bam=f"{OUTDIR}/alignments/full/{{sample}}.bam",
        bai=f"{OUTDIR}/alignments/full/{{sample}}.bam.bai",
        junction=f"{OUTDIR}/alignments/full/{{sample}}.Chimeric.out.junction",
    output:
        fastq=f"{OUTDIR}/{CLASSIFY_MT_RULE_DIR}/{{sample}}.high_confidence_mt.fastq.gz",
        ambiguous=f"{OUTDIR}/{CLASSIFY_MT_RULE_DIR}/{{sample}}.ambiguous_mt.fastq.gz",
        evidence=f"{OUTDIR}/{CLASSIFY_MT_RULE_DIR}/{{sample}}.mt_evidence.fastq.gz",
        tsv=f"{OUTDIR}/{CLASSIFY_MT_RULE_DIR}/{{sample}}.mt_read_classification.tsv",
        summary=f"{OUTDIR}/{CLASSIFY_MT_RULE_DIR}/{{sample}}.mt_read_summary.json",
    params:
        names=MT_NAMES,
        mapq=CFG["mapping"]["minimum_mapq_full_genome"],
        include_low_mapq=bool_flag("mt_realign", "include_low_mapq", "--include-low-mapq"),
        include_multimappers=bool_flag("mt_realign", "include_multimappers", "--include-multimappers"),
        include_supplementary=bool_flag("mt_realign", "include_supplementary", "--include-supplementary"),
        include_secondary=bool_flag("mt_realign", "include_secondary", "--include-secondary"),
        include_chimeric=bool_flag("mt_realign", "include_chimeric_mt_reads", "--include-chimeric-mt-reads"),
        include_mates=bool_flag("mt_realign", "include_mates_of_mt_evidence_reads", "--include-mates-of-mt-evidence-reads"),
        write_tsv=bool_flag("mt_realign", "write_read_classification_tsv", "--write-read-classification-tsv"),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/classify_mt_reads.py --bam {input.bam} --mt-contig-names {params.names} "
        "--min-mapq {params.mapq} --chimeric-junction {input.junction} "
        "--high-confidence-fastq {output.fastq} --ambiguous-fastq {output.ambiguous} "
        "--mt-evidence-fastq {output.evidence} --classification {output.tsv} --summary {output.summary} "
        "{params.include_low_mapq} {params.include_multimappers} {params.include_supplementary} "
        "{params.include_secondary} {params.include_chimeric} {params.include_mates} {params.write_tsv}"


rule select_whole_genome_mt_star:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/star_full",
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        fastq=f"{OUTDIR}/{WHOLE_GENOME_MT_STAR_RULE_DIR}/{{sample}}.high_confidence_mt.fastq.gz",
        ambiguous=f"{OUTDIR}/{WHOLE_GENOME_MT_STAR_RULE_DIR}/{{sample}}.ambiguous_mt.fastq.gz",
        evidence=f"{OUTDIR}/{WHOLE_GENOME_MT_STAR_RULE_DIR}/{{sample}}.mt_evidence.fastq.gz",
        tsv=f"{OUTDIR}/{WHOLE_GENOME_MT_STAR_RULE_DIR}/{{sample}}.mt_read_classification.tsv",
        summary=f"{OUTDIR}/{WHOLE_GENOME_MT_STAR_RULE_DIR}/{{sample}}.mt_read_summary.json",
        log=f"{OUTDIR}/alignments/full_stream/{{sample}}.Log.final.out",
    params:
        names=MT_NAMES,
        sample=lambda wildcards: wildcards.sample,
        prefix=lambda wildcards: f"{OUTDIR}/alignments/full_stream/.star_tmp/{wildcards.sample}.",
        tmp=lambda wildcards: f"{OUTDIR}/alignments/full_stream/.collate_tmp/{wildcards.sample}",
        star_options=star_option_string("mapping"),
        min_mt_mapq=CFG["mapping"].get("whole_genome_min_mt_mapq", 0),
        min_mt_aligned_fraction=CFG["mapping"].get("whole_genome_min_mt_aligned_fraction", 0.5),
        ambiguous_mapq_below=CFG["mapping"].get("whole_genome_ambiguous_mapq_below", 10),
        competing_nuclear_aligned_fraction=CFG["mapping"].get("whole_genome_competing_nuclear_aligned_fraction", 0.5),
        keep_ambiguous=bool_flag("mapping", "keep_ambiguous_mt_nuclear_reads", "--keep-ambiguous-mt-nuclear"),
    threads:
        max(2, CFG["mapping"].get("star_threads", 8))
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.evidence}) $(dirname {output.log}) $(dirname {params.prefix}) {params.tmp} && "
        "READ2=$(if python scripts/fastq_gz_has_records.py {input.r2}; then printf ' {input.r2}'; fi) && "
        "STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.r1} $READ2 "
        "--readFilesCommand 'gzip -cd' --outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM Unsorted --outStd BAM_Unsorted {params.star_options} --twopassMode Basic | "
        "samtools collate -@ {threads} -u -O -T {params.tmp}/collate - | "
        "python scripts/select_whole_genome_mt_from_sam.py --sample {params.sample} --mt-contig-names {params.names} "
        "--input-format bam "
        "--mt-evidence-fastq {output.evidence} --high-confidence-fastq {output.fastq} "
        "--ambiguous-fastq {output.ambiguous} --classification {output.tsv} --summary {output.summary} "
        "--min-mt-mapq {params.min_mt_mapq} --min-mt-aligned-fraction {params.min_mt_aligned_fraction} "
        "--ambiguous-mapq-below {params.ambiguous_mapq_below} "
        "--competing-nuclear-aligned-fraction {params.competing_nuclear_aligned_fraction} {params.keep_ambiguous} && "
        "cp {params.prefix}Log.final.out {output.log}"


rule select_nuclear_unmapped_star:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/star_nuclear",
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        fastq=f"{OUTDIR}/{NUCLEAR_UNMAPPED_STAR_RULE_DIR}/{{sample}}.high_confidence_mt.fastq.gz",
        ambiguous=f"{OUTDIR}/{NUCLEAR_UNMAPPED_STAR_RULE_DIR}/{{sample}}.ambiguous_mt.fastq.gz",
        evidence=f"{OUTDIR}/{NUCLEAR_UNMAPPED_STAR_RULE_DIR}/{{sample}}.mt_evidence.fastq.gz",
        tsv=f"{OUTDIR}/{NUCLEAR_UNMAPPED_STAR_RULE_DIR}/{{sample}}.mt_read_classification.tsv",
        summary=f"{OUTDIR}/{NUCLEAR_UNMAPPED_STAR_RULE_DIR}/{{sample}}.mt_read_summary.json",
        log=f"{OUTDIR}/alignments/nuclear/{{sample}}.Log.final.out",
    params:
        prefix=lambda wildcards: f"{OUTDIR}/alignments/nuclear/.star_tmp/{wildcards.sample}.",
        sample=lambda wildcards: wildcards.sample,
        star_options=star_option_string("mapping", "star_nuclear_options"),
    threads:
        CFG["mapping"]["star_threads"]
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.evidence}) $(dirname {output.log}) $(dirname {params.prefix}) && "
        "READ2=$(if python scripts/fastq_gz_has_records.py {input.r2}; then printf ' {input.r2}'; fi) && "
        "STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.r1} $READ2 "
        "--readFilesCommand 'gzip -cd' --outFileNamePrefix {params.prefix} --outSAMtype None "
        "--outReadsUnmapped Fastx {params.star_options} && "
        "cp {params.prefix}Log.final.out {output.log} && "
        "python scripts/collect_nuclear_unmapped_fastq.py --sample {params.sample} --source star_nuclear_unmapped "
        "--mate1 {params.prefix}Unmapped.out.mate1 --mate2 {params.prefix}Unmapped.out.mate2 "
        "--mt-evidence-fastq {output.evidence} --high-confidence-fastq {output.fastq} "
        "--ambiguous-fastq {output.ambiguous} --classification {output.tsv} --summary {output.summary} "
        "--threads {threads}"


rule select_nuclear_unmapped_minimap2:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/minimap2_nuclear.mmi",
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        fastq=f"{OUTDIR}/{NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR}/{{sample}}.high_confidence_mt.fastq.gz",
        ambiguous=f"{OUTDIR}/{NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR}/{{sample}}.ambiguous_mt.fastq.gz",
        evidence=f"{OUTDIR}/{NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR}/{{sample}}.mt_evidence.fastq.gz",
        tsv=f"{OUTDIR}/{NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR}/{{sample}}.mt_read_classification.tsv",
        summary=f"{OUTDIR}/{NUCLEAR_UNMAPPED_MINIMAP2_RULE_DIR}/{{sample}}.mt_read_summary.json",
    params:
        sample=lambda wildcards: wildcards.sample,
        preset=CFG["mapping"].get("first_pass_minimap2_preset", "sr"),
        extra=CFG["mapping"].get("first_pass_minimap2_extra", ""),
    threads:
        CFG["mapping"].get("first_pass_minimap2_threads", CFG["mt_realign"].get("minimap2_threads", 4))
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.evidence}) && "
        "READ2=$(if python scripts/fastq_gz_has_records.py {input.r2}; then printf ' {input.r2}'; fi) && "
        "minimap2 -t {threads} -ax {params.preset} {params.extra} {input.idx} {input.r1} $READ2 | "
        "python scripts/select_nuclear_unmapped_from_sam.py --sample {params.sample} --source minimap2_nuclear_unmapped "
        "--mt-evidence-fastq {output.evidence} --high-confidence-fastq {output.fastq} "
        "--ambiguous-fastq {output.ambiguous} --classification {output.tsv} --summary {output.summary}"


rule select_whole_genome_mt_minimap2:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/minimap2_full.mmi",
        r1=f"{OUTDIR}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        fastq=f"{OUTDIR}/{WHOLE_GENOME_MT_MINIMAP2_RULE_DIR}/{{sample}}.high_confidence_mt.fastq.gz",
        ambiguous=f"{OUTDIR}/{WHOLE_GENOME_MT_MINIMAP2_RULE_DIR}/{{sample}}.ambiguous_mt.fastq.gz",
        evidence=f"{OUTDIR}/{WHOLE_GENOME_MT_MINIMAP2_RULE_DIR}/{{sample}}.mt_evidence.fastq.gz",
        tsv=f"{OUTDIR}/{WHOLE_GENOME_MT_MINIMAP2_RULE_DIR}/{{sample}}.mt_read_classification.tsv",
        summary=f"{OUTDIR}/{WHOLE_GENOME_MT_MINIMAP2_RULE_DIR}/{{sample}}.mt_read_summary.json",
    params:
        sample=lambda wildcards: wildcards.sample,
        preset=CFG["mapping"].get("first_pass_minimap2_preset", "sr"),
        extra=CFG["mapping"].get("first_pass_minimap2_extra", ""),
        mt_names=MT_NAMES,
        min_mt_mapq=CFG["mapping"].get("whole_genome_min_mt_mapq", 0),
        min_mt_aligned_fraction=CFG["mapping"].get("whole_genome_min_mt_aligned_fraction", 0.5),
        ambiguous_mapq_below=CFG["mapping"].get("whole_genome_ambiguous_mapq_below", 10),
        competing_nuclear_aligned_fraction=CFG["mapping"].get("whole_genome_competing_nuclear_aligned_fraction", 0.5),
        keep_ambiguous=bool_flag("mapping", "keep_ambiguous_mt_nuclear_reads", "--keep-ambiguous-mt-nuclear"),
    threads:
        CFG["mapping"].get("first_pass_minimap2_threads", CFG["mt_realign"].get("minimap2_threads", 4))
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.evidence}) && "
        "READ2=$(if python scripts/fastq_gz_has_records.py {input.r2}; then printf ' {input.r2}'; fi) && "
        "minimap2 -t {threads} -ax {params.preset} --secondary=yes {params.extra} {input.idx} {input.r1} $READ2 | "
        "python scripts/select_whole_genome_mt_from_sam.py --sample {params.sample} "
        "--mt-contig-names {params.mt_names} --min-mt-mapq {params.min_mt_mapq} "
        "--min-mt-aligned-fraction {params.min_mt_aligned_fraction} "
        "--ambiguous-mapq-below {params.ambiguous_mapq_below} "
        "--competing-nuclear-aligned-fraction {params.competing_nuclear_aligned_fraction} "
        "{params.keep_ambiguous} --mt-evidence-fastq {output.evidence} "
        "--high-confidence-fastq {output.fastq} --ambiguous-fastq {output.ambiguous} "
        "--classification {output.tsv} --summary {output.summary}"


rule realign_mt_reads:
    input:
        idx=f"{WORKDIR}/indexes/{SPECIES}/minimap2_mt_{{rotation}}.mmi",
        fastq=f"{OUTDIR}/mt_reads/{{sample}}.mt_evidence.fastq.gz",
    output:
        bam=f"{OUTDIR}/mt_minimap2/{{rotation}}/{{sample}}.bam",
        bai=f"{OUTDIR}/mt_minimap2/{{rotation}}/{{sample}}.bam.bai",
    params:
        preset=CFG["mt_realign"].get("minimap2_preset", "sr"),
        extra=CFG["mt_realign"].get("minimap2_extra", ""),
    threads:
        CFG["mt_realign"].get("minimap2_threads", 4)
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "mkdir -p $(dirname {output.bam}) && "
        "minimap2 -t {threads} -ax {params.preset} --secondary=yes {params.extra} {input.idx} {input.fastq} | "
        "samtools sort -@ {threads} -o {output.bam} && "
        "samtools index {output.bam}"


rule parse_split_alignments:
    input:
        bam=f"{OUTDIR}/mt_minimap2/{{rotation}}/{{sample}}.bam",
    output:
        candidates=f"{OUTDIR}/deletions/rotated/{{rotation}}/{{sample}}.candidate_deletion_reads.tsv",
        filtered=f"{OUTDIR}/deletions/rotated/{{rotation}}/{{sample}}.filtered_deletion_reads.tsv",
        summary=f"{OUTDIR}/deletions/rotated/{{rotation}}/{{sample}}.deletion_summary.tsv",
    params:
        sample=lambda wildcards: wildcards.sample,
        species=SPECIES,
        mt_length=MT_LENGTH,
        padding=CFG["mt_realign"].get("circular_padding", 0),
        rotation_start=lambda wildcards: ROTATION_STARTS[wildcards.rotation],
        rotation_name=lambda wildcards: wildcards.rotation,
        min_anchor=CFG["junctions"]["min_anchor_length"],
        min_del=CFG["junctions"]["min_deletion_size"],
        max_del=CFG["junctions"]["max_deletion_size"],
        min_mapq=CFG["mt_realign"].get("minimap2_min_mapq", 0),
        min_aligned_fraction=CFG["mt_realign"].get("min_segment_aligned_fraction", 0.15),
        max_soft_clip_fraction=CFG["mt_realign"].get("max_soft_clip_fraction", 0.9),
        max_query_overlap=CFG["mt_realign"].get("max_query_overlap_bp", 10),
        max_query_gap=CFG["mt_realign"].get("max_query_gap_bp", 20),
        include_secondary=bool_flag("mt_realign", "minimap2_include_secondary", "--include-secondary"),
        include_supplementary=bool_flag("mt_realign", "minimap2_include_supplementary", "--include-supplementary"),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/call_minimap2_deletions.py --sample {params.sample} --species {params.species} "
        "--bam {input.bam} --mt-length {params.mt_length} "
        "--rotation-start {params.rotation_start} --rotation-name {params.rotation_name} "
        "--min-anchor-length {params.min_anchor} --min-deletion-size {params.min_del} --max-deletion-size {params.max_del} "
        "--min-mapq {params.min_mapq} --min-segment-aligned-fraction {params.min_aligned_fraction} "
        "--max-soft-clip-fraction {params.max_soft_clip_fraction} --max-query-overlap-bp {params.max_query_overlap} "
        "--max-query-gap-bp {params.max_query_gap} {params.include_secondary} {params.include_supplementary} "
        "--candidates {output.candidates} --filtered {output.filtered} --summary {output.summary}"


rule cluster_junctions:
    input:
        rotated_sample_outputs(f"{OUTDIR}/deletions/rotated/{{rotation}}/{{sample}}.filtered_deletion_reads.tsv")
    output:
        all_reads=f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        clusters=f"{OUTDIR}/junctions/junction_clusters.unannotated.tsv",
        id_map=f"{OUTDIR}/junctions/junction_id_map.tsv",
    params:
        slop=CFG["junctions"]["breakpoint_slop_bp"],
        min_support=CFG["junctions"]["min_split_read_support"],
        mt_length=MT_LENGTH,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/consolidate_deletions.py --slop {params.slop} --min-support {params.min_support} "
        "--mt-length {params.mt_length} "
        "--all-reads {output.all_reads} --clusters {output.clusters} --id-map {output.id_map} {input}"


rule estimate_breakpoint_reference_support:
    input:
        clusters=f"{OUTDIR}/junctions/junction_clusters.unannotated.tsv",
        all_reads=f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        bams=rotated_sample_outputs(f"{OUTDIR}/mt_minimap2/{{rotation}}/{{sample}}.bam"),
    output:
        clusters=f"{OUTDIR}/junctions/junction_clusters.with_reference_support.tsv",
        support=f"{OUTDIR}/analysis/breakpoint_reference_support.tsv",
    params:
        mt_length=MT_LENGTH,
        rotation_starts=",".join(f"{name}:{start}" for name, start in ROTATION_STARTS.items()),
        window=CFG["junctions"].get("reference_support_window_bp", 20),
        min_mapq=CFG["mt_realign"].get("minimap2_min_mapq", 0),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/estimate_breakpoint_reference_support.py --clusters {input.clusters} "
        "--all-reads {input.all_reads} --bam {input.bams} --rotation-starts {params.rotation_starts} "
        "--mt-length {params.mt_length} --window-bp {params.window} --min-mapq {params.min_mapq} "
        "--out-clusters {output.clusters} --out-reference-support {output.support}"


rule annotate_junctions:
    input:
        clusters=f"{OUTDIR}/junctions/junction_clusters.with_reference_support.tsv",
        features=f"{OUTDIR}/annotations/mt_features.tsv",
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
    output:
        clusters=f"{OUTDIR}/junctions/junction_clusters.tsv",
    params:
        mt_length=MT_LENGTH,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/annotate_junctions.py --clusters {input.clusters} --features {input.features} "
        "--mt-length {params.mt_length} --config {input.config} --output {output.clusters}"


rule build_matrices:
    input:
        samples=lambda wildcards: checkpoints.resolve_samples.get().output.samples,
        clusters=f"{OUTDIR}/junctions/junction_clusters.tsv",
        id_map=f"{OUTDIR}/junctions/junction_id_map.tsv",
        all_reads=f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
        counts=sample_outputs(f"{OUTDIR}/qc/{{sample}}/fragment_counts.tsv"),
        mt_summaries=sample_outputs(f"{OUTDIR}/mt_reads/{{sample}}.mt_read_summary.json"),
    output:
        exact_raw=f"{OUTDIR}/matrices/exact_deletion_raw_counts.tsv",
        exact_mtpm=f"{OUTDIR}/matrices/exact_deletion_support_per_million_mt_reads.tsv",
        affected_raw=f"{OUTDIR}/matrices/affected_feature_raw_counts.tsv",
        affected_mtpm=f"{OUTDIR}/matrices/affected_feature_support_per_million_mt_reads.tsv",
        impact_raw=f"{OUTDIR}/matrices/feature_impact_class_raw_counts.tsv",
        impact_mtpm=f"{OUTDIR}/matrices/feature_impact_class_support_per_million_mt_reads.tsv",
        per_gene=f"{OUTDIR}/analysis/per_gene_affected_burden.tsv",
        burden=f"{OUTDIR}/analysis/deletion_burden.tsv",
        exact_comparison=f"{OUTDIR}/analysis/exact_deletion_comparison.tsv",
        affected_comparison=f"{OUTDIR}/analysis/affected_feature_comparison.tsv",
        impact_comparison=f"{OUTDIR}/analysis/feature_impact_class_comparison.tsv",
        size_tests=f"{OUTDIR}/analysis/deletion_size_distribution_tests.tsv",
        size_bin_summary=f"{OUTDIR}/analysis/deletion_size_bin_summary.tsv",
        factorial_model_summary=f"{OUTDIR}/analysis/factorial_model_summary.tsv",
        metadata_assoc=f"{OUTDIR}/analysis/deletion_metadata_associations.tsv",
        qc_summary=f"{OUTDIR}/analysis/qc_summary.tsv",
    params:
        group=CFG["dataset"].get("primary_group_column", ""),
        group_columns=",".join(CFG["dataset"].get("group_columns", [])),
        normalization_denominator=CFG["analysis"].get("normalization_denominator", "total_usable_reads"),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/analyze_deletions.py --samples {input.samples} --clusters {input.clusters} "
        "--id-map {input.id_map} --all-reads {input.all_reads} --config {input.config} "
        "--group-column {params.group} --group-columns {params.group_columns} "
        "--normalization-denominator {params.normalization_denominator} "
        "--fragment-counts {input.counts} --mt-summaries {input.mt_summaries} "
        "--out-exact-raw {output.exact_raw} --out-exact-mtpm {output.exact_mtpm} "
        "--out-affected-raw {output.affected_raw} --out-affected-mtpm {output.affected_mtpm} "
        "--out-impact-class-raw {output.impact_raw} --out-impact-class-mtpm {output.impact_mtpm} "
        "--out-per-gene-burden {output.per_gene} --out-burden {output.burden} "
        "--out-exact-comparison {output.exact_comparison} --out-affected-comparison {output.affected_comparison} "
        "--out-impact-class-comparison {output.impact_comparison} --out-size-tests {output.size_tests} "
        "--out-size-bin-summary {output.size_bin_summary} --out-factorial-model-summary {output.factorial_model_summary} "
        "--out-metadata-associations {output.metadata_assoc} --out-qc-summary {output.qc_summary}"


rule plot_results:
    input:
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
        samples=lambda wildcards: checkpoints.resolve_samples.get().output.samples,
        features=f"{OUTDIR}/annotations/mt_features.tsv",
        all_reads=f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        clusters=f"{OUTDIR}/junctions/junction_clusters.tsv",
        burden=f"{OUTDIR}/analysis/deletion_burden.tsv",
        exact_mtpm=f"{OUTDIR}/matrices/exact_deletion_support_per_million_mt_reads.tsv",
        affected_raw=f"{OUTDIR}/matrices/affected_feature_raw_counts.tsv",
        affected_mtpm=f"{OUTDIR}/matrices/affected_feature_support_per_million_mt_reads.tsv",
        impact_mtpm=f"{OUTDIR}/matrices/feature_impact_class_support_per_million_mt_reads.tsv",
        per_gene=f"{OUTDIR}/analysis/per_gene_affected_burden.tsv",
        exact_comparison=f"{OUTDIR}/analysis/exact_deletion_comparison.tsv",
    output:
        burden=f"{OUTDIR}/plots/deletion_burden_by_sample.pdf",
        unique=f"{OUTDIR}/plots/unique_exact_deletions_by_sample.pdf",
        burden_factorial=f"{OUTDIR}/plots/deletion_burden_factorial_interaction.pdf",
        unique_factorial=f"{OUTDIR}/plots/unique_exact_deletions_factorial_interaction.pdf",
        size_unweighted=f"{OUTDIR}/plots/deletion_size_distribution_unweighted.pdf",
        size_weighted=f"{OUTDIR}/plots/deletion_size_distribution_support_weighted.pdf",
        size_weighted_log=f"{OUTDIR}/plots/deletion_size_distribution_support_weighted_log_y.pdf",
        size_small=f"{OUTDIR}/plots/deletion_size_distribution_small.pdf",
        size_medium=f"{OUTDIR}/plots/deletion_size_distribution_medium.pdf",
        size_large=f"{OUTDIR}/plots/deletion_size_distribution_large.pdf",
        rainfall_left=f"{OUTDIR}/plots/deletion_rainfall_left_breakpoint.pdf",
        rainfall_right=f"{OUTDIR}/plots/deletion_rainfall_right_breakpoint.pdf",
        rainfall_midpoint=f"{OUTDIR}/plots/deletion_rainfall_midpoint.pdf",
        breakpoint_pair_map=f"{OUTDIR}/plots/breakpoint_pair_support_map.pdf",
        affected_support=f"{OUTDIR}/plots/affected_feature_support.pdf",
        affected_counts=f"{OUTDIR}/plots/affected_feature_counts.pdf",
        affected_proportions=f"{OUTDIR}/plots/affected_feature_proportions.pdf",
        impact=f"{OUTDIR}/plots/feature_impact_classes.pdf",
        per_gene=f"{OUTDIR}/plots/per_gene_affected_burden.pdf",
        recurrence=f"{OUTDIR}/plots/exact_deletion_recurrence.pdf",
        exact_pca=f"{OUTDIR}/plots/exact_deletion_pca.pdf",
        exact_mds=f"{OUTDIR}/plots/exact_deletion_bray_curtis_mds.pdf",
        affected_pca=f"{OUTDIR}/plots/affected_feature_pca.pdf",
        affected_mds=f"{OUTDIR}/plots/affected_feature_bray_curtis_mds.pdf",
    params:
        group=CFG["dataset"].get("primary_group_column", ""),
        rainfall_min_support_per_million=CFG.get("plots", {}).get("rainfall_min_support_per_million", 1.0),
        rainfall_max_points_per_group=CFG.get("plots", {}).get("rainfall_max_points_per_group", 0),
        mt_length=MT_LENGTH,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/plot_deletion_results.py --samples {input.samples} --features {input.features} "
        "--config {input.config} --mt-length {params.mt_length} "
        "--all-reads {input.all_reads} --clusters {input.clusters} --burden {input.burden} "
        "--exact-mtpm {input.exact_mtpm} --affected-raw {input.affected_raw} --affected-mtpm {input.affected_mtpm} "
        "--impact-class-mtpm {input.impact_mtpm} --per-gene-burden {input.per_gene} "
        "--exact-comparison {input.exact_comparison} --group-column {params.group} "
        "--out-burden {output.burden} --out-unique-count {output.unique} "
        "--out-burden-factorial {output.burden_factorial} --out-unique-factorial {output.unique_factorial} "
        "--out-size-unweighted {output.size_unweighted} --out-size-weighted {output.size_weighted} "
        "--out-size-weighted-log {output.size_weighted_log} --out-size-small {output.size_small} "
        "--out-size-medium {output.size_medium} --out-size-large {output.size_large} "
        "--out-rainfall-left {output.rainfall_left} --out-rainfall-right {output.rainfall_right} "
        "--out-rainfall-midpoint {output.rainfall_midpoint} --out-breakpoint-pair-map {output.breakpoint_pair_map} "
        "--out-affected-support {output.affected_support} "
        "--out-affected-counts {output.affected_counts} --out-affected-proportions {output.affected_proportions} "
        "--out-impact-class {output.impact} --out-per-gene {output.per_gene} --out-exact-recurrence {output.recurrence} "
        "--out-exact-pca {output.exact_pca} --out-exact-mds {output.exact_mds} "
        "--out-affected-pca {output.affected_pca} --out-affected-mds {output.affected_mds} "
        "--rainfall-min-support-per-million {params.rainfall_min_support_per_million} "
        "--rainfall-max-points-per-group {params.rainfall_max_points_per_group}"


rule make_report:
    input:
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
        samples=lambda wildcards: checkpoints.resolve_samples.get().output.samples,
        features=f"{OUTDIR}/annotations/mt_features.tsv",
        clusters=f"{OUTDIR}/junctions/junction_clusters.tsv",
        junction_reads=f"{OUTDIR}/junctions/all_samples.filtered_junction_reads.tsv",
        qc_summary=f"{OUTDIR}/analysis/qc_summary.tsv",
        burden=f"{OUTDIR}/analysis/deletion_burden.tsv",
        exact_comparison=f"{OUTDIR}/analysis/exact_deletion_comparison.tsv",
        affected_comparison=f"{OUTDIR}/analysis/affected_feature_comparison.tsv",
        impact_comparison=f"{OUTDIR}/analysis/feature_impact_class_comparison.tsv",
        size_tests=f"{OUTDIR}/analysis/deletion_size_distribution_tests.tsv",
        size_bin_summary=f"{OUTDIR}/analysis/deletion_size_bin_summary.tsv",
        factorial_model_summary=f"{OUTDIR}/analysis/factorial_model_summary.tsv",
        metadata_assoc=f"{OUTDIR}/analysis/deletion_metadata_associations.tsv",
        per_gene=f"{OUTDIR}/analysis/per_gene_affected_burden.tsv",
        known_sequence_summary=f"{OUTDIR}/analysis/known_sequence_search_summary.tsv",
        known_sequence_hits=f"{OUTDIR}/analysis/known_sequence_search_hits.tsv",
        plots=[
            f"{OUTDIR}/plots/deletion_burden_by_sample.pdf",
            f"{OUTDIR}/plots/unique_exact_deletions_by_sample.pdf",
            f"{OUTDIR}/plots/deletion_burden_factorial_interaction.pdf",
            f"{OUTDIR}/plots/unique_exact_deletions_factorial_interaction.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_unweighted.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_support_weighted.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_support_weighted_log_y.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_small.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_medium.pdf",
            f"{OUTDIR}/plots/deletion_size_distribution_large.pdf",
            f"{OUTDIR}/plots/deletion_rainfall_left_breakpoint.pdf",
            f"{OUTDIR}/plots/deletion_rainfall_right_breakpoint.pdf",
            f"{OUTDIR}/plots/deletion_rainfall_midpoint.pdf",
            f"{OUTDIR}/plots/breakpoint_pair_support_map.pdf",
            f"{OUTDIR}/plots/affected_feature_support.pdf",
            f"{OUTDIR}/plots/affected_feature_counts.pdf",
            f"{OUTDIR}/plots/affected_feature_proportions.pdf",
            f"{OUTDIR}/plots/feature_impact_classes.pdf",
            f"{OUTDIR}/plots/per_gene_affected_burden.pdf",
            f"{OUTDIR}/plots/exact_deletion_recurrence.pdf",
            f"{OUTDIR}/plots/exact_deletion_pca.pdf",
            f"{OUTDIR}/plots/exact_deletion_bray_curtis_mds.pdf",
            f"{OUTDIR}/plots/affected_feature_pca.pdf",
            f"{OUTDIR}/plots/affected_feature_bray_curtis_mds.pdf",
        ],
    output:
        html=f"{OUTDIR}/.report/index.html",
        read_list_manifest=f"{OUTDIR}/.report/read_lists/manifest.tsv",
    params:
        title=CFG["dataset"].get("title", DATASET),
        group=CFG["dataset"].get("primary_group_column", ""),
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/make_deletion_report.py --title '{params.title}' --config {input.config} "
        "--samples {input.samples} --features {input.features} --qc-summary {input.qc_summary} "
        "--clusters {input.clusters} --junction-reads {input.junction_reads} --burden {input.burden} --exact-comparison {input.exact_comparison} "
        "--affected-comparison {input.affected_comparison} --impact-class-comparison {input.impact_comparison} "
        "--size-tests {input.size_tests} --size-bin-summary {input.size_bin_summary} "
        "--factorial-model-summary {input.factorial_model_summary} --metadata-associations {input.metadata_assoc} "
        "--per-gene-burden {input.per_gene} --known-sequence-summary {input.known_sequence_summary} "
        "--known-sequence-hits {input.known_sequence_hits} "
        "--read-list-manifest {output.read_list_manifest} "
        "--plots {input.plots} --output {output.html}"


rule make_deliverables:
    input:
        report=f"{OUTDIR}/.report/index.html",
        read_list_manifest=f"{OUTDIR}/.report/read_lists/manifest.tsv",
        config=DATASET_CONFIG if DATASET_CONFIG else RESOLVED_CONFIG,
        clusters=f"{OUTDIR}/junctions/junction_clusters.tsv",
        burden=f"{OUTDIR}/analysis/deletion_burden.tsv",
        exact_comparison=f"{OUTDIR}/analysis/exact_deletion_comparison.tsv",
        affected_comparison=f"{OUTDIR}/analysis/affected_feature_comparison.tsv",
        impact_comparison=f"{OUTDIR}/analysis/feature_impact_class_comparison.tsv",
        known_sequence_summary=f"{OUTDIR}/analysis/known_sequence_search_summary.tsv",
        known_sequence_hits=f"{OUTDIR}/analysis/known_sequence_search_hits.tsv",
        exact_mtpm=f"{OUTDIR}/matrices/exact_deletion_support_per_million_mt_reads.tsv",
        affected_mtpm=f"{OUTDIR}/matrices/affected_feature_support_per_million_mt_reads.tsv",
        burden_plot=f"{OUTDIR}/plots/deletion_burden_by_sample.pdf",
        rainfall_left=f"{OUTDIR}/plots/deletion_rainfall_left_breakpoint.pdf",
        rainfall_right=f"{OUTDIR}/plots/deletion_rainfall_right_breakpoint.pdf",
        rainfall_midpoint=f"{OUTDIR}/plots/deletion_rainfall_midpoint.pdf",
        breakpoint_pair_map=f"{OUTDIR}/plots/breakpoint_pair_support_map.pdf",
        affected_support=f"{OUTDIR}/plots/affected_feature_support.pdf",
        recurrence=f"{OUTDIR}/plots/exact_deletion_recurrence.pdf",
    output:
        complete=f"{DELIVERABLES_DIR}/DELIVERABLES_COMPLETE.txt",
    params:
        results_dir=OUTDIR,
        outdir=DELIVERABLES_DIR,
        dataset=DATASET,
    conda:
        "envs/mitochondrial-deletions.yaml"
    shell:
        "python scripts/make_deliverables.py --results-dir {params.results_dir} --dataset {params.dataset} "
        "--config {input.config} --defaults config/defaults.yaml "
        "--output-dir {params.outdir} --complete {output.complete}"
