# Mitochondrial Deletion Discovery

Configurable Snakemake workflow for comparing sequencing datasets for mitochondrial deletion evidence. The workflow supports RNA and DNA read inputs through dataset configuration; sample names, accessions, and biological conclusions are not hard-coded into rules or scripts.

Current dataset configs:

- rat aging muscle GPA dataset from NCBI BioProject `PRJNA793055`;
- local human common mtDNA deletion dataset with HDFn, KSS-95, and KSS-96 FASTQs;
- local human nanopore mtDNA dataset using an uncompressed FASTQ staged by the workflow;
- local matched human bulk sequencing dataset.

## Conceptual Model

The main report is focused on mitochondrial circular-remap deletion calls. The first-pass genome alignment is used as read selection and provenance, not as a reported biological result stream. The default selection mode is competitive whole-genome assignment: reads are aligned with nuclear chromosomes and mtDNA present together, then reads whose best/selected evidence is mitochondrial are passed to mitochondrial remapping. Nuclear-only unmapped-read selection and the older full-genome mitochondrial-evidence scanner remain available as configuration-driven sensitivity/reproducibility modes.

**Mitochondrial circular-remap results.** Selected reads are remapped with minimap2 to mitochondrial-only references and converted back to the original mtDNA coordinate system. Query order and alignment strand define a directed retained adjacency `L -> R`; the inferred deleted interval is the forward circular arc from retained base `L` to retained base `R`. The reciprocal `R -> L` adjacency is a different deletion model and is not collapsed by choosing the shorter arc. By default these results are normalized per million usable reads after read preparation. The denominator can be changed with `analysis.normalization_denominator`; supported values are `total_usable_reads` and `mt_evidence_reads`.

**Mitochondrial-evidence reads** are the reads retained after first-pass genome assignment because their best or selected alignment evidence is mitochondrial. They are the reads written to the remap-input FASTQs. They are the input to circular remapping and can optionally be used as the per-million denominator, but they are not the default denominator. This is different from the local breakpoint reference-support denominator described below.

**Supplementary configured sequence searches.** Dataset configs can define literal breakpoint-spanning sequences to search in the retained remap-input FASTQs. This is useful for sanity-checking named deletions such as the human common mtDNA deletion, but it only detects the configured motifs and is not a replacement for the remapped split-read caller.

Main stages:

1. Stage or download FASTQs.
2. Optionally trim with `fastp`.
3. Use a first-pass genome alignment to select reads whose best/selected evidence is mitochondrial.
4. Optionally run configured literal sequence searches over retained remap-input FASTQs.
5. Remap retained reads with minimap2 to normal and rotated mitochondrial references.
6. Infer deletion-like events from minimap2 split/supplementary alignments.
7. Convert rotated coordinates back to original mtDNA coordinates.
8. Canonicalize circular breakpoint pairs and deduplicate support across rotations.
9. Annotate exact deletions by affected mitochondrial features and configured deletion target matches.
10. Build exact-deletion and affected-feature matrices, statistics, plots, and report.

## Result Levels

The report is organized around the following result levels.

**Exact deletions** are directed coordinate-level inferred deletion models. They have alignment-directed left and right retained-flanking breakpoints, deleted size, wrapping status, complement diagnostics, direction and rotation status, support, normalized support, sample/group labels, and optional configured target labels such as the human common mtDNA deletion.

For each exact deletion, the workflow also estimates local reference-spanning support at the two breakpoints. This asks a narrow question: in the same mitochondrial remap stream, how many primary alignments span the left and right breakpoint neighborhoods without requiring a deletion split? The workflow counts local spanning depth in the normal and rotated mitochondrial remaps and uses the larger count for each breakpoint, which avoids summing the same evidence twice across rotations. The report gives the left and right reference-spanning counts, the smaller of those two counts, and a local split-support fraction:

`split-supporting reads / (split-supporting reads + minimum local reference-spanning reads)`

This is a local alignment-support metric, not the denominator used for the main per-million plots. For RNA data, it should not be interpreted as mtDNA heteroplasmy. For DNA data, it is a local breakpoint-support summary rather than a complete heteroplasmy model unless the dataset and coverage assumptions justify that interpretation. It is most interpretable when reads are long enough to span the configured breakpoint windows and when both breakpoint neighborhoods have coverage. The denominator is calculated for remap-called exact deletions because the numerator and reference-spanning counts come from the same remapped read set. Configured sequence searches are kept as supplementary literal motif checks; their counts are not converted into this denominator because a motif hit alone does not define the comparable non-deletion spanning-read population.

**Affected-feature categories** are biology-level events. For each deletion, the workflow determines which annotated mitochondrial genes or features overlap the deleted interval. Feature names come from the reference annotation, are sorted by genomic order, and are joined with `+`, for example `MT-ATP6+MT-CO3+MT-ND3`. This makes group comparisons more stable when breakpoints vary slightly but affect the same genes.

The annotation step reduces raw GTF rows to one biological feature per gene/name before assigning affected-feature labels. This avoids counting separate gene, transcript, exon, and CDS records for the same biological feature. Dataset configs can add noncoding mitochondrial regions under `analysis.mt_regions`, such as control-region/D-loop intervals, direct-repeat windows, origins, or other coordinate intervals that are biologically useful but absent from the GTF.

Dataset configs can also add `annotations.feature_aliases` when a reference annotation uses accession-like or otherwise unhelpful mitochondrial feature names. Aliases affect report labels, affected-feature categories, plots, and comparison tables while keeping the raw reference annotation available in the generated annotation files and resolved config.

The report also includes collapsed impact classes such as single-feature, two-feature, multi-feature mixed, tRNA/rRNA-involved, D-loop-involved, and intergenic. Exact labels are retained in tables even when plots collapse or limit categories for readability. Long affected-feature plots show the top categories plus an `Other categories` bin; full labels remain available in the tables.

## Circular Mitochondrial Remapping

Mitochondrial DNA is circular, but short-read aligners use linear reference sequences. A deletion-like junction can therefore be represented differently depending on where the artificial linear reference begins.

The workflow handles this by using two non-duplicated mitochondrial references:

- the normal mtDNA sequence;
- a second mtDNA sequence rotated by half the mitochondrial genome.

The primary caller does **not** use a duplicated reference such as `mtDNA + mtDNA`, because the repeated sequence creates avoidable secondary placement. A doubled or padded reference can be useful as an optional diagnostic for artificial linear-boundary behavior when it is paired with an explicit central-copy coordinate policy.

After minimap2 remapping:

- coordinates from the rotated reference are converted back to the original mtDNA coordinate system without reversing the directed adjacency;
- query order and strand define which retained flank occurs before and after the junction;
- the deleted interval is the forward circular arc from the directed left breakpoint to the directed right breakpoint, excluding both retained breakpoint bases;
- reciprocal `L -> R` and `R -> L` junctions remain distinct exact deletions;
- the same read supporting the same directed deletion in both rotations is counted once in matrices, burden summaries, statistics, and reports;
- same-read reciprocal conflicts across rotations are retained as ambiguous evidence and excluded from primary summaries by default;
- final exact-deletion IDs are generated from directed coordinates and deleted size.

Example:

- A read supporting `8472 -> 13448` implies a deleted interval of about 5 kb between those retained flanks.
- A read supporting `13448 -> 8472` implies the complementary origin-spanning interval of about 11.6 kb.
- The two models share an unordered breakpoint pair but are not biologically interchangeable. If alignments cannot distinguish them reliably, the evidence is marked ambiguous rather than resolved by interval length.

Circular handling is core correctness logic, not a plotting trick. Unit tests cover strand-normalized directed junctions, coordinate conversion, separate reciprocal models, long and origin-spanning arcs, wrapping interval annotation, direction conflicts, and rotated-reference deduplication.

The deleted bases are absent from a junction-spanning read. Their identity is inferred from the new adjacency of the retained flanks. For example, a read containing reference sequence ending at retained base `L` followed immediately by sequence beginning at retained base `R` supports the `L|R` junction and the forward `L -> R` deleted arc. This is sequence evidence for an inferred deletion model, not by itself proof that the source molecule was a viable mtDNA genome.

## Assay Identity And Interpretation

Dataset configuration should state `dataset.read_technology`, `dataset.molecule_type`, and `dataset.library_strategy`. Reports use these fields to select applicable assumptions; they do not infer DNA versus RNA from sample names or mapper presets. Use `unknown` when the assay has not been confirmed.

For nanopore data, long anchors can support direct junction inspection, but base errors, homopolymers, supplementary or alternative placements, ligation artifacts, chimeric reads, and concatemers can affect calls. For Illumina data, short split anchors can be difficult to place uniquely around repeats and NUMTs; the current workflow does not call a deletion from mate distance alone.

For RNA-derived data, transcript processing, reverse transcription, and template switching can create deletion-like split alignments, and read abundance does not measure mtDNA heteroplasmy. For DNA-derived data, split evidence is closer to genome-molecule evidence but can reflect NUMTs, PCR or ligation chimeras, mapping ambiguity, and sampling. In either case, coordinate evidence should not be described as a confirmed biological deletion without appropriate context and validation.

## RNA-Specific Transcript Artifact Handling

Mitochondrial RNA-seq split alignments can reflect RNA processing or alignment artifacts rather than mtDNA deletions. For RNA datasets, the workflow treats minimap2 split alignments as evidence that may support a deletion model, not as automatic biological deletions.

Expected adjacent transcript junctions are configured by feature name, for example
`MT-ATP8` to `MT-ATP6` or `MT-ND4L` to `MT-ND4`. These labels are applied after
feature aliases are resolved, so species-specific raw annotation names can be
converted to readable mitochondrial feature names first. With the default
`junctions.exclude_expected_transcript_junctions: true`, transcript-compatible
split reads are excluded from deletion burden, exact-deletion, affected-feature,
and group-comparison summaries. Their counts are reported in QC tables so
the filtering is visible.

Current filters and annotations include:

- minimum anchor length on both sides of a split;
- minimum and maximum deletion size;
- minimum split-read support after clustering;
- mapping-quality and supplementary/secondary-alignment settings;
- query overlap and query gap limits between split segments;
- expected transcript-junction annotation from `annotations.expected_adjacent_transcripts`;
- optional exclusion of configured transcript-compatible junctions from deletion summaries with `junctions.exclude_expected_transcript_junctions`;
- affected interval features, breakpoint-nearest features, size class, and configured deletion target matching.

Important configurable settings:

- `junctions.min_anchor_length`
- `junctions.arc_assignment`: defaults to `alignment_directed`, where stored split-alignment query order determines the directed deleted interval.
- `junctions.alignment_pairing_mode`: defaults to all compatible segment pairs within one physical read or paired-end mate. Paired-end mates are never merged into the same alignment chain.
- `junctions.ambiguous_direction_policy`: defaults to excluding reciprocal same-read conflicts from primary summaries while retaining them for QC.
- `junctions.min_deletion_size`
- `junctions.max_deletion_size`
- `junctions.breakpoint_slop_bp`
- `junctions.min_split_read_support`
- `junctions.exclude_expected_transcript_junctions`
- `mt_realign.minimap2_min_mapq`
- `mt_realign.minimap2_index_extra`
- `mt_realign.minimap2_extra`
- `mt_realign.min_segment_aligned_fraction`
- `mt_realign.max_soft_clip_fraction`
- `mt_realign.max_query_overlap_bp`
- `mt_realign.max_query_gap_bp`
- `analysis.mt_regions`
- `analysis.known_deletions`
- `analysis.known_sequence_searches`: supplementary literal motif checks. If a search also provides `left_breakpoint`/`right_breakpoint` or has coordinate-like text such as `mtDNA_8471_13449`, the workflow uses that as an additional configured deletion target for labeling nearby remap calls.
- `analysis.effect_size_pseudocount_per_million`
- `annotations.feature_aliases`

Future improvements to consider:

- richer mitochondrial transcript-processing site models;
- explicit artifact tiers such as high-confidence deletion-like, processing-compatible, ambiguous, and likely artifact;
- read-level evidence summaries for top calls;
- optional candidate-junction template realignment;
- optional confirmation with BWA-MEM2 or another short-read split aligner.

`mt_realign.minimap2_index_extra` is for index-sensitive minimap2 settings such
as `-k` and `-w`. Put those settings there so the `.mmi` index is built with the
intended seed parameters. Use `mt_realign.minimap2_extra` for mapping-time
options that do not affect index construction.

## Mapper Choices And Read Types

The current default workflow uses competitive whole-genome first-pass assignment followed by minimap2 mitochondrial remapping. In default short-read mode, STAR maps reads to the full genome including mtDNA and streams the unsorted alignment output directly into the mitochondrial-read selector. The selector uses read-name collation, not a full coordinate BAM plus name-sort, so the default path keeps the small remap-input FASTQs and STAR logs rather than large whole-genome BAMs. For long-read mode, minimap2 can map reads to the full genome including mtDNA and retain reads whose best evidence is mitochondrial. STAR or minimap2 first-pass output is read-selection/provenance, not a separate reported biological deletion stream.

For short-read first-pass assignment, HISAT2 is worth evaluating because it is fast and may be sufficient if the goal is to assign reads competitively between nuclear and mitochondrial references. It is not the default because this repository currently has tested workflow rules for STAR and minimap2 first-pass selection, while HISAT2 would add another index/output convention to validate.

For long-read inputs, minimap2 is the natural first-pass mapper. Set `mapping.first_pass_aligner: minimap2` and choose a preset appropriate to the chemistry and assay. For RNA, use a splice-aware transcriptomic preset when nuclear mapping needs to recognize introns. For DNA, use a genomic long-read preset. The mitochondrial remap step can also use minimap2 for long reads, but should use a long-read preset and less short-read-specific split-segment assumptions.

Phase 2 currently uses direct split/supplementary read alignments to support inferred deletion models. It does not infer deletions from mate-pair distance alone. For paired short reads, mate information can help QC or future confidence scoring, but the reported deletion calls require breakpoint-spanning split evidence after circular coordinate handling.

The first-pass selection mode is controlled by `mapping.first_pass_read_selection`:

- `whole_genome_mt_best` maps against the full genome including mtDNA and passes reads with mitochondrial best/selected evidence to mitochondrial remapping. This is the default. With `mapping.first_pass_aligner: star`, this path streams STAR output through read-name collation into the selector and does not create full-genome BAMs.
- `nuclear_unmapped_reads` maps against a nuclear-only reference and passes unmapped reads to mitochondrial remapping. This is retained for strict depletion-style sensitivity checks, but it can discard real mitochondrial reads with NUMT-like nuclear alignments.
- `mt_evidence_reads` is an alternative mode that maps against the full genome and scans the BAM/chimeric output for mitochondrial evidence before remapping.

The first-pass aligner is controlled by `mapping.first_pass_aligner`. In `whole_genome_mt_best` mode, the implemented choices are `star` for short-read RNA-oriented workflows and `minimap2` for long-read, DNA, or mapper-sensitivity experiments. The `mapping.keep_ambiguous_mt_nuclear_reads` setting controls whether reads whose primary evidence is mitochondrial but whose MAPQ or secondary alignments suggest nuclear ambiguity are retained for remapping.

For nanopore/ONT FASTQ inputs, use a dataset config with local FASTQ paths, `qc.trim_reads: false`, and `mapping.first_pass_aligner: minimap2`. Local uncompressed `.fastq` files are accepted; the staging rule writes a gzipped copy under the dataset results directory. A typical ONT configuration is:

```yaml
qc:
  run_fastqc: false
  trim_reads: false

mapping:
  first_pass_read_selection: whole_genome_mt_best
  first_pass_aligner: minimap2
  first_pass_minimap2_preset: splice
  keep_ambiguous_mt_nuclear_reads: true

mt_realign:
  minimap2_preset: map-ont
```

## Offline Reference Files

The workflow can run without internet access when references and sample FASTQs are provided locally. In a dataset config, set local paths instead of URLs:

```yaml
references:
  human:
    genome_path: /path/to/genome.fa.gz
    annotation_path: /path/to/annotation.gtf.gz
    mt_reference_name: Revised Cambridge Reference Sequence (rCRS)
    mt_reference_accession: NC_012920.1
    mt_contig_names: [MT, chrM, MTDNA]
    mt_length: 16569
```

Local `.gz` files are decompressed into the workflow reference directory. URL-based downloads remain supported for testing and for users who want the workflow to fetch references.
The optional `mt_reference_name` and `mt_reference_accession` fields are shown in the report so collaborators can confirm the mitochondrial coordinate standard used for deletion coordinates.

Local sample FASTQs can be provided in the sample table with `fastq_1` and optional `fastq_2` columns. Both `.fastq.gz` and uncompressed `.fastq` are supported for local inputs.

Resolved SRA metadata is cached under `metadata/cache/`. Once a cache file exists, metadata resolution uses it before trying NCBI again, so reruns after configuration changes do not require live NCBI metadata access.

## Report Outputs

The main report answers:

- do groups differ in total deletion burden?
- do groups differ in the number of distinct deletions?
- in two-factor datasets, do age, treatment, or their interaction explain sample-level deletion outcomes?
- do groups differ in deletion size distribution?
- do groups differ in specific exact deletions?
- do groups differ in which genes/features are affected?
- do samples cluster by exact deletion profile?
- do samples cluster by affected-feature profile?
- are continuous metadata variables associated with deletion burden or diversity?

Main plots include:

- total deletion burden by sample and group;
- distinct exact deletions by sample and group;
- age-by-treatment interaction plots when `age` and `treatment` metadata are present;
- deletion size distributions, including log-y and size-restricted views;
- deletion rainfall plots shown as full-size per-group figures by left breakpoint, right breakpoint, and circular deleted-interval midpoint;
- breakpoint-pair support maps showing which deletion starts pair with which deletion ends;
- group-split pooled breakpoint support-density plots showing where deletion endpoints accumulate, with binned support split into left and right breakpoint bars behind a circular-smoothed total-support curve;
- affected-feature normalized support and within-group proportions;
- collapsed feature-impact classes;
- per-gene affected burden;
- group-colored exact deletion recurrence;
- PCA and Bray-Curtis MDS for exact deletions and affected-feature categories, without static sample labels or centroids.

Plot group order and colors are chosen once per report and reused across plots. For two-factor designs with `age` and `treatment`, groups are ordered by age and then treatment, with control-like groups first within each age. Control-like groups use subdued colors and treatment groups use red-family colors when possible.

The deletion rainfall plots, breakpoint-pair support map, and pooled breakpoint support-density plots are visual displays, not separate calling steps. They use the same display rule. By default `plots.rainfall_min_support_per_million: 0.0`, so low-abundance datasets are not hidden by a normalized-support cutoff. `plots.rainfall_max_points_per_group: 300` limits each group panel to the highest-support displayed exact deletions so dense datasets remain readable. These settings affect only the plots; the exact-deletion tables, matrices, comparisons, and read-list links retain the full analyzed call set subject to their own table-display filters.

Large result tables in the HTML report are searchable, sortable, and paged. Small tables are shown directly without search controls.

Heatmaps are intentionally not part of the main report.

## Repository Layout

- `Snakefile` - workflow rules.
- `config/defaults.yaml` - general workflow defaults.
- `config/datasets/*.yaml` - dataset-specific configuration.
- `envs/mitochondrial-deletions.yaml` - conda environment used by Snakemake rules.
- `scripts/` - Python workflow scripts.
- `docs/workflow_methods_and_assumptions.md` - detailed workflow stages, coordinate semantics, assumptions, and assay-specific interpretation.
- `docs/directed_circular_deletion_workflow.md` - normative circular-arc, alignment-chain, clustering, reporting, and validation requirements.
- `tests/` - focused unit tests for metadata, directed circular coordinates, deletion evidence, reporting helpers, and parsing.
- `planning/` - planning notes and reference analysis material.
- `AGENTS.md` - local instructions for Codex or other coding agents.

Generated files are intentionally ignored by git. See `.gitignore`.

## Install

Use the listed channel order with strict channel priority so compiled packages
come from compatible conda-forge/bioconda builds:

```bash
conda config --set channel_priority strict
conda env create -f envs/mitochondrial-deletions.yaml
conda activate mitochondrial-deletions
```

The environment is expected to provide the workflow, Python analysis stack, and
external command-line tools. Verify a new install with:

```bash
python -m pytest -q
snakemake --version
minimap2 --version
samtools --version | head -1
fasterq-dump --version | head -2
fastp --version
seqkit version
```

Then run a workflow dry-run before downloading data or building results:

```bash
snakemake --configfile config/datasets/rat_aging_muscle.yaml --dry-run --cores 1
```

Snakemake can also create per-rule environments under `.snakemake/conda/` when `--use-conda` is used.

## Run The Rat GPA Dataset

Dry-run first:

```bash
snakemake --use-conda --cores 8 \
  --resources download=2 \
  --configfile config/datasets/rat_aging_muscle.yaml \
  --dry-run
```

Then run:

```bash
snakemake --use-conda --cores 8 \
  --resources download=2 \
  --configfile config/datasets/rat_aging_muscle.yaml
```

Main report:

```text
results/rat_aging_muscle/rat_aging_muscle_deliverables/index.html
```

## Run The Human Common Deletion Dataset

The local common-deletion FASTQs are expected under:

```text
BCLConvert_07_25_2024_15_02_59Z-752430691/
```

The sample table is:

```text
metadata/human_common_deletion.samples.tsv
```

These can be raw BCLConvert FASTQs. The normal workflow stages local FASTQs, runs FastQC, trims with `fastp`, maps trimmed reads, then performs minimap2 mitochondrial remapping.

Dry-run first:

```bash
snakemake --use-conda --cores 8 \
  --configfile config/datasets/human_common_deletion.yaml \
  --dry-run
```

Then run:

```bash
snakemake --use-conda --cores 8 \
  --configfile config/datasets/human_common_deletion.yaml
```

Main report:

```text
results/human_common_deletion/human_common_deletion_deliverables/index.html
```

## Run The Human Nanopore Dataset

The included nanopore example expects this local single-end FASTQ:

```text
SUP_ONT_rCS.fastq
```

The sample table and config are:

```text
metadata/human_nanopore.samples.tsv
config/datasets/human_nanopore.yaml
```

This config disables FastQC and fastp trimming, stages the local uncompressed FASTQ as gzip, uses minimap2 for competitive whole-genome mitochondrial read selection, then remaps retained reads to normal and rotated mitochondrial references with `map-ont`.

Dry-run first:

```bash
snakemake results/human_nanopore/human_nanopore_deliverables/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/human_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --dry-run
```

Then run:

```bash
snakemake results/human_nanopore/human_nanopore_deliverables/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/human_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete
```

Main report:

```text
results/human_nanopore/human_nanopore_deliverables/index.html
```

## Advanced: Rerun From Existing Trimmed FASTQs

Use this only when downloads and trimming are already complete and you do not want Snakemake to backtrack into raw FASTQ staging or trimming rules.

Dry-run first:

```bash
snakemake results/rat_aging_muscle/rat_aging_muscle_deliverables/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/rat_aging_muscle.yaml \
  --config workflow_start_from=trimmed \
  --rerun-triggers mtime \
  --dry-run
```

Check that the dry-run does not include `prepare_reads`, `fastqc_raw`, or `trim_reads`. Then run the same command without `--dry-run`.

## Key Outputs

Final deliverables:

- `results/<dataset>/<dataset>_deliverables/index.html`
- `results/<dataset>/<dataset>_deliverables/tables/`
- `results/<dataset>/<dataset>_deliverables/plots/`
- `results/<dataset>/<dataset>_deliverables/matrices/`
- `results/<dataset>/<dataset>_deliverables/config/resolved_config.yaml`

Important machine-readable outputs:

- `junctions/junction_clusters.tsv` - alignment-directed exact deletions with annotation, direction status, complement diagnostics, rotation agreement, and schema version.
- `junctions/all_samples.filtered_junction_reads.tsv` - filtered remap read-level rows before matrix-level normal/rotated deduplication.
- `junctions/ambiguous_direction_reads.tsv` - reciprocal-direction conflicts retained for audit and excluded from primary summaries by default.
- `<dataset>_deliverables/tables/run_methods.tsv` - machine-readable resolved settings used by the run.
- `<dataset>_deliverables/tables/data_dictionary.tsv` - delivered table/column inventory and definitions.
- `analysis/breakpoint_reference_support.tsv` - local reference-spanning read counts at exact-deletion breakpoints and the resulting local split-support fraction.
- `matrices/exact_deletion_raw_counts.tsv`
- `matrices/exact_deletion_support_per_million_mt_reads.tsv`
- `matrices/affected_feature_raw_counts.tsv`
- `matrices/affected_feature_support_per_million_mt_reads.tsv`
- `analysis/deletion_burden.tsv`
- `analysis/exact_deletion_comparison.tsv`
- `analysis/affected_feature_comparison.tsv`
- `analysis/feature_impact_class_comparison.tsv`
- `analysis/deletion_size_distribution_tests.tsv`
- `analysis/deletion_size_bin_summary.tsv`
- `analysis/factorial_model_summary.tsv` when a two-factor age-by-treatment design is available
- `analysis/per_gene_affected_burden.tsv`
- `analysis/qc_summary.tsv`
- `analysis/known_sequence_search_summary.tsv` when configured sequence searches are present
- `analysis/known_sequence_search_hits.tsv` when configured sequence searches are present

The `*_per_million_mt_reads.tsv` filenames are stable output names and do not identify the normalization denominator by themselves. Check the `normalization_denominator` and `normalization_reads` columns in `analysis/deletion_burden.tsv` and the report method table to determine whether the run used total usable reads or retained mitochondrial-evidence reads as the denominator.

For usability, the HTML report embeds a filtered exact-deletions table by default: exact deletions with at least 50 supporting reads are shown, and configured deletion-target matches are always retained. The complete unfiltered exact-deletion table is delivered as `tables/exact_deletions.tsv`. Adjust this display-only behavior with `report.exact_deletion_table` in the configuration.

When read-level evidence is available, the HTML report links read-count cells to sidecar TSVs in `read_lists/`. Exact-deletion support counts link to the reads supporting that deletion, configured deletion-target remap counts link to the reads supporting all nearby remap calls assigned to that target, and configured sequence-search counts link to reads containing the configured literal motif.

The dataset name is included in the deliverables folder name so copied report folders from different datasets do not overwrite each other.

## Testing

```bash
conda activate mitochondrial-deletions
python -m pytest -q
```

## About `fasterq.tmp.*` Folders

`fasterq-dump` can create temporary folders named like `fasterq.tmp.<host>.<pid>` if it is not given an explicit temp path. The workflow wrapper passes an explicit temp directory inside each accession staging directory under `results/<dataset>/fastq/.<accession>.fasterq/tmp`. Old top-level `fasterq.tmp.*` folders are leftover scratch folders and should not be committed.

## Notes For Making This A Git Repository

Recommended to commit:

- `Snakefile`, `scripts/`, `tests/`;
- `config/`;
- `envs/`;
- `README.md`, `AGENTS.md`, and intentional planning docs.

Recommended not to commit:

- `.snakemake/`;
- `results/`;
- `resources/` if it contains downloaded references or indexes;
- `BCLConvert_*/`;
- FASTQs, BAMs, BAI files, indexes, and SRA caches;
- top-level `fasterq.tmp.*` folders;
- local agent/session folders such as `.agents/` and `.codex/`.
