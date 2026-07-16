# Mitochondrial Deletion Discovery

Configurable Snakemake workflow for comparing sequencing datasets for mitochondrial deletion evidence. The workflow supports RNA and DNA read inputs through dataset configuration; sample names, accessions, and biological conclusions are not hard-coded into rules or scripts.

Dataset configs included in the repository:

- rat aging muscle GPA dataset from NCBI BioProject `PRJNA793055`;
- local human common mtDNA deletion dataset with HDFn, KSS-95, and KSS-96 FASTQs;
- local human nanopore mtDNA dataset using an uncompressed FASTQ staged by the workflow;
- local matched human bulk sequencing dataset.

## Conceptual Model

The workflow reports coordinate-focused mitochondrial deletion evidence after circular canonicalization, filtering, physical-observation deduplication, annotation, and quality profiling. The default first-pass selection mode is competitive whole-genome assignment: reads are aligned with nuclear chromosomes and mtDNA present together, then reads whose best/selected evidence is mitochondrial are passed to mitochondrial remapping. `nuclear_unmapped_reads` and `mt_evidence_reads` are configuration-driven alternatives for sensitivity and reproducibility checks.

**Mitochondrial circular-remap results.** Selected reads are remapped with minimap2 to mitochondrial-only references and converted back to the original mtDNA coordinate system. Query order and alignment strand define a directed retained adjacency `L -> R`; the inferred deleted interval is the forward circular arc from retained base `L` to retained base `R`. The reciprocal `R -> L` adjacency is a different deletion model and is not collapsed by choosing the shorter arc. By default these results are normalized per million usable reads after read preparation. The denominator can be changed with `analysis.normalization_denominator`; supported values are `total_usable_reads` and `mt_evidence_reads`.

**Mitochondrial-evidence reads** are the reads retained after first-pass genome assignment because their best or selected alignment evidence is mitochondrial. They are the reads written to the remap-input FASTQs. They are the input to circular remapping and can optionally be used as the per-million denominator, but they are not the default denominator. This is different from the local breakpoint reference-support denominator described below.

**Short-read RNA dual-caller evidence.** When `quality.short_read_rna_dual_caller.enabled: true`, mitochondrial-to-mitochondrial records from STAR `Chimeric.out.junction` are evaluated alongside minimap2 remap evidence. STAR supplies an additional short-read junction-detection route; it does not replace competitive whole-genome read selection or circular mitochondrial remapping. STAR and minimap2 rows are converted to the same directed circular coordinate model, and the same physical observation supporting the same event in both callers is counted once. STAR-Fusion is not used.

**Supplementary configured sequence searches.** Dataset configs can define literal breakpoint-spanning sequences to search in the retained remap-input FASTQs. This is useful for sanity-checking named deletions such as the human common mtDNA deletion, but it only detects the configured motifs and is not a replacement for the remapped split-read caller.

Main stages:

1. Stage or download FASTQs.
2. Optionally trim with `fastp`.
3. Use a first-pass genome alignment to select reads whose best/selected evidence is mitochondrial.
4. Optionally run configured literal sequence searches over retained remap-input FASTQs.
5. Remap retained reads with minimap2 to normal and rotated mitochondrial references.
6. Infer deletion-like events from minimap2 split/supplementary alignments and, when configured for short-read RNA, STAR chimeric alignments.
7. Convert all evidence to directed original-reference mtDNA coordinates.
8. Canonicalize circular breakpoints and deduplicate physical observations across rotations and callers.
9. Annotate exact deletions, assign evidence-quality tiers, and record caller concordance and quality flags.
10. Build separate stringent, standard, and exploratory matrices, statistics, plots, and reports from the same canonical event set.

## Result Levels

The report is organized around the following result levels.

**Exact deletions** are directed coordinate-level inferred deletion models. They have alignment-directed left and right retained-flanking breakpoints, deleted size, wrapping status, complement diagnostics, direction and rotation status, support, normalized support, sample/group labels, reference-specific major/minor replication-arc context, and optional configured target labels such as the human common mtDNA deletion.

For each exact deletion with minimap2 remap evidence, the workflow also estimates local reference-spanning support at the two breakpoints. This asks a narrow question: in the same mitochondrial remap stream, how many primary alignments span the left and right breakpoint neighborhoods without requiring a deletion split? The workflow counts local spanning depth in the normal and rotated mitochondrial remaps and uses the larger count for each breakpoint, which avoids summing the same evidence twice across rotations. STAR-only calls are marked unavailable because a STAR chimeric numerator and minimap2 remap denominator would not be comparable. The report gives the left and right reference-spanning counts, the smaller of those two counts, and a local split-support fraction:

`split-supporting reads / (split-supporting reads + minimum local reference-spanning reads)`

This is a local alignment-support metric, not the denominator used for the main per-million plots. For RNA data, it should not be interpreted as mtDNA heteroplasmy. For DNA data, it is a local breakpoint-support summary rather than a complete heteroplasmy model unless the dataset and coverage assumptions justify that interpretation. It is most interpretable when reads are long enough to span the configured breakpoint windows and when both breakpoint neighborhoods have coverage. The denominator is calculated for remap-called exact deletions because the numerator and reference-spanning counts come from the same remapped read set. Configured sequence searches are kept as supplementary literal motif checks; their counts are not converted into this denominator because a motif hit alone does not define the comparable non-deletion spanning-read population.

**Affected-feature categories** are deterministic interval annotations, not exact-deletion identities. For each deletion, the workflow determines which annotated mitochondrial genes or features overlap the deleted interval. Feature names come from the reference annotation, are sorted by genomic order, and are joined with `+`, for example `MT-ATP6+MT-CO3+MT-ND3`. This makes group comparisons more stable when breakpoints vary slightly but affect the same genes.

The annotation step reduces raw GTF rows to one biological feature per gene/name before assigning affected-feature labels. This avoids counting separate gene, transcript, exon, and CDS records for the same biological feature. Dataset configs can add noncoding mitochondrial regions under `analysis.mt_regions`, such as control-region/D-loop intervals, direct-repeat windows, origins, or other coordinate intervals that are biologically useful but absent from the GTF.

Reference-specific `references.<species>.replication_arcs` entries define the mitochondrial **minor arc** and **major arc** coordinate convention used for reporting. Every exact deletion records `replication_arc_context` as `minor_arc_only`, `major_arc_only`, or `major_and_minor_arcs`, together with the number of deleted reference bases overlapping each arc. These replication-arc annotations are separate from `analysis.mt_regions`: they do not enter affected-feature labels, and they do not choose the alignment-directed deleted circular interval. The configured report table states the boundary coordinates and whether an arc wraps the artificial coordinate origin.

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

The phrase *alignment-directed deleted arc* refers to the circular interval inferred from one read-supported `L -> R` adjacency. *Major arc* and *minor arc* refer instead to fixed, reference-specific regions bounded by configured replication-origin landmarks. A called deletion can lie only in the major arc, only in the minor arc, or overlap both; the fixed arc labels never resolve a reciprocal-direction ambiguity.

Example:

- A read supporting `8472 -> 13448` implies a deleted interval of about 5 kb between those retained flanks.
- A read supporting `13448 -> 8472` implies the complementary origin-spanning interval of about 11.6 kb.
- The two models share an unordered breakpoint pair but are not biologically interchangeable. If alignments cannot distinguish them reliably, the evidence is marked ambiguous rather than resolved by interval length.

Circular handling is core correctness logic, not a plotting trick. Unit tests cover strand-normalized directed junctions, coordinate conversion, separate reciprocal models, long and origin-spanning arcs, wrapping interval annotation, direction conflicts, and rotated-reference deduplication.

The deleted bases are absent from a junction-spanning read. Their identity is inferred from the new adjacency of the retained flanks. For example, a read containing reference sequence ending at retained base `L` followed immediately by sequence beginning at retained base `R` supports the `L|R` junction and the forward `L -> R` deleted arc. This is sequence evidence for an inferred deletion model, not by itself proof that the source molecule was a viable mtDNA genome.

## Assay Identity And Interpretation

Dataset configuration should state `dataset.read_technology`, `dataset.molecule_type`, `dataset.assay_type`, and `dataset.library_strategy`. Reports use these fields to select applicable assumptions; they do not infer DNA versus RNA or bulk versus single-cell design from sample names or mapper presets. Use `unknown` when an assay property has not been confirmed.

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

Filters and annotations include:

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
- `quality.enabled`: enables canonical quality evidence tables and profile reports.
- `quality.minimum_supported_observations` and `quality.minimum_strong_observations`: define the observation-count components of evidence tiers.
- `quality.primary_report_profile`: identifies the profile presented as the primary interpretation in the report index.
- `quality.report_profiles.<profile>.include_tiers`: defines cumulative profile membership without changing stable exact-deletion IDs.
- `quality.short_read_rna_dual_caller.enabled`: enables STAR chimeric evidence in addition to minimap2 remap evidence for explicitly configured short-read RNA datasets.
- `quality.short_read_rna_dual_caller.star_min_anchor_length`, `star_max_query_overlap_bp`, `star_max_query_gap_bp`, `star_require_gene_anchors`, and `star_exclude_same_gene`: define STAR-specific candidate filters.

The quality layer records caller provenance, physical-observation support, within-sample replication, cross-caller corroboration, alignment geometry, breakpoint dispersion, rotation support, transcript compatibility, and applicable quality flags. These fields remain visible in the canonical tables even when a report profile excludes the event.

`mt_realign.minimap2_index_extra` is for index-sensitive minimap2 settings such
as `-k` and `-w`. Put those settings there so the `.mmi` index is built with the
intended seed parameters. Use `mt_realign.minimap2_extra` for mapping-time
options that do not affect index construction. Minimap2 indexes are named from
their preset and index-sensitive settings, so datasets with incompatible seed
profiles cannot silently reuse the same whole-genome or mitochondrial index.

## Mapper Choices And Read Types

The default workflow uses competitive whole-genome first-pass assignment followed by minimap2 mitochondrial remapping. In short-read mode, STAR maps reads to the full genome including mtDNA and streams the unsorted alignment output directly into the mitochondrial-read selector. The selector uses read-name collation, not a full coordinate BAM plus name-sort, so the path keeps the small remap-input FASTQs and STAR logs rather than large whole-genome BAMs. When the short-read RNA dual-caller quality stream is enabled, STAR chimeric junction records also contribute candidate evidence after STAR-specific filtering and canonicalization. For long-read mode, minimap2 maps reads to the full genome including mtDNA and retains reads whose best evidence is mitochondrial.

For long-read inputs, minimap2 is the natural first-pass mapper. Set `mapping.first_pass_aligner: minimap2` and choose a preset appropriate to the chemistry and assay. For RNA, use a splice-aware transcriptomic preset when nuclear mapping needs to recognize introns. For DNA, use a genomic long-read preset. The mitochondrial remap step can also use minimap2 for long reads, but should use a long-read preset and less short-read-specific split-segment assumptions.

Deletion evidence requires a valid breakpoint-spanning split or gapped alignment within one read. The workflow does not infer deletions from mate-pair distance alone. Paired-end identifiers are collapsed to fragment-level observations after candidate generation, but mates are never joined to create a deletion call.

The first-pass selection mode is controlled by `mapping.first_pass_read_selection`:

- `whole_genome_mt_best` maps against the full genome including mtDNA and passes reads with mitochondrial best/selected evidence to mitochondrial remapping. This is the default. With `mapping.first_pass_aligner: star`, this path streams STAR output through read-name collation into the selector and does not create full-genome BAMs.
- `nuclear_unmapped_reads` maps against a nuclear-only reference and passes unmapped reads to mitochondrial remapping. This mode provides a strict depletion-style sensitivity check but can discard real mitochondrial reads with NUMT-like nuclear alignments.
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
  first_pass_minimap2_index_extra: -k15 -w5
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
    replication_arc_boundary_basis: MitoBreak analytical major/minor arc convention using O_H=407 and O_L=5747 on the rCRS coordinate system.
    replication_arcs:
      - name: minor_arc
        display_name: Minor arc
        start: 408
        end: 5746
      - name: major_arc
        display_name: Major arc
        start: 5747
        end: 407
```

Local `.gz` files are decompressed into the workflow reference directory. A configured local path takes precedence over the corresponding URL. URL-based downloads remain supported when the local path field is omitted.
The optional `mt_reference_name` and `mt_reference_accession` fields are shown in the report so collaborators can confirm the mitochondrial coordinate standard used for deletion coordinates. When `replication_arcs` are configured, `replication_arc_boundary_basis` should state how the origin landmarks were chosen, and each arc should include an explicit `boundary_definition` documenting the coordinate convention.

Local sample FASTQs can be provided in the sample table with `fastq_1` and optional `fastq_2` columns. Both `.fastq.gz` and uncompressed `.fastq` are supported for local inputs.

The repository-root `inputs/` directory is reserved for local source data and reference files. It is excluded from Git; tracked sample tables and dataset configs may point to files there without adding the input data to the repository.

Resolved SRA metadata is cached under `metadata/cache/`. Once a cache file exists, metadata resolution uses it before trying NCBI again, so reruns after configuration changes do not require live NCBI metadata access.

## Report Outputs

Open `results/<dataset>/quality/report/index.html` to choose among the three evidence profiles. `standard` is the primary interpretation, `stringent` restricts the same canonical event set to strong evidence, and `exploratory` adds review-tier events. Each profile recalculates its own matrices, comparisons, plots, and ordinations.

The profile reports answer:

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
- deletion size distributions, including log-y and size-restricted views, with HTML bar/bin mouseovers for the size interval, group, plotted value, and contributing read count;
- deletion rainfall plots shown as full-size per-group figures by left breakpoint, right breakpoint, and circular deleted-interval midpoint, with interactive support filtering and point mouseovers in the HTML report;
- circular breakpoint-chord plots joining the directed breakpoints of each threshold-eligible exact deletion, with interactive normalized-support and observation controls in the HTML report;
- circular exact-deletion group-comparison plots with replicate-significance, exploratory replicate-p, and technical read-depth views plus optional effect, support, and direction refinements;
- breakpoint-pair support maps showing which deletion starts pair with which deletion ends;
- group-split pooled breakpoint support-density plots showing where deletion endpoints accumulate, with binned plotted support split into left and right breakpoint bars behind a circular-smoothed total-support curve;
- affected-feature normalized support and within-group proportions, with HTML bar mouseovers for category, group, and plotted value;
- collapsed feature-impact classes;
- per-gene affected burden and group-colored exact deletion recurrence, with HTML bar mouseovers for the feature or exact deletion, group, and plotted value;
- PCA and Bray-Curtis MDS for exact deletions and affected-feature categories, without static sample labels or centroids.
- mitochondrial gene-pair PCA for datasets with the configured short-read RNA STAR evidence stream.

Plot group order and colors are chosen once per report and reused across plots. For two-factor designs with `age` and `treatment`, groups are ordered by age and then treatment, with control-like groups first within each age. Control-like groups use subdued colors and treatment groups use red-family colors when possible.

The deletion rainfall plots, breakpoint-pair support map, and pooled breakpoint support-density plots are visual displays, not separate calling steps. They use the same display rule. By default `plots.rainfall_min_support_per_million: 0.0` and `plots.rainfall_max_points_per_group: 0`, so all eligible exact deletions are loaded and low-abundance calls are not hidden by a normalized-support cutoff or a fixed 300-call cap. The HTML rainfall and breakpoint-pair views provide logarithmic minimum-support sliders and mouseovers with the directed coordinates, support, affected features, arc context, and configured matches. The PDF remains a static all-call snapshot. These settings affect only the plots; the exact-deletion tables, matrices, comparisons, and read-list links retain the full analyzed call set subject to their own table-display filters.

In pooled breakpoint support-density plots, bar height and the smoothed curve use the configured plotted support metric (normally deletion support per million usable reads). A density bin's call count is the number of distinct exact-deletion calls contributing an endpoint, not the number of supporting reads. HTML hover metadata reports that call count, plotted support, and raw supporting observations separately; this distinction is important when one common deletion contributes many observations to a bin.

The baseline circular breakpoint-chord PDFs use the same support threshold and any explicitly configured per-group count cap as the rainfall plots. Their HTML views load every exact deletion passing the support threshold before an optional cap. The logarithmic support slider controls the minimum normalized support. The observation selector is linked in both directions: its `Auto` value reports the lowest raw count among calls passing the support and size filters; choosing a numeric observation cutoff moves the support slider to the lowest normalized support among calls meeting that observation and size cutoff; moving the slider returns the selector to `Auto`. Rainfall and breakpoint-pair views use the same interaction. Circular comparison plots load delivered exact-deletion comparison rows with at least one supporting observation across the two compared groups and provide named views for replicate-level BH significance, exploratory unadjusted replicate p-values, and technical read-depth enrichment. Zero-versus-zero rows remain in the complete comparison TSV but are not drawn. Display refinements do not create statistical significance.

Mitochondrial coordinate plots use a shared annotation palette: D-loop/control region coral, protein-coding genes green, rRNA cyan, and tRNA purple. Circular plots place coordinate 1 at 12 o'clock and increase coordinates clockwise. Feature and chord mouseovers expose names, coordinates, directed deletion IDs, support, arc annotations, and applicable comparison statistics.

Large result tables in the HTML report are searchable, sortable, and paged. Small tables are shown directly without search controls.

Heatmaps are intentionally not part of the main report.

## Repository Layout

- `Snakefile` - workflow rules.
- `config/defaults.yaml` - general workflow defaults.
- `config/datasets/*.yaml` - dataset-specific configuration.
- `envs/mitochondrial-deletions.yaml` - conda environment used by Snakemake rules.
- `scripts/` - Python workflow scripts.
- `report_assets/` - CSS and JavaScript inlined into generated interactive reports.
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

Use the dataset-specific dry-run command below before downloading data or building results. The explicit deliverables target and `--rerun-triggers mtime` keep the planned job list tied to the complete requested dataset output.

Snakemake can also create per-rule environments under `.snakemake/conda/` when `--use-conda` is used.

## Run The Rat GPA Dataset

This dataset contains 22 single-end Illumina RNA-seq runs from quadriceps muscle,
split across 18- and 34-month animals and Control and GPA treatment groups. The
fixed sample table is `metadata/rat_aging_muscle.samples.tsv`; it uses direct ENA
FASTQ URLs, so this dataset requires approximately 34 GB of network downloads for
a clean run.

Dry-run first:

```bash
snakemake results/rat_aging_muscle/rat_aging_muscle_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --resources download=2 \
  --configfile config/datasets/rat_aging_muscle.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds \
  --dry-run
```

Then run:

```bash
snakemake results/rat_aging_muscle/rat_aging_muscle_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --resources download=2 \
  --configfile config/datasets/rat_aging_muscle.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds
```

Deliverable reports:

```text
results/rat_aging_muscle/rat_aging_muscle_deliverables/index.html
results/rat_aging_muscle/rat_aging_muscle_deliverables_light/index.html
```

## Run The Human Common Deletion Dataset

The local common-deletion FASTQs are expected under:

```text
inputs/human_common_deletion/source_bundle/
```

The sample table is:

```text
metadata/human_common_deletion.samples.tsv
```

These can be raw BCLConvert FASTQs. The normal workflow stages local FASTQs, runs FastQC, trims with `fastp`, maps trimmed reads, then performs minimap2 mitochondrial remapping.

Dry-run first:

```bash
snakemake results/human_common_deletion/human_common_deletion_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --configfile config/datasets/human_common_deletion.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds \
  --dry-run
```

Then run:

```bash
snakemake results/human_common_deletion/human_common_deletion_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --configfile config/datasets/human_common_deletion.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds
```

Deliverable reports:

```text
results/human_common_deletion/human_common_deletion_deliverables/index.html
results/human_common_deletion/human_common_deletion_deliverables_light/index.html
```

## Run The Matched Human Bulk RNA-seq Dataset

This dataset contains one paired-end Illumina bulk RNA-seq sample matched to the
Nanopore dataset. Its sample table and config are:

```text
metadata/human_bulkseq_matched_nanopore.samples.tsv
config/datasets/human_bulkseq_matched_nanopore.yaml
```

Dry-run first:

```bash
snakemake results/human_bulkseq_matched_nanopore/human_bulkseq_matched_nanopore_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --configfile config/datasets/human_bulkseq_matched_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds \
  --dry-run
```

Then run:

```bash
snakemake results/human_bulkseq_matched_nanopore/human_bulkseq_matched_nanopore_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda --cores 8 \
  --configfile config/datasets/human_bulkseq_matched_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds
```

Deliverable reports:

```text
results/human_bulkseq_matched_nanopore/human_bulkseq_matched_nanopore_deliverables/index.html
results/human_bulkseq_matched_nanopore/human_bulkseq_matched_nanopore_deliverables_light/index.html
```

## Run The Human Nanopore Dataset

The nanopore dataset expects this local single-end FASTQ:

```text
inputs/human_nanopore/fastq/SUP_ONT_rCS.fastq
```

The sample table and config are:

```text
metadata/human_nanopore.samples.tsv
config/datasets/human_nanopore.yaml
```

This config disables FastQC and fastp trimming, stages the local uncompressed FASTQ as gzip, uses minimap2 for competitive whole-genome mitochondrial read selection, then remaps retained reads to normal and rotated mitochondrial references with `map-ont`.

Dry-run first:

```bash
snakemake results/human_nanopore/human_nanopore_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/human_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds \
  --dry-run
```

Then run:

```bash
snakemake results/human_nanopore/human_nanopore_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/human_nanopore.yaml \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --printshellcmds
```

Deliverable reports:

```text
results/human_nanopore/human_nanopore_deliverables/index.html
results/human_nanopore/human_nanopore_deliverables_light/index.html
```

## Advanced: Rerun From Existing Trimmed FASTQs

Use this only when downloads and trimming are already complete and you do not want Snakemake to backtrack into raw FASTQ staging or trimming rules.

Dry-run first:

```bash
snakemake results/rat_aging_muscle/rat_aging_muscle_deliverables_light/DELIVERABLES_COMPLETE.txt \
  --use-conda \
  --cores 8 \
  --configfile config/datasets/rat_aging_muscle.yaml \
  --config workflow_start_from=trimmed \
  --rerun-triggers mtime \
  --dry-run
```

Check that the dry-run does not include `prepare_reads`, `fastqc_raw`, or `trim_reads`. Then run the same command without `--dry-run`.

## Key Outputs

Report and canonical evidence outputs:

- `results/<dataset>/<dataset>_deliverables/index.html` - self-contained deliverable selector for the stringent, standard, and exploratory packages.
- `results/<dataset>/<dataset>_deliverables/profiles/<profile>/` - packaged report, tables, matrices, plots, and read lists for one profile.
- `results/<dataset>/<dataset>_deliverables/shared/` - packaged canonical evidence and provenance shared by all profiles.
- `results/<dataset>/<dataset>_deliverables.zip` - portable ZIP archive of the complete full deliverable folder.
- `results/<dataset>/<dataset>_deliverables_light/index.html` - shareable report selector with cluster-level tables, matrices, plots, methods, and configuration, excluding read lists and observation-level audit tables.
- `results/<dataset>/<dataset>_deliverables_light.zip` - portable ZIP archive of the complete light deliverable folder.
- `results/<dataset>/quality/report/index.html` - profile selector and profile counts.
- `results/<dataset>/quality/shared/` - canonical source, observation, cluster, tier, profile-membership, and resolved-configuration tables.
- `results/<dataset>/quality/profiles/<profile>/.report/index.html` - one complete profile report.
- `results/<dataset>/quality/profiles/<profile>/junctions/` - profile-filtered canonical exact deletions and observations.
- `results/<dataset>/quality/profiles/<profile>/plots/` - profile-specific plots.
- `results/<dataset>/quality/profiles/<profile>/matrices/` - profile-specific matrices.
- `results/<dataset>/quality/profiles/<profile>/analysis/` - profile-specific statistics and summaries.

Important machine-readable outputs:

- `quality/shared/source_candidates.tsv` - caller candidates and explicit pass/fail reasons.
- `quality/shared/canonical_observations.tsv` - cross-rotation and cross-caller deduplicated physical observations.
- `quality/shared/canonical_clusters.tsv` - stable exact deletion IDs, support, annotations, caller status, evidence tier, and quality flags.
- `quality/shared/ambiguous_direction_observations.tsv` - reciprocal-direction conflicts retained for audit.
- `quality/shared/report_profile_membership.tsv` - stable cluster membership in every report profile.
- `quality/profiles/<profile>/matrices/gene_pair_support_per_million.tsv` - short-read RNA gene-pair aggregation when applicable.
- `quality/shared/breakpoint_reference_support.tsv` - local remap reference-spanning support when applicable.

The `*_per_million_mt_reads.tsv` filenames are stable output names and do not identify the normalization denominator by themselves. Check the `normalization_denominator` and `normalization_reads` columns in `quality/profiles/<profile>/analysis/deletion_burden.tsv` and the report method table to determine whether the run used total usable reads or retained mitochondrial-evidence reads as the denominator.

For usability, each HTML profile report embeds at most 500 exact deletions by default. Configured deletion-target matches are prioritized, followed by the highest-support calls. No absolute supporting-read threshold is applied by default. The complete profile call set is in `quality/profiles/<profile>/junctions/junction_clusters.tsv`; the canonical all-tier set is in `quality/shared/canonical_clusters.tsv`. Adjust this display-only behavior with `report.exact_deletion_table` in the configuration.

When read-level evidence is available, the HTML report links read-count cells to sidecar TSVs in `read_lists/`. Exact-deletion support counts link to the reads supporting that deletion, configured deletion-target remap counts link to the reads supporting all nearby remap calls assigned to that target, and configured sequence-search counts link to reads containing the configured literal motif.

## Testing

```bash
conda activate mitochondrial-deletions
python -m pytest -q
```

## About `fasterq.tmp.*` Folders

`fasterq-dump` can create temporary folders named like `fasterq.tmp.<host>.<pid>` if it is not given an explicit temp path. The workflow wrapper passes an explicit temp directory inside each accession staging directory under `results/<dataset>/fastq/.<accession>.fasterq/tmp`. Top-level `fasterq.tmp.*` folders are scratch output and should not be committed.

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
