# Workflow Methods And Assumptions

## Purpose

This document explains how the workflow converts sequencing reads into
coordinate-focused mitochondrial deletion evidence, which assumptions are made at
each stage, and how interpretation changes with read technology and molecule type.
The generated HTML report presents the resolved settings and applicable caveats for
each run.

## Evidence Levels

The workflow distinguishes four levels:

1. **Split-alignment evidence** is an alignment pattern supporting a directed
   retained adjacency.
2. **Inferred deletion model** is the directed circular interval absent between those
   retained flanks under a deletion interpretation.
3. **Exact deletion** is a clustered coordinate-level inferred model with read
   support and rotation provenance.
4. **Affected-feature category** is an annotation derived from overlap between the
   directed deleted interval and configured reference features.

None of these terms, by itself, establishes that the source was a viable biological
mtDNA deletion molecule.

## Workflow Path

### 1. Inputs and metadata

Samples can use local FASTQ files or accessions. Metadata resolution prefers cached
local metadata. Dataset configuration identifies the reference, grouping columns,
read technology, molecule type, assay type, and library strategy. Unknown assay
properties must be reported as `unknown`, not inferred from sample names or mapper
presets.

### 2. Read preparation

Reads are staged and optionally trimmed with `fastp`. Reports state whether trimming
actually occurred and retain total usable reads for normalization. A rerun beginning
from existing intermediate data must not silently backtrack into download or
trimming rules.

### 3. First-pass read selection

The default `whole_genome_mt_best` mode competitively maps reads with nuclear and
mitochondrial references present together, then retains reads with selected
mitochondrial best evidence. This is read selection and provenance, not a separate
reported deletion-calling stream.

Available selection modes also include:

- `nuclear_unmapped_reads` retains reads not aligned to a nuclear-only reference;
- `mt_evidence_reads` retains the alternative full-genome mitochondrial-evidence scan.

NUMTs and ambiguous nuclear/mitochondrial placements remain relevant even after
selection.

### 4. Circular mitochondrial remapping

Retained reads are mapped with minimap2 to a normal mitochondrial reference and a
reference rotated by approximately half the mitochondrial genome. Both references
contain one nonduplicated mtDNA copy. Coordinates from every rotation are converted
back to the configured standard coordinate system.

The two rotations reduce dependence on an artificial linear origin. They do not
make reciprocal junction directions equivalent.

### 5. Directed split-alignment calling

Usable segments are ordered on the query. Same-strand query order is normalized into
a forward-reference retained adjacency `L -> R`:

- for plus-strand records, the earlier query segment ends at retained base `L` and
  the later segment begins at retained base `R`;
- for minus-strand records, the later query segment ends at retained base `L` and
  the earlier segment begins at retained base `R`.

The inferred deleted interval is the circular forward arc from `L` to `R`, excluding
both retained breakpoint bases. `R -> L` is the complementary deletion model and has
a different exact-deletion identifier.

The caller evaluates compatible query-segment pairs within one physical read
sequence. For paired-end data, SAM read1/read2 flags define separate alignment
chains even when both mates have the same query name. Supplementary records and
SA-tagged primary or secondary records remain eligible by default. Ordinary
single-segment records cannot form a split chain and are discarded before
in-memory grouping. These policies, MAPQ, anchor length, aligned fraction, soft
clipping, and
query overlap/gap thresholds are configuration-driven and should be evaluated in
sensitivity analyses.

### 6. Ambiguity handling

Reciprocal directions supported by the same read and unordered breakpoint pair,
including coordinates within configured circular breakpoint slop, are classified as
ambiguous. The default policy retains these rows in an audit table and excludes them
from primary clustering, matrices, burden, plots, and tests.

Different high-quality reads can support genuinely different reciprocal models.
Those models remain separate exact deletions rather than being collapsed by
breakpoint-pair identity.

### 7. Physical-observation consolidation and clustering

For configured short-read RNA datasets, STAR mitochondrial chimeric records and
minimap2 remap records enter one canonical evidence table after caller-specific
filters. The same sample, physical read or fragment, and directed breakpoint model
is counted once when it is observed by both callers. This prevents caller agreement
from inflating abundance while retaining caller provenance as corroborating
evidence. Paired-end mates are not joined to invent a breakpoint-spanning chain;
observations with the same fragment name are collapsed only after an individual mate
has supplied valid junction evidence.

Mate-placement fields are descriptive rather than calling criteria. When the
retained intermediates do not contain the required mate alignments, the canonical
observation table records `not_available_from_retained_intermediates` instead of an
empty value or an inferred concordance result.
For single-end libraries, no mate exists; mate evidence is neither required nor used
to accept an individually supported split or chimeric observation.

Directed left and right breakpoints are clustered within configured circular slop.
The same sample/read supporting the same directed cluster in multiple rotations is
counted once. Rotation support and whether an event has single- or multiple-rotation
evidence are reported. STAR-derived coordinates are already in the configured
full-genome mitochondrial coordinate system and do not acquire artificial rotation
support.

An unordered `breakpoint_pair_id` is retained only for diagnostics. It is not the
exact-deletion identifier.

### 8. Annotation, evidence tiers, and report profiles

The directed interval determines:

- affected features;
- fully removed and partially overlapped features;
- nearest and breakpoint-overlapping features;
- origin-spanning and control-region involvement;
- overlap with the configured major and minor mitochondrial replication arcs;
- size class;
- configured deletion-target matches.

Major/minor replication-arc overlap is a reference annotation applied after the
directed deleted interval has been inferred. `replication_arc_context` distinguishes
`major_arc_only`, `minor_arc_only`, and `major_and_minor_arcs`, while
`major_arc_deleted_bp` and `minor_arc_deleted_bp` give the corresponding overlap
lengths. Replication arcs are configured under
`references.<species>.replication_arcs`; they are not affected features and do not
participate in reciprocal interval selection.

Expected mitochondrial transcript junctions are configuration-driven. When enabled,
matching transcript-compatible evidence is excluded from primary deletion summaries
but retained in QC.

Every canonical exact deletion receives one evidence tier. With the default
thresholds:

- `strong` has at least two distinct physical observations and either repeated
  support within a sample or a physical observation corroborated by both callers;
- `supported` has at least two distinct observations but lacks the additional
  within-sample or cross-caller condition;
- `review` has fewer than two distinct observations;
- `rejected` contains a configured expected-transcript junction or unresolved
  reciprocal-direction ambiguity.

These are evidence-management tiers, not declarations that a call is a confirmed
biological mtDNA deletion. Thresholds and profile membership are configuration
fields. The workflow creates three independent analysis/report views from the same
stable exact-deletion identifiers:

- `stringent` includes `strong`;
- `standard` includes `strong` and `supported` and is the primary report;
- `exploratory` includes `strong`, `supported`, and `review`.

Matrices, PCA/MDS coordinates, statistics, tables, and plots are recalculated within
each profile. Ordination axes therefore must not be compared as though they were the
same axes across profiles. The report index gives cluster and physical-observation
counts for every profile without duplicating the canonical source tables.

### 9. Normalization and local reference support

Primary support can be normalized per million total usable reads or per million
first-pass mitochondrial-evidence reads. The selected denominator is reported and
does not validate a candidate junction.

For exact deletions with minimap2 evidence, local reference-spanning support counts
primary remap alignments covering the undeleted reference around each breakpoint.
The local split-support fraction is a breakpoint-neighborhood alignment metric, not
automatically mtDNA heteroplasmy. STAR-only calls are marked unavailable because a
STAR chimeric numerator cannot be compared directly with a minimap2 remap
reference-support denominator.

### 9.1 Circular deletion displays

Circular breakpoint-chord plots are report views over the canonical exact-deletion
objects; they are not additional callers. Each chord joins the alignment-directed
left and right retained breakpoints. Coordinate 1 is at 12 o'clock and coordinates
increase clockwise. The outer annotation ring and the linear mitochondrial feature
tracks use the same feature classes and colors: D-loop/control region coral,
protein-coding genes green, rRNA cyan, and tRNA purple.

The baseline chord PDFs use the rainfall support threshold and any explicitly
configured per-group count cap. The HTML location view loads every call passing
that support threshold before the count cap. Its logarithmic support slider controls
minimum normalized support. Rainfall HTML views load the complete eligible call set
by default and provide the same support slider plus point mouseovers; rainfall PDFs
remain static all-call snapshots. The
observation selector can add an absolute raw-evidence cutoff; `Auto` reports the
lowest raw count among calls retained by the support slider, and moving that slider
returns the selector to `Auto`. Reset restores the baseline PDF display set.

Pooled breakpoint support-density plots use the configured plotted support metric
for bar heights and circular smoothing. The endpoint-count values shown in HTML
hover metadata count distinct exact-deletion calls in each bin; they are not read
counts. Raw supporting observations and plotted support are reported as separate
fields so high-support calls can be distinguished from bins containing many
low-support calls.

Circular group-comparison views load delivered exact-deletion comparison rows with
at least one supporting observation across the two compared groups. Zero-versus-zero
rows remain in the complete comparison TSV but do not produce chords.
The replicate-significant view applies replicate-level BH q <= 0.05 and is the
appropriate preset for biological group conclusions. The unadjusted replicate-p
view is exploratory. The read-depth BH view is technical count evidence rather than
biological-replicate significance. Optional effect-size, supporting-observation, and
direction controls refine only what is drawn. HTML mouseovers expose the underlying
coordinates, exact-deletion ID, support, arc annotation, and applicable statistics;
the static PDF remains the portable baseline view.

### 10. Short-read RNA caller roles

STAR maps against the full nuclear-plus-mitochondrial genome. Its chimeric junction
records can recover short breakpoint-spanning RNA evidence that is difficult for a
mitochondrial-only minimap2 remap to represent. STAR candidates must have two
mitochondrial anchors, pass configured query geometry and size limits, and, when
configured, anchor to distinct annotated mitochondrial genes. Same-gene chimeric
records, missing gene anchors, expected transcript junctions, and other failed rows
remain visible in the source-candidate audit with their filter reason.

Minimap2 remains the circular mitochondrial remap caller and supplies rotation-aware
split-alignment evidence and local reference-spanning denominators. STAR does not
replace remapping, and STAR-Fusion is not run. Exact coordinate identity remains the
directed breakpoint model; `star_gene_pair_label` and
`breakpoint_flanking_gene_pair_label` are annotations used for interpretation and
gene-pair aggregation, not substitutes for exact deletion IDs.

The mitochondrial gene-pair PCA is generated only when the short-read RNA STAR
evidence stream is enabled. It aggregates retained observations by STAR annotated
gene pair where available and uses breakpoint-flanking gene pairs for other evidence
within those dual-caller datasets. Like every ordination, it is rebuilt separately
for each report profile.

## Coordinate Convention

Breakpoints are 1-based retained flanking bases:

- non-wrapping `L -> R` deletes `L + 1` through `R - 1`;
- wrapping `L -> R` deletes `L + 1` through the mtDNA end and position `1` through
  `R - 1`;
- deleted size excludes both retained breakpoint bases;
- `wraps_origin` describes the directed deleted interval, not the physical location
  of the read or an unordered breakpoint pair.

The *deleted circular arc* in the caller is event-specific and follows the retained
adjacency `L -> R`. The mitochondrial *major arc* and *minor arc* are fixed
reference regions bounded by configured replication-origin landmarks. These are
different uses of the word "arc." Major/minor arc labels describe where the inferred
deleted bases fall; they cannot determine whether `L -> R` or `R -> L` is supported.

## Assumptions

Primary results assume that:

- the configured mitochondrial reference and coordinate standard are appropriate;
- accepted split segments belong to one original molecule or library fragment;
- query order and strand are interpreted consistently by STAR or minimap2, their
  parsers, and the caller;
- alternative placements, repeats, and NUMTs do not better explain the retained
  evidence;
- the read adjacency is not a basecalling, ligation, PCR, reverse-transcription,
  template-switch, or other library artifact;
- directed breakpoints within clustering slop represent the same coordinate-level
  event;
- read support measures detected evidence rather than viability, molecule frequency,
  or heteroplasmy;
- absence of evidence is not evidence of absence.

The report presents applicable alternative explanations in a configuration-driven
table organized by sequencing technology, molecule type, and assay type. The table
distinguishes technical artifacts and alignment ambiguity from biological RNA
processing that can create deletion-like evidence without an mtDNA deletion. These
entries are interpretation considerations, not classifications of every reported
call as artifactual.

## Technology-Specific Interpretation

### Nanopore

Long reads can provide long anchors and direct junction-spanning evidence. Base
errors, homopolymers, supplementary alignments, alternative placements, chimeric
reads, ligation artifacts, and concatemers can affect breakpoint precision and call
validity. Alignment-chain and secondary-alignment policies materially affect the
candidate set.

### Illumina

Short split anchors can be difficult to place uniquely around repeats and NUMTs.
Read length and breakpoint position affect sensitivity. The workflow does
not infer deletions from mate distance alone. PCR duplication and library chimeras
can inflate support unless addressed by the experimental design or a validated
deduplication method.

## Molecule-Type-Specific Interpretation

### RNA

Mitochondrial transcript processing, polycistronic RNA structure, alignment behavior,
reverse transcription, and template switching can create deletion-like junctions.
Configured expected-transcript filtering covers only the supplied model. RNA support
does not directly measure mtDNA genome frequency or heteroplasmy.

### DNA

DNA reads are closer to direct genome-molecule evidence but can reflect NUMTs,
PCR or ligation chimeras, mapping ambiguity, and sampling. Control-region-lacking
molecules may be defective, complemented, transient, or artifactual; viability must
not be inferred solely from split alignments.

### Unknown

When molecule type is unknown, the report states that DNA/RNA-specific biological
interpretation is limited. It must not select one set of assumptions heuristically.

## Validation Options

For prioritized calls, inspect the existing BAM records and compare:

- query spans and CIGARs;
- strand and primary/supplementary/secondary flags;
- MAPQ, alignment score, edit distance, and `SA`/minimap2 tags;
- directed-junction agreement across rotations.

Synthetic templates provide a direct diagnostic:

- `L|R` for the reported directed deletion;
- `R|L` for the complementary model;
- a wild-type circular-origin template.

A doubled or padded mitochondrial reference can diagnose artificial linear-boundary
behavior but creates repeated-reference placement ambiguity. It is supplementary and
does not replace the normal plus rotated primary design.

## Result Provenance

The deliverable entry point is
`results/<dataset>/<dataset>_deliverables/index.html`. It links self-contained
stringent, standard, and exploratory packages. Each profile package contains its
report, tables, matrices, plots, and deletion read lists. Shared canonical evidence,
profile membership, tier summaries, and resolved configuration are included under
the deliverable `shared/` and `config/` directories. A matching
`<dataset>_deliverables.zip` containing the top-level full deliverable directory is
created by default.

The workflow also creates
`results/<dataset>/<dataset>_deliverables_light/index.html` and a matching
`<dataset>_deliverables_light.zip` by default. The light package retains each
profile's report, cluster-level exact-deletion and analysis tables, matrices, plots,
methods, configuration, and tier summaries. It excludes per-deletion read lists,
source-candidate and canonical-observation tables, ID maps, and intermediate audit
tables. Read-list links and instructions are removed from the light reports. The ZIP
contains the top-level light deliverable directory. Both archives can be opened with
standard ZIP tools.

The corresponding workflow working outputs remain under
`results/<dataset>/quality/`: `report/index.html` is the working profile selector,
`shared/` contains canonical outputs, and `profiles/<profile>/` contains the
profile-specific analysis tree.

The embedded exact-deletion report table has no absolute supporting-read threshold
by default. It displays all calls when at most 500 are present; otherwise it
prioritizes configured targets and then the highest-support calls up to the 500-row
display cap. The profile-specific `tables/exact_deletions.tsv` remains unfiltered.

Canonical deliverables include:

- `shared/source_candidates.tsv`, including passing and failed source rows;
- `shared/canonical_observations.tsv`, one row per deduplicated physical
  observation and exact-deletion assignment;
- `shared/canonical_clusters.tsv`, one row per canonical exact deletion;
- `shared/ambiguous_direction_observations.tsv`;
- `shared/report_profile_membership.tsv`;
- `shared/quality_tier_summary.tsv`;
- `shared/resolved_quality_config.yaml`;
- profile-specific `junctions/`, `matrices/`, `analysis/`, `plots/`, and `.report/`
  directories.

Canonical observation and cluster tables retain directed and pre-cluster
coordinates, caller and rotation provenance, strand and query geometry, support,
quality tier and flags, exact IDs, feature annotations, and applicable alignment
metrics. Empty or unavailable fields mean the metric does not apply to that caller
or cannot be derived from the available intermediate evidence; they must not be
treated as zero.
