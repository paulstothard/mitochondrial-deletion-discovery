# Directed Circular Deletion Arc Workflow Fix

## Status

Implemented on 2026-07-10. The corrected workflow uses result schema
`2.0-alignment-directed-arcs`. The human nanopore results were rebuilt downstream
from the existing normal and half-rotation BAMs after a dry-run confirmed that no
download, trimming, first-pass mapping, or mitochondrial-remapping rules were
scheduled.

Validation completed during implementation:

- 52 unit/regression tests passed;
- reciprocal directions remain separate and directed arcs longer than half the
  mitochondrial genome are retained;
- secondary alignments are excluded from primary calls by default;
- same-read reciprocal conflicts within configured breakpoint slop are retained as
  ambiguous evidence and excluded by default;
- tables, matrices, plots, HTML report, read lists, methods provenance, and data
  dictionary were regenerated;
- regenerated location, size, and affected-feature plots and the rendered report were
  inspected;
- the previous nanopore deliverable was preserved as
  `human_nanopore_legacy_shortest_arc_deliverables_20260710`.

Candidate-specific synthetic `L|R`/`R|L` template remapping and doubled-reference
diagnostics were not launched automatically across all corrected calls. They remain
orthogonal validation steps for a configured set of prioritized candidates; applying
them to every call without an explicit selection and repeated-reference placement
policy would be expensive and potentially misleading. This does not affect the
directed-arc correction, but it means the regenerated coordinate evidence should not
be described as biological validation of every event.

## Scope

Correct the workflow so that a circular mitochondrial deletion interval is assigned
from the directed breakpoint junction supported by split-alignment query order and
strand. The workflow must not replace that interval with the shorter of the two
possible circular arcs.

The work also needs to make the complete analysis path and its assumptions visible
in reports and documentation. The workflow supports nanopore and Illumina data and
may be used with DNA or RNA libraries, so explanations and caveats must be selected
from resolved configuration rather than inferred from a dataset name.

This plan separates two issues:

1. Correctness of circular arc assignment and downstream interpretation.
2. Confidence in the split alignments used as deletion evidence.

The arc correction should be implemented and tested first. Alignment-evidence
hardening should be a separate, reviewable change so changes in event meaning can be
distinguished from changes in call count.

## Current Problem

The minimap2 caller currently derives a directed breakpoint pair from query order
and strand, but then canonicalizes the pair by choosing the shorter circular arc.
That operation can replace the alignment-directed deleted interval with its
complement. Reciprocal junctions are subsequently merged, annotated, plotted, and
reported as the same exact deletion.

As a result:

- an origin-spanning interval can be reported even when the split alignment implied
  the non-origin-spanning complementary interval;
- affected-feature and control-region annotations can describe the assumed shorter
  arc rather than the alignment-directed arc;
- final cluster tables and plots cannot, by themselves, establish which directed
  junction the supporting alignments contained;
- reciprocal directed junctions that represent different deletion models can be
  combined.

## Corrected Deletion Semantics

### Directed junction

After segments are ordered on the query and normalized for strand, a junction
`L -> R` means that retained reference base `L` is adjacent to retained reference
base `R` in the read-supported model.

Using the workflow's existing flanking-base coordinate convention:

- the deleted interval is the circular forward arc from `L` to `R`;
- the deleted bases are exclusive of the retained breakpoint bases;
- a non-wrapping interval contains `L + 1` through `R - 1`;
- an origin-spanning interval contains `L + 1` through the mitochondrial genome end
  and position `1` through `R - 1`;
- `R -> L` is the complementary deletion model, not an interchangeable
  representation of `L -> R`.

### Canonicalization

Canonicalization should mean:

- conversion from each rotated reference into the configured standard
  mitochondrial coordinate system;
- consistent 1-based retained-flanking-base coordinates;
- stable directed deletion identifiers;
- deduplication of equivalent evidence across rotations.

Canonicalization must not mean selecting the shorter circular arc.

### Ambiguity

If usable alignments or reference rotations support both `L -> R` and `R -> L`, the
workflow should retain that conflict explicitly. It should not resolve the conflict
using interval length. Ambiguous events should be excluded from primary deletion
burden by default while remaining available in QC and evidence tables.

## Implementation Plan

### 1. Introduce an explicit result schema version

- Add a workflow result-schema version to the resolved configuration and report.
- Identify the new arc convention as `alignment_directed`.
- Treat existing shortest-arc results as legacy results with a different semantic
  version.
- Prevent legacy and directed exact-deletion matrices from being compared as though
  their column identifiers have the same meaning.

### 2. Correct the minimap2 deletion caller

- Preserve the existing query-order and strand normalization as the basis for the
  directed junction, after adding targeted tests for both strands.
- Remove shortest-arc selection from `deletion_from_segments()`.
- Calculate `deleted_size`, `wraps_origin`, and `deleted_interval` directly from the
  directed pair.
- Generate `exact_deletion_id` from directed breakpoints and directed deleted size.
- Populate the currently empty `reported_deleted_size` field or replace it with
  clearer directed-arc fields.
- Add diagnostic fields:
  - `arc_assignment_method`;
  - `complement_deleted_size`;
  - `complement_wraps_origin`;
  - `direction_status`;
  - `alignment_chain_id` or equivalent provenance identifier.
- Preserve enough per-segment evidence to audit a call without reconstructing caller
  internals:
  - query start and end;
  - reference start and end before and after rotation conversion;
  - CIGAR;
  - strand;
  - primary, supplementary, and secondary flags;
  - MAPQ;
  - available alignment-score, edit-distance, `SA`, and minimap2 tags.

### 3. Harden split-alignment evidence in a separate change

- Pair adjacent query segments instead of testing every pair of compatible segments.
- Prefer coherent primary/supplementary chains.
- Do not include secondary alignments in primary calls by default.
- Keep secondary-alignment sensitivity analysis configuration-driven and label its
  results clearly.
- Require consistent read identity, strand, query ordering, and acceptable query
  overlap/gap.
- Add configurable alignment-score and MAPQ requirements suitable for the selected
  read technology.
- Record why every candidate was accepted, rejected, or classified as ambiguous.
- Do not infer deletions from Illumina mate distance alone; retain the existing
  requirement for breakpoint-spanning split evidence unless a separately validated
  caller is added.

### 4. Correct clustering and cross-rotation deduplication

- Cluster directed `L -> R` events only with other directed `L -> R` events within
  configured breakpoint slop.
- Do not merge `L -> R` with `R -> L`.
- Deduplicate a read across the normal and half-rotated references only when the
  normalized directed junction agrees.
- Record rotation agreement and disagreement counts.
- Add an unordered `breakpoint_pair_id` only as a diagnostic grouping key; do not use
  it as the exact-deletion identifier.
- Use circular breakpoint distance near the coordinate origin when applying slop.
- Recompute the representative directed interval after clustering and retain each
  read's original directed coordinates.

### 5. Update annotation, matrices, analysis, and plotting

- Derive affected features, fully removed features, partial overlaps, nearest
  features, size class, control-region involvement, and configured-target matches
  from the directed deleted interval.
- Keep ambiguous-direction evidence out of primary exact-deletion and
  affected-feature matrices by default.
- Make inclusion of ambiguous evidence an explicit sensitivity option.
- Use directed coordinates for rainfall plots, breakpoint-pair plots, recurrence
  labels, interval midpoints, and origin-spanning markers.
- Replace wording such as "shorter canonical interval" with
  "alignment-directed circular interval."
- Ensure exact deletions and affected-feature categories remain distinct result
  levels.

### 6. Preserve a labelled legacy mode only if needed

- If reproducibility of historical results is required, provide an explicit
  `legacy_shortest_arc` option.
- Never make legacy mode the default.
- Add a prominent report warning when legacy mode is active.
- Do not describe legacy shortest-arc output as alignment-directed.
- Prefer a separate result directory or semantic version so legacy and corrected
  outputs cannot overwrite one another silently.

## Report Transparency Requirements

Every generated report should contain a concise but complete methods and assumptions
section built from the resolved configuration. It must explain what this particular
run did, not merely describe all capabilities the workflow might have used.

### 1. Run identity and provenance

Report:

- workflow commit or version when available;
- result-schema version;
- dataset title and sample count;
- configured species;
- mitochondrial reference name, accession, length, and contig names;
- annotation source;
- read technology;
- molecule type;
- library strategy when supplied;
- paired or single-end layout;
- local versus downloaded inputs without exposing unnecessary private paths;
- resolved configuration and a link to its delivered copy.

Add explicit configuration fields where needed rather than inferring assay properties
from sample names or mapper presets. Suggested fields are:

- `dataset.read_technology`: for example `nanopore`, `illumina`, or `other`;
- `dataset.molecule_type`: `dna`, `rna`, or `unknown`;
- `dataset.library_strategy`: optional free text or a controlled value;
- `dataset.pcr_amplified`: `true`, `false`, or `unknown` when relevant.

### 2. Workflow path used in the run

Describe the actual configured stages in order:

1. How reads were staged or downloaded.
2. Whether reads were trimmed and which quality decisions were applied.
3. Which first-pass aligner and selection mode were used.
4. What "mitochondrial-evidence reads" means for that selection mode.
5. Which mitochondrial references and rotations were used for remapping.
6. Which minimap2 preset, index-sensitive settings, and mapping-time settings were
   active.
7. How split segments were selected and joined into candidate evidence.
8. How direction, deleted interval, circular wrapping, clustering, rotation
   deduplication, and ambiguity were handled.
9. Which transcript-processing and other artifact filters were applied.
10. How exact deletions were annotated and converted into affected-feature
    categories.
11. Which denominator was used for normalized support.
12. Which calls were included or excluded from primary plots and statistical tests.

The text must change when configuration changes. For example, a report must not say
that reads were trimmed, secondary alignments were excluded, or ambiguous calls were
filtered unless those actions occurred in that run.

### 3. Explicit arc-assignment explanation

Include a small diagram or text example explaining:

- breakpoint coordinates identify two positions on a circle;
- query order and strand establish the directed retained adjacency;
- the forward circular arc from the directed left breakpoint to the directed right
  breakpoint is the inferred deleted interval;
- the complementary arc is a different hypothesis;
- origin spanning means the directed interval passes through the configured
  coordinate origin;
- interval coordinates refer to retained flanking bases, producing the documented
  one-base-excluded size convention.

The report should provide these fields in the exact-deletion table or linked evidence
table:

- directed left and right breakpoints;
- deleted size and interval;
- wraps-origin status;
- complement size;
- arc-assignment method;
- direction-consensus status;
- support by reference rotation;
- supporting read-list link.

### 4. Evidence terminology

Use the following distinction consistently:

- **Split-alignment evidence**: one or more alignments that support a directed
  breakpoint adjacency.
- **Inferred deletion model**: the directed circular interval implied by that
  adjacency after filtering and coordinate conversion.
- **Exact deletion**: a clustered coordinate-level inferred deletion model.
- **Affected-feature category**: a deterministic annotation derived from the
  directed deleted interval.
- **Biological deletion**: a biological interpretation that requires appropriate
  filtering, context, and ideally orthogonal validation.

Reports must not imply that every minimap2 split or supplementary alignment is a
confirmed biological mtDNA deletion.

### 5. Assumptions disclosed in every report

State whether the primary analysis assumes that:

- the configured mitochondrial reference and coordinate standard are appropriate;
- split segments belong to the same original molecule or library fragment;
- query order and strand have been interpreted consistently by the aligner and
  caller;
- retained alignments are sufficiently unique and are not better explained by NUMTs,
  repeats, or alternative placements;
- the alignment chain represents a physical adjacency rather than a sequencing,
  basecalling, ligation, PCR, reverse-transcription, or template-switch artifact;
- breakpoint clustering within the configured slop represents the same exact event;
- read support is a measure of detected evidence and not automatically molecule
  frequency or heteroplasmy;
- the configured normalization denominator answers the intended comparison question;
- absence of detected evidence is not proof that a deletion is absent.

### 6. Read-technology-specific explanation

For nanopore data, explain:

- long reads may span deletion junctions with long anchors;
- elevated base-error rates, homopolymers, supplementary alignments, and alternative
  placements can affect breakpoint precision;
- chimeric reads, ligation artifacts, concatemers, and library preparation can create
  apparent junctions;
- secondary-alignment handling and alignment-chain validation materially affect the
  candidate set.

For Illumina data, explain:

- short anchors provide less unique placement around repeats and NUMTs;
- a split read can support a breakpoint junction, while mate distance alone is not a
  call in the current workflow;
- PCR duplication and library chimeras can inflate support unless handled by the
  experimental design or a validated deduplication step;
- breakpoint precision and sensitivity depend strongly on read length and breakpoint
  position within the read.

Do not present nanopore-specific or Illumina-specific caveats when the corresponding
technology was not used.

### 7. Molecule-type-specific explanation

For RNA-derived data, explain:

- split alignments may arise from transcript processing, polycistronic RNA structure,
  splicing-like alignment behavior, reverse transcription, or template switching;
- configured expected transcript-junction filtering is limited to the supplied model;
- RNA abundance and transcript stability affect support;
- deletion-supporting RNA reads do not directly measure mtDNA heteroplasmy or genome
  copy fraction.

For DNA-derived data, explain:

- detected reads are closer to direct genome-molecule evidence but can still reflect
  NUMTs, PCR or ligation chimeras, mapping ambiguity, and sampling effects;
- local split-support fraction is a breakpoint-neighborhood alignment metric, not a
  complete heteroplasmy estimate unless additional dataset assumptions are justified;
- control-region-lacking molecules may represent low-abundance defective molecules or
  artifacts and should not be declared viable or nonviable from alignment evidence
  alone.

If molecule type is unknown, the report should state that assay-specific biological
interpretation is limited rather than selecting DNA or RNA language heuristically.

### 8. Filtering, normalization, and detection limits

Display the resolved values for the most interpretation-sensitive settings:

- minimum anchor length;
- minimum and maximum directed deletion size;
- minimum split-read support;
- breakpoint clustering slop;
- minimum MAPQ;
- secondary and supplementary alignment policy;
- query overlap and gap limits;
- minimum aligned fraction and maximum soft-clipping fraction;
- expected transcript-junction filtering;
- ambiguous-direction policy;
- normalized-support denominator;
- plot-only display thresholds and point caps.

Clearly separate call filters from plot-only display filters. State that normalization
changes comparability of support values but does not validate a candidate junction.

### 9. QC and warning summaries

Add report warnings when applicable:

- legacy shortest-arc mode was used;
- arc direction is ambiguous for any displayed event;
- only one rotation supports an event;
- secondary alignments contributed to primary calls;
- MAPQ 0 alignments contributed to primary calls;
- expected transcript filtering was not configured for an RNA dataset;
- molecule type or read technology is unknown;
- a group contains too few samples for inferential statistics;
- configured normalization counts are missing or zero;
- a result was produced under a schema older than the report code expects.

Warnings must describe what occurred and its interpretive consequence without making
dataset-specific biological conclusions.

### 10. Machine-readable provenance and data dictionary

- Deliver a machine-readable run-methods table containing resolved settings and short
  definitions.
- Deliver a data dictionary for canonical deletion, evidence, matrix, analysis, and
  QC columns.
- Include the result-schema version and arc-assignment method in canonical tables.
- Document which columns are workflow metadata, biological metadata, coordinates,
  evidence counts, normalization denominators, or derived annotations.
- Make report read lists include the directed and pre-cluster coordinates, strand,
  rotation, direction status, and enough alignment provenance for audit.

## Repository Documentation Updates

### README

- Replace the shortest-arc description with directed-junction semantics.
- Explain normal and half-rotated remapping and why rotation does not remove
  directional information.
- Document optional doubled-reference validation as a diagnostic, not the primary
  caller.
- Describe synthetic `L|R`, `R|L`, and wild-type junction validation.
- Add a compact technology-by-molecule-type interpretation table.
- Document the result-schema and legacy-mode policy.

### Configuration documentation

- Document all assay identity, arc assignment, ambiguity, alignment-chain, and
  evidence-filter settings.
- Identify which minimap2 settings affect index construction.
- Provide general nanopore-DNA, nanopore-RNA, Illumina-DNA, and Illumina-RNA example
  configurations without hard-coded sample names or biological conclusions.
- State defaults and the scientific consequence of changing them.

### Output documentation

- Document the retained-flanking-base coordinate convention with examples.
- Define every canonical deletion and evidence field.
- Explain exact deletions versus affected-feature categories.
- Explain raw support, normalized support, local reference-spanning support, and
  configured literal sequence searches as distinct quantities.
- Document how old shortest-arc outputs can be recognized.

### Circular-coordinate audit

- Extend the audit beyond coordinate round trips and wrapping arithmetic.
- Audit query order and strand for directed junctions.
- Verify direction agreement across rotations.
- Include worked reciprocal examples and show that they remain separate.
- Replace any conclusion that treats internally consistent shortest-arc wrapping as
  proof of the alignment-directed deleted interval.

## Test Plan

### Unit tests

- A plus-strand `13500 -> 1000` junction remains origin spanning.
- A plus-strand `1000 -> 13500` junction remains non-wrapping and is not merged with
  its reciprocal.
- Reverse-strand representations of each junction produce the same directed call as
  the corresponding plus-strand representation.
- Normal and half-rotated coordinates normalize to the same directed junction.
- A directed deletion larger than half the mitochondrial genome is retained when it
  passes configured size filters.
- Complement size is calculated but does not change the call.
- Same-read support across rotations is deduplicated only when direction agrees.
- Conflicting directions are marked ambiguous.
- Circular breakpoint slop works across coordinate `1`.
- Origin-spanning affected-feature annotation uses both correct interval pieces.
- Directed exact-deletion identifiers differ for reciprocal junctions.

### Caller integration tests

Create small synthetic circular references and reads covering:

- a non-wrapping deletion shorter than half the genome;
- a non-wrapping deletion longer than half the genome;
- an origin-spanning deletion;
- both read strands;
- wild-type reads spanning the artificial linear origin;
- alternative secondary placements;
- conflicting alignment chains.

Verify normal and rotated remapping, caller output, clustering, annotation, and
deduplication together.

### Report and plotting tests

- Generate reports for representative nanopore-DNA, nanopore-RNA, Illumina-DNA, and
  Illumina-RNA configurations.
- Verify that only applicable technology and molecule-type caveats appear.
- Verify that the workflow-path narrative matches resolved configuration.
- Verify that plots use directed biological coordinate columns rather than workflow
  metadata.
- Inspect origin-spanning plots and feature annotations visually.
- Confirm group colors and order remain consistent.
- Confirm legends and labels are not clipped or misleading.
- Confirm empty panels state the actual reason data are unavailable.
- Confirm report text does not describe ambiguous or unfiltered evidence as confirmed
  biological deletions.
- Confirm the data dictionary matches current table columns.

## Validation Plan For Existing Nanopore Results

### Stage 1: Re-call existing alignments

- Re-run the corrected caller against the existing normal and half-rotation BAMs.
- Do not download, trim, or remap raw reads for the initial arc comparison.
- Compare old shortest-arc identifiers with new directed identifiers.
- Summarize how many old origin-spanning calls were:
  - alignment-directed origin-spanning calls;
  - complementary arcs created by shortest-arc canonicalization;
  - directionally ambiguous;
  - supported by only one reference rotation.

### Stage 2: Inspect alignment evidence

- Inspect high-support and biologically consequential origin-spanning candidates.
- Verify query spans, CIGARs, strand, flags, MAPQ, alignment scores, and alignment-chain
  tags in both rotation BAMs.
- Confirm that supporting reads contain the same directed retained adjacency.
- Determine whether secondary alignments or incompatible chains generated candidate
  pairs.

### Stage 3: Synthetic junction validation

For each prioritized breakpoint pair, construct:

- an `L|R` directed-junction template;
- an `R|L` complementary-junction template;
- a wild-type circular-origin template.

Require continuous junction-spanning alignment with adequate anchors and a meaningful
score advantage over competing templates. Keep thresholds configuration-driven.

### Stage 4: Optional doubled-reference diagnostic

- Use a doubled or padded mitochondrial reference to test artificial linear-boundary
  behavior.
- Apply an explicit central-copy or equivalent coordinate policy.
- Track secondary placement caused by repeated reference sequence.
- Treat this as supplementary validation; do not replace the normal plus half-rotation
  primary design solely with a doubled reference.

## Rerun And Delivery Plan

- Use existing mitochondrial remap BAMs for the first corrected downstream rerun.
- Before recommending or executing a Snakemake rerun, dry-run the exact command.
- Verify that its job list does not include unintended downloads, read preparation,
  trimming, first-pass mapping, or mitochondrial remapping.
- Write corrected results to a versioned or otherwise isolated output location until
  validation is complete.
- Regenerate tables, matrices, statistics, plots, report, read lists, configuration,
  data dictionary, and deliverables.
- Test code and inspect the generated report and plots before accepting the work.
- Preserve the previous report as a clearly labelled legacy artifact for comparison.

## Acceptance Criteria

The correction is complete only when:

1. No production code chooses a deleted arc because it is shorter.
2. Directed reciprocal junctions remain distinct exact deletions.
3. Plus and minus strands and all configured rotations yield consistent directed
   coordinates for synthetic truth cases.
4. Ambiguous direction is represented explicitly and excluded from primary summaries
   by default.
5. Affected features and origin-spanning status are derived from the directed arc.
6. Primary calls use a documented, configurable alignment-chain policy.
7. Reports state the actual workflow path, settings, assumptions, and assay-specific
   limitations for the run.
8. Reports distinguish coordinate evidence, inferred deletion models, and biological
   conclusions.
9. Machine-readable tables contain the provenance needed to audit arc assignment.
10. The full test suite passes and generated plots and reports have been inspected.
11. The corrected nanopore rerun is demonstrated not to invoke unintended upstream
    download or trimming work.

## Proposed Implementation Sequence

1. Add directed-arc unit tests that fail under the current shortest-arc behavior.
2. Correct caller semantics and result identifiers.
3. Correct clustering, deduplication, annotation, matrices, and plots.
4. Add schema versioning and legacy-result safeguards.
5. Add report methods, assumptions, warnings, and machine-readable provenance.
6. Update README, configuration documentation, output documentation, and audit.
7. Re-call existing BAMs and inspect corrected nanopore outputs.
8. Harden alignment-chain selection in a separate reviewed change.
9. Run synthetic-junction and optional doubled-reference validation.
10. Produce and inspect final versioned deliverables.
