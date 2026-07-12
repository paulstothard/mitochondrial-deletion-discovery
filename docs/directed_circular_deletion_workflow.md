# Directed Circular Deletion Workflow Specification

## Purpose

This document defines the workflow requirements for inferring coordinate-focused
mitochondrial deletion evidence from split alignments on a circular reference. It
applies to Nanopore and Illumina reads and to DNA, RNA, or unspecified molecule
types. Dataset-specific biological conclusions are outside its scope.

The workflow must keep two questions separate:

1. Which directed circular interval is implied by a split-alignment chain?
2. How much confidence should be assigned to that alignment evidence?

Arc assignment is coordinate logic. MAPQ, secondary-alignment eligibility, anchor
thresholds, transcript-junction filtering, and artifact assessment are configurable
evidence policies.

## Result Schema

The mitochondrial remap caller uses schema
`2.1-alignment-directed-arcs-mate-aware`. The canonical quality evidence layer uses
schema `3.0-quality-evidence-multi-caller` for its shared observation and cluster
tables. Both schemas use the same directed circular coordinate convention.

Exact-deletion identifiers encode the directed left breakpoint, right breakpoint,
and deleted size under the coordinate convention defined below.

## Core Invariants

### Query order and strand normalization

Split segments are ordered by query coordinates and normalized to a
forward-reference retained adjacency. For two compatible same-strand segments:

- on the plus strand, the earlier segment ends at retained reference base `L` and
  the later segment begins at retained reference base `R`;
- on the minus strand, the later segment ends at retained reference base `L` and
  the earlier segment begins at retained reference base `R`;
- the read-supported retained adjacency is `L -> R`.

### Directed deleted arc

For mitochondrial length `N`, the inferred deleted interval is the circular forward
arc from `L` to `R`, excluding both retained breakpoint bases.

- If `R > L`, deleted bases are `L + 1` through `R - 1` and deleted size is
  `R - L - 1`.
- If `R <= L`, deleted bases are `L + 1` through `N` plus `1` through `R - 1`, and
  deleted size is `N - L + R - 1`.

`R -> L` is the complementary deletion model. It is not an interchangeable
representation of `L -> R`, even when it is shorter.

### Coordinate conversion

All calls are converted from each rotated reference to the configured standard
mitochondrial coordinate system before clustering. For a rotated reference beginning
at standard coordinate `rotation_start`:

```text
standard = ((rotated_position + rotation_start - 2) % mt_length) + 1
```

Breakpoints use 1-based retained-flanking-base coordinates. Deleted interval and
size exclude both breakpoint bases.

### Physical read identity

Alignment chains represent one physical read sequence.

- SAM `read1` and `read2` flags place paired-end mates in separate chains even when
  both mates share a query name.
- Unpaired short reads and long reads use their query name as the chain identity.
- The caller must never infer a deletion from paired-end insert distance or by
  joining an alignment from mate 1 to an alignment from mate 2.
- When both mates support the same exact deletion, support is deduplicated at the
  fragment/base-query-name level.

### Rotation handling

Normal and rotated references provide alternate coordinate origins for the same
mitochondrial molecule. After coordinate conversion:

- equivalent directed calls may be deduplicated across rotations;
- reciprocal directed calls must remain distinct;
- `rotation_support`, `rotation_count`, and `rotation_agreement` record provenance;
- a plotting coordinate must never overwrite a biological breakpoint coordinate.

### Ambiguous direction

If one fragment and unordered breakpoint pair support both `L -> R` and `R -> L`,
the evidence is marked as a reciprocal-direction conflict. The default primary
summary excludes such rows while retaining them in an audit table. Interval length
must not resolve the conflict.

## Evidence Policy

The following settings are configuration-driven and recorded in the resolved
configuration and report:

- `junctions.arc_assignment`;
- `junctions.alignment_pairing_mode`;
- `junctions.ambiguous_direction_policy`;
- minimum and maximum deletion size;
- breakpoint clustering slop and minimum read support;
- minimum anchor length;
- minimum MAPQ;
- minimum aligned fraction and maximum soft-clipped fraction;
- maximum query overlap and gap;
- eligibility of secondary and supplementary alignments;
- expected transcript-junction filtering.

The default evaluates all compatible segment pairs within one physical
read or mate. Supplementary records and SA-tagged primary or secondary records are
eligible by default; ordinary single-segment records cannot form a split chain and
are not retained by the caller.
Adjacent-only pairing and secondary-excluded calling are sensitivity modes and must
be labelled as such when used.

Every accepted read-level row retains enough provenance to audit the call:

- physical fragment and alignment-chain identifiers;
- stored query spans;
- raw and converted reference spans;
- CIGAR, strand, and SAM flags;
- primary, supplementary, and secondary status;
- MAPQ, alignment score, edit distance, `SA`, and minimap2 tags when available;
- rotation name and offset;
- directed and complementary interval properties.

## Evidence Sources And Canonical Observations

Minimap2 mitochondrial-remap split alignments are an evidence source for every
configured read type. Short-read RNA datasets can additionally enable STAR
mitochondrial-to-mitochondrial chimeric records with
`quality.short_read_rna_dual_caller.enabled`.

Each caller applies its own alignment and query-geometry filters before accepted
rows enter one canonical observation table. STAR records are expressed directly in
the configured mitochondrial coordinate system. Minimap2 records are converted from
normal or rotated mitochondrial references before consolidation.

A canonical observation represents one physical read for Nanopore or single-end
Illumina and one sequenced fragment for paired-end Illumina. The same sample,
physical observation, and directed breakpoint model detected by multiple callers or
rotations is counted once. Caller and rotation agreement are retained as provenance
rather than additional molecule support.

Canonical observations record, where available:

- caller-specific and combined evidence status;
- read or fragment identity and library layout;
- strand, CIGARs, query spans, anchors, aligned lengths, clipping, gap, overlap,
  query coverage, and distance from the junction boundary to a read end;
- segment MAPQ, alignment score, edit distance, and primary, secondary, or
  supplementary status;
- raw caller coordinates, canonical coordinates, and alignment-pattern identity;
- rotation support, direction ambiguity, transcript compatibility, nuclear
  competition status, and quality flags;
- paired-end collapse and mate-context availability.

Unavailable caller metrics are empty or carry an explicit availability status; they
are not interpreted as zero. Mate-placement fields use
`not_available_from_retained_intermediates` when the retained alignment products do
not contain the records required to calculate them.

## Canonical Deletion Objects

A canonical exact-deletion object includes at least:

- exact deletion ID and unordered diagnostic breakpoint-pair ID;
- left and right retained breakpoints;
- deleted size, interval, and origin-wrapping status;
- complementary size and wrapping status;
- arc-assignment and direction status;
- read support and normalized support;
- affected, fully removed, and partially overlapped features;
- flanking or nearest features;
- size class and configured known-deletion match;
- sample, group, and rotation metadata;
- evidence tier and quality flags;
- STAR, minimap2, both-caller, and combined distinct-observation support;
- within-sample replication, supporting-sample count, alignment-pattern diversity,
  breakpoint dispersion, anchor, error, coverage, chain-complexity, and
  multiple-hypothesis summaries;
- local remap reference-spanning support when minimap2 evidence provides a
  like-for-like split-support numerator.

Exact deletions remain coordinate-level events. Affected-feature categories are
derived labels based on overlap with the directed deleted interval and must not
replace exact-deletion identity.

## Clustering And Deduplication

Clustering uses circular breakpoint distance and respects direction.

- `L -> R` clusters only with compatible `L -> R` evidence within configured slop.
- `L -> R` does not cluster with `R -> L`.
- Representative deleted size and interval are recomputed from representative
  breakpoints; coordinates and size are not summarized independently.
- Per-read coordinates are preserved separately from cluster representatives.
- Support is deduplicated by sample and physical fragment within an exact cluster.

## Annotation And Reporting

Feature effects, control-region overlap, size classes, configured-target matches,
rainfall positions, interval midpoints, and origin markers are derived from the
directed interval.

Each report records:

- result schema, workflow configuration, reference identity, and sample metadata;
- read technology, molecule type, and library strategy;
- first-pass selection and mitochondrial-remap settings;
- arc assignment, alignment-chain identity, pairing mode, and ambiguity policy;
- alignment and deletion thresholds;
- normalization denominator;
- warnings for MAPQ 0 eligibility, secondary-alignment eligibility, single-rotation
  support, and reciprocal-direction conflicts;
- a plain-language explanation that split alignments support an inferred coordinate
  model rather than proving a viable biological deletion.

Reports distinguish:

- mitochondrial-evidence reads selected for remapping;
- deletion-supporting split-alignment reads;
- local reference-spanning reads;
- exact coordinate deletions;
- affected-feature categories.

Every canonical cluster receives one evidence-management tier:

- `strong`: at least the configured strong observation count plus within-sample
  replication or cross-caller corroboration of a physical observation;
- `supported`: at least the configured supported observation count without the
  additional strong-evidence condition;
- `review`: fewer than the configured supported observation count;
- `rejected`: a configured expected-transcript junction or unresolved reciprocal
  direction ambiguity.

These tiers describe evidence handling, not biological confirmation. The default
report profiles are cumulative views over stable exact-deletion IDs:

- `stringent`: strong;
- `standard`: strong and supported;
- `exploratory`: strong, supported, and review.

Each profile rebuilds its own tables, matrices, normalization, statistics,
ordinations, and plots. PCA and MDS axes are profile-specific. The report index at
`results/<dataset>/quality/report/index.html` links the profile reports and records
their retained cluster and observation counts. Shared canonical tables are under
`results/<dataset>/quality/shared/`, and profile outputs are under
`results/<dataset>/quality/profiles/<profile>/`.

A complete workflow target creates
`results/<dataset>/<dataset>_deliverables/index.html`. This self-contained selector
links packaged profile directories containing each report and its tables, matrices,
plots, and read lists, plus the shared canonical evidence and resolved configuration.

Short-read RNA reports include a gene-pair matrix and PCA only when the STAR evidence
stream is enabled. Gene-pair labels annotate or aggregate exact coordinate events;
they do not replace exact-deletion identity or create additional candidates.

## Assay-Specific Interpretation

### Nanopore

Long reads can provide long junction anchors. Base errors, homopolymers, alternative
placements, chimeric reads, ligation artifacts, and concatemers can affect evidence.

### Illumina

Short split anchors can be difficult to place uniquely around repeats and NUMTs.
Calls require breakpoint-spanning evidence and are not inferred from mate distance.
PCR duplication and library chimeras can inflate support.

### RNA

Transcript processing, polycistronic RNA structure, reverse transcription, template
switching, and spliced alignment can produce deletion-like junctions. Expected
transcript-junction filtering is configurable. RNA support does not directly measure
mtDNA heteroplasmy or genome copy fraction.

### DNA

DNA reads are closer to genome-molecule evidence but can reflect NUMTs, PCR or
ligation chimeras, mapping ambiguity, and sampling. Split-read support alone does not
establish molecule viability or heteroplasmy.

### Unknown molecule type

When molecule type is unspecified, reports state that DNA- and RNA-specific
interpretation is limited rather than selecting assumptions heuristically.

## Validation Requirements

Targeted tests for circular coordinate handling, query ordering, mate identity,
clustering, and deletion IDs cover:

- plus- and minus-strand stored BAM order;
- non-wrapping and origin-wrapping directed arcs;
- arcs longer and shorter than half the mitochondrial genome;
- reciprocal junction separation;
- normal/rotated coordinate equivalence and circular slop near the origin;
- paired mates sharing a query name;
- prohibition of cross-mate candidate generation;
- fragment-level support deduplication;
- arbitrary unconfigured discovery breakpoints;
- representative size recomputation after clustering;
- empty and single-rotation outputs.

For a prioritized candidate, validation may additionally compare the existing BAM
records with synthetic `L|R`, complementary `R|L`, and wild-type circular-junction
templates. Doubled-reference mapping is supplementary because repeated sequence can
introduce equivalent placements.

Rendered-output validation for tables, matrices, normalization, grouping, plotting,
and reports confirms that biological columns are used,
labels and legends are readable, empty panels state the real limitation, and report
text matches the resolved configuration.

## Rerun Safety

Before rerunning from existing intermediates:

1. Refresh metadata only when required, preferring local cached metadata.
2. Dry-run the exact Snakemake command.
3. Verify that the job list contains no unintended download, trimming, first-pass
   mapping, or mitochondrial-remapping rules.
4. Restrict allowed rules when upstream scratch inputs are absent but remap BAMs are
   complete.
5. Preserve configured literal sequence-search outputs unless those searches are
   intentionally being recomputed.

## Acceptance Criteria

The workflow satisfies this specification when:

1. Deleted intervals follow strand-normalized split-alignment query order rather
   than interval length.
2. Minus-strand segment roles are reversed to express the retained adjacency on
   the forward reference.
3. Paired mates cannot form one split-alignment chain.
4. Reciprocal models remain distinct or are explicitly marked ambiguous.
5. Rotation conversion and cross-rotation deduplication preserve direction.
6. Tables and reports expose the assumptions and provenance needed to audit a call.
7. Behavior remains configuration-driven across Nanopore, Illumina, RNA, and DNA.
8. Cross-caller detection of one physical observation increases provenance without
   increasing molecule support.
9. Stable exact-deletion IDs are shared across report profiles, while profile
   matrices and ordinations contain only the retained observations.
10. Report text and empty panels follow resolved assay, caller, layout, and metric
    availability settings.
