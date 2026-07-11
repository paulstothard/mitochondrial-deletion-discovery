# Circular Coordinate And Merge Audit

> Scope limitation: this audit validates coordinate conversion, round trips,
> flanking-base conventions, and representative-coordinate consistency. Its tables
> do not independently validate directed arc assignment. Directed arc behavior is
> defined and tested by the current workflow specification and caller tests.

Dataset used for concrete examples: `human_common_deletion`. This is a targeted audit of coordinate conversion and merge behavior, not a full deletion report.

## Code Paths Checked

- Rotated reference construction: `scripts/make_rotated_mt_reference.py:13-17`.
- Rotated-to-standard coordinate conversion: `scripts/circular_deletions.py:11-15`.
- Deleted-size and interval convention: `scripts/circular_deletions.py:18-38`.
- Split-read deletion calling and storage of raw versus converted coordinates: `scripts/call_minimap2_deletions.py:51-62` and `scripts/call_minimap2_deletions.py:85-123`.
- Merge/deduplication of rotated-reference calls: `scripts/consolidate_deletions.py:11-88`.
- Plot-only origin-crossing y-coordinate for the breakpoint-pair map: `scripts/plot_deletion_results.py:1249-1252`.

## Coordinate Convention

The workflow uses 1-based mitochondrial coordinates after reading BAM alignments. `pysam.reference_start` is 0-based, so the caller adds 1. `pysam.reference_end` is 0-based exclusive, which is equivalent to the 1-based inclusive end of the aligned segment. The stored `left_breakpoint` and `right_breakpoint` are the retained flanking bases around the deleted interval, not the first and last deleted bases.

Because of that flanking-breakpoint convention, deleted size is `right - left - 1` for non-wrapping deletions and `mt_length - left + right - 1` for origin-spanning deletions. Deleted intervals are stored as 1-based closed intervals between the breakpoints, for example `left + 1` through `right - 1`. If an external notation defines breakpoints as first/last deleted bases, it will differ by one base from this workflow's flanking-base notation.

## Offset Direction

Mitochondrial genome length in this audit: `16569`. Rotation starts: `{'normal': 1, 'half': 8285}`.

A rotated reference with `rotation_start = X` begins with standard coordinate `X`. Therefore rotated position 1 converts to standard position `X`. The conversion is:

`standard = ((rotated_position + rotation_start - 2) % mt_length) + 1`

The inverse used for the round-trip checks is:

`rotated = ((standard_position - rotation_start) % mt_length) + 1`

See `position_roundtrip.tsv` for positions near the standard origin, standard genome end, and offset boundary.

## Targeted Results

- Worked examples written: `5` rows in `worked_examples.tsv`.
- Table-level coordinate checks are in `table_checks.tsv`.
- A merge check that recomputes interval properties from representative breakpoints
  is in `corrected_merge_table_checks.tsv`.
- Read-to-cluster coordinate comparison: 926 of 3040 read-level rows differ from the merged representative coordinates; see `read_vs_cluster_coordinate_differences.tsv`.

## Audit Conclusions

- Wrong offset direction: no evidence in the round-trip tests or worked examples.
- Off-by-one error: no internal off-by-one inconsistency was found under the workflow's flanking-breakpoint convention. The convention itself must be stated clearly because it differs from deletion-size formulas that treat start/end as deleted-base coordinates.
- Sorting start/end incorrectly: no evidence that converted coordinates are blindly sorted. Origin-spanning calls preserve `left_breakpoint > right_breakpoint`.
- Negative deletion lengths: no negative sizes found in the checked tables.
- Deleted-size consistency: the source cluster table has `65` size mismatches because
  left breakpoint, right breakpoint, and size were summarized independently. The
  representative-breakpoint recomputation has `0` size mismatches.
- Origin-spanning classification: assigned from converted standard coordinates using `right <= left`, after conversion.
- Merging before coordinate conversion: no evidence. Merge inputs already contain standard coordinates produced by the caller.
- Double-counting support across rotations: the merge code deduplicates by `(sample, read_id)` within a breakpoint cluster. The same read can still support different clusters if it produces distinct breakpoint pairs outside the configured slop.
- Overwriting true coordinates with plotting coordinates: no evidence. The breakpoint-pair support map creates `adjusted_right_breakpoint` as a plotting-only column and does not overwrite `right_breakpoint`.
- Creating origin-spanning calls due to incorrect wrapping: no evidence from these checks.

## Important Reporting/Plotting Note

The merge step writes representative merged coordinates to the standard `left_breakpoint`, `right_breakpoint`, `deleted_size`, and `wraps_origin` columns in `all_samples.filtered_junction_reads.tsv`. Per-read converted coordinates are retained separately as `read_left_breakpoint`, `read_right_breakpoint`, and `read_deleted_size`. With nonzero `breakpoint_slop_bp`, read-level rows can differ slightly from their cluster representative; this is normal provenance, but downstream exact-deletion plots and tables should use the representative columns. Report plots also rejoin the merged cluster table and recompute `deleted_size` from representative circular breakpoints before plotting.
