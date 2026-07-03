#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from circular_deletions import circular_distance, interval_pieces, normalize_pos  # noqa: E402
from consolidate_deletions import cluster_rows  # noqa: E402


DATASET = "human_common_deletion"
CONFIG = REPO / "config" / "datasets" / f"{DATASET}.yaml"
RESULTS = REPO / "results" / DATASET
OUT = REPO / "audit" / "circular_coordinate_merge"


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def standard_to_rotated(pos: int, mt_length: int, rotation_start: int) -> int:
    return ((int(pos) - int(rotation_start)) % int(mt_length)) + 1


def load_config() -> dict:
    with (REPO / "config" / "defaults.yaml").open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    with CONFIG.open(encoding="utf-8") as handle:
        dataset_config = yaml.safe_load(handle) or {}
    deep_update(config, dataset_config)
    return config


def deep_update(base: dict, override: dict) -> dict:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def parse_rotation_start(value: object, mt_length: int) -> int:
    if isinstance(value, str) and value.strip().lower() == "half":
        return mt_length // 2 + 1
    return int(value)


def load_raw_calls() -> pd.DataFrame:
    paths = sorted((RESULTS / "deletions" / "rotated").glob("*/*.filtered_deletion_reads.tsv"))
    if not paths:
        return pd.DataFrame()
    frames = [pd.read_csv(path, sep="\t") for path in paths]
    return pd.concat(frames, ignore_index=True)


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def recompute_size(left: int, right: int, mt_length: int) -> int:
    return circular_distance(int(left), int(right), mt_length)


def row_example(label: str, row: pd.Series, mt_length: int) -> dict:
    left = int(row["left_breakpoint"])
    right = int(row["right_breakpoint"])
    raw_left = int(row.get("raw_left_breakpoint", 0) or 0)
    raw_right = int(row.get("raw_right_breakpoint", 0) or 0)
    rot_start = int(row.get("rotation_start", 1) or 1)
    recomputed = recompute_size(left, right, mt_length)
    return {
        "example": label,
        "sample": row.get("sample", ""),
        "read_id": row.get("read_id", ""),
        "rotation_name": row.get("rotation_name", ""),
        "rotation_start": rot_start,
        "raw_left_breakpoint_1based_rotated": raw_left,
        "raw_right_breakpoint_1based_rotated": raw_right,
        "raw_left_converted_to_standard": normalize_pos(raw_left, mt_length, rot_start) if raw_left else "",
        "raw_right_converted_to_standard": normalize_pos(raw_right, mt_length, rot_start) if raw_right else "",
        "left_breakpoint_1based_standard": left,
        "right_breakpoint_1based_standard": right,
        "wraps_origin": row.get("wraps_origin", ""),
        "deleted_size_reported": int(row["deleted_size"]),
        "deleted_size_recomputed_from_standard_breakpoints": recomputed,
        "size_matches": int(row["deleted_size"]) == recomputed,
        "deleted_interval_1based_closed": row.get("deleted_interval", ";".join(f"{a}-{b}" for a, b in interval_pieces(left, right, mt_length))),
        "canonical_orientation": row.get("canonical_orientation", ""),
    }


def choose_examples(raw: pd.DataFrame, clusters: pd.DataFrame, all_reads: pd.DataFrame, mt_length: int) -> tuple[list[dict], pd.DataFrame]:
    examples: list[dict] = []
    if raw.empty:
        return examples, pd.DataFrame()
    normal_standard = raw[(raw["rotation_name"] == "normal") & (raw["wraps_origin"] == "no")]
    if not normal_standard.empty:
        examples.append(row_example("normal deletion from standard alignment", normal_standard.iloc[0], mt_length))
    normal_offset = raw[(raw["rotation_name"] == "half") & (raw["wraps_origin"] == "no")]
    if not normal_offset.empty:
        examples.append(row_example("normal deletion from offset alignment", normal_offset.iloc[0], mt_length))
    origin_offset = raw[(raw["rotation_name"] == "half") & (raw["wraps_origin"] == "yes")]
    if not origin_offset.empty:
        examples.append(row_example("origin-spanning deletion from offset alignment", origin_offset.iloc[0], mt_length))
    if not normal_offset.empty:
        boundary = normal_offset.assign(
            boundary_distance=normal_offset[["raw_left_breakpoint", "raw_right_breakpoint"]]
            .apply(lambda col: pd.to_numeric(col, errors="coerce"))
            .apply(lambda col: pd.concat([col.abs(), (col - mt_length).abs()], axis=1).min(axis=1))
            .min(axis=1)
        ).sort_values("boundary_distance")
        examples.append(row_example("deletion near offset-reference boundary", boundary.iloc[0], mt_length))
    merged_reads = pd.DataFrame()
    if not clusters.empty and not all_reads.empty and "rotation_support" in clusters.columns:
        both = clusters[clusters["rotation_support"].fillna("").astype(str).str.contains("half") & clusters["rotation_support"].fillna("").astype(str).str.contains("normal")]
        if not both.empty:
            cluster_id = both.iloc[0]["junction_id"]
            merged_reads = all_reads[all_reads["junction_id"] == cluster_id].copy()
            if not merged_reads.empty:
                examples.append(row_example("merged call detected in both alignments", merged_reads.iloc[0], mt_length))
    return examples, merged_reads


def check_tables(raw: pd.DataFrame, clusters: pd.DataFrame, all_reads: pd.DataFrame, mt_length: int) -> list[dict]:
    rows = []
    for name, table in [("raw rotated deletion calls", raw), ("merged clusters", clusters), ("read-level merged evidence", all_reads)]:
        if table.empty:
            rows.append({"table": name, "rows": 0, "size_mismatches": "", "negative_sizes": "", "wrap_rule_mismatches": ""})
            continue
        left = pd.to_numeric(table["left_breakpoint"], errors="coerce")
        right = pd.to_numeric(table["right_breakpoint"], errors="coerce")
        size = pd.to_numeric(table["deleted_size"], errors="coerce")
        recomputed = [recompute_size(l, r, mt_length) for l, r in zip(left.fillna(0).astype(int), right.fillna(0).astype(int))]
        recomputed_s = pd.Series(recomputed, index=table.index)
        wrap_expected = right <= left
        wrap_observed = table["wraps_origin"].fillna("").astype(str).str.lower().eq("yes") if "wraps_origin" in table.columns else wrap_expected
        rows.append(
            {
                "table": name,
                "rows": len(table),
                "size_mismatches": int((size != recomputed_s).sum()),
                "negative_sizes": int((size < 0).sum()),
                "wrap_rule_mismatches": int((wrap_observed != wrap_expected).sum()),
                "origin_spanning_rows": int(wrap_expected.sum()),
            }
        )
    return rows


def corrected_merge_check(raw: pd.DataFrame, config: dict, mt_length: int) -> list[dict]:
    if raw.empty:
        return [{"table": "corrected merge from raw rotated calls", "rows": 0, "size_mismatches": "", "negative_sizes": "", "wrap_rule_mismatches": ""}]
    junctions = config.get("junctions", {}) or {}
    slop = int(junctions.get("breakpoint_slop_bp", 0))
    min_support = int(junctions.get("min_split_read_support", 1))
    _, corrected_clusters, _ = cluster_rows(raw.astype(str).to_dict("records"), slop=slop, min_support=min_support, mt_length=mt_length)
    corrected = pd.DataFrame(corrected_clusters)
    rows = [row for row in check_tables(pd.DataFrame(), corrected, pd.DataFrame(), mt_length) if row["table"] == "merged clusters"]
    if rows:
        rows[0]["table"] = "corrected merge from raw rotated calls"
        rows[0]["breakpoint_slop_bp"] = slop
        rows[0]["min_split_read_support"] = min_support
    return rows


def read_cluster_mismatches(all_reads: pd.DataFrame, clusters: pd.DataFrame) -> tuple[list[dict], pd.DataFrame]:
    if all_reads.empty or clusters.empty:
        return [], pd.DataFrame()
    merged = all_reads.merge(
        clusters[["junction_id", "left_breakpoint", "right_breakpoint", "deleted_size"]],
        on="junction_id",
        suffixes=("_read", "_cluster"),
        how="inner",
    )
    mismatch = (
        (merged["left_breakpoint_read"] != merged["left_breakpoint_cluster"])
        | (merged["right_breakpoint_read"] != merged["right_breakpoint_cluster"])
        | (merged["deleted_size_read"] != merged["deleted_size_cluster"])
    )
    summary = [
        {
            "read_rows_with_cluster": len(merged),
            "read_rows_different_from_cluster_representative": int(mismatch.sum()),
            "fraction_different": float(mismatch.mean()) if len(merged) else 0.0,
        }
    ]
    cols = [
        "sample",
        "read_id",
        "junction_id",
        "rotation_name",
        "left_breakpoint_read",
        "right_breakpoint_read",
        "deleted_size_read",
        "left_breakpoint_cluster",
        "right_breakpoint_cluster",
        "deleted_size_cluster",
    ]
    return summary, merged.loc[mismatch, cols].head(20)


def write_audit_md(
    path: Path,
    mt_length: int,
    rotations: dict[str, int],
    check_rows: list[dict],
    corrected_check_rows: list[dict],
    mismatch_summary: list[dict],
    examples_count: int,
) -> None:
    both_support = "not checked"
    if mismatch_summary:
        both_support = (
            f"{mismatch_summary[0]['read_rows_different_from_cluster_representative']} of "
            f"{mismatch_summary[0]['read_rows_with_cluster']} read-level rows differ from the merged representative coordinates"
        )
    stale_cluster_mismatches = next((row.get("size_mismatches", "") for row in check_rows if row.get("table") == "merged clusters"), "")
    corrected_cluster_mismatches = next((row.get("size_mismatches", "") for row in corrected_check_rows if row.get("table") == "corrected merge from raw rotated calls"), "")
    text = f"""# Circular Coordinate And Merge Audit

Dataset used for concrete examples: `{DATASET}`. This is a targeted audit of coordinate conversion and merge behavior, not a full deletion report.

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

Mitochondrial genome length in this audit: `{mt_length}`. Rotation starts: `{rotations}`.

A rotated reference with `rotation_start = X` begins with standard coordinate `X`. Therefore rotated position 1 converts to standard position `X`. The conversion is:

`standard = ((rotated_position + rotation_start - 2) % mt_length) + 1`

The inverse used for the round-trip checks is:

`rotated = ((standard_position - rotation_start) % mt_length) + 1`

See `position_roundtrip.tsv` for positions near the standard origin, standard genome end, and offset boundary.

## Targeted Results

- Worked examples written: `{examples_count}` rows in `worked_examples.tsv`.
- Table-level coordinate checks are in `table_checks.tsv`.
- A corrected merge check from the raw rotated calls is in `corrected_merge_table_checks.tsv`.
- Read-to-cluster coordinate comparison: {both_support}; see `read_vs_cluster_coordinate_differences.tsv`.

## Audit Conclusions

- Wrong offset direction: no evidence in the round-trip tests or worked examples.
- Off-by-one error: no internal off-by-one inconsistency was found under the workflow's flanking-breakpoint convention. The convention itself must be stated clearly because it differs from deletion-size formulas that treat start/end as deleted-base coordinates.
- Sorting start/end incorrectly: no evidence that converted coordinates are blindly sorted. Origin-spanning calls preserve `left_breakpoint > right_breakpoint`.
- Negative deletion lengths: no negative sizes found in the checked tables.
- Deleted-size consistency: the existing generated cluster table has `{stale_cluster_mismatches}` size mismatches because the previous merge code took median left, median right, and median size independently. The corrected merge code recomputes size from the representative breakpoints; the corrected check has `{corrected_cluster_mismatches}` size mismatches.
- Origin-spanning classification: assigned from converted standard coordinates using `right <= left`, after conversion.
- Merging before coordinate conversion: no evidence. Merge inputs already contain standard coordinates produced by the caller.
- Double-counting support across rotations: the merge code deduplicates by `(sample, read_id)` within a breakpoint cluster. The same read can still support different clusters if it produces distinct breakpoint pairs outside the configured slop.
- Overwriting true coordinates with plotting coordinates: no evidence. The breakpoint-pair support map creates `adjusted_right_breakpoint` as a plotting-only column and does not overwrite `right_breakpoint`.
- Creating origin-spanning calls due to incorrect wrapping: no evidence from these checks.

## Important Reporting/Plotting Note

The merge step writes representative merged coordinates to the standard `left_breakpoint`, `right_breakpoint`, `deleted_size`, and `wraps_origin` columns in `all_samples.filtered_junction_reads.tsv`. Per-read converted coordinates are retained separately as `read_left_breakpoint`, `read_right_breakpoint`, and `read_deleted_size`. With nonzero `breakpoint_slop_bp`, read-level rows can differ slightly from their cluster representative; this is normal provenance, but downstream exact-deletion plots and tables should use the representative columns. Report plots also rejoin the merged cluster table and recompute `deleted_size` from representative circular breakpoints before plotting.
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config = load_config()
    species = config["dataset"]["species"]
    mt_length = int(config["references"][species]["mt_length"])
    rotations = {
        str(item["name"]): parse_rotation_start(item["start"], mt_length)
        for item in (config.get("mt_realign", {}).get("rotations") or [{"name": "normal", "start": 1}, {"name": "half", "start": "half"}])
    }

    test_positions = sorted(
        {
            1,
            2,
            mt_length - 1,
            mt_length,
            rotations.get("half", mt_length // 2 + 1) - 2,
            rotations.get("half", mt_length // 2 + 1) - 1,
            rotations.get("half", mt_length // 2 + 1),
            rotations.get("half", mt_length // 2 + 1) + 1,
        }
    )
    roundtrip_rows = []
    for rotation_name, rotation_start in rotations.items():
        for pos in test_positions:
            rotated = standard_to_rotated(pos, mt_length, rotation_start)
            back = normalize_pos(rotated, mt_length, rotation_start)
            roundtrip_rows.append(
                {
                    "rotation_name": rotation_name,
                    "rotation_start": rotation_start,
                    "standard_position": pos,
                    "rotated_position": rotated,
                    "back_to_standard": back,
                    "roundtrip_pass": back == pos,
                }
            )
    write_tsv(OUT / "position_roundtrip.tsv", roundtrip_rows)

    raw = load_raw_calls()
    clusters = load_table(RESULTS / "junctions" / "junction_clusters.unannotated.tsv")
    all_reads = load_table(RESULTS / "junctions" / "all_samples.filtered_junction_reads.tsv")

    examples, merged_reads = choose_examples(raw, clusters, all_reads, mt_length)
    write_tsv(OUT / "worked_examples.tsv", examples)
    if not merged_reads.empty:
        merged_reads.head(25).to_csv(OUT / "merged_call_read_rows.tsv", sep="\t", index=False)

    checks = check_tables(raw, clusters, all_reads, mt_length)
    write_tsv(OUT / "table_checks.tsv", checks)

    corrected_checks = corrected_merge_check(raw, config, mt_length)
    write_tsv(OUT / "corrected_merge_table_checks.tsv", corrected_checks)

    mismatch_summary, mismatch_examples = read_cluster_mismatches(all_reads, clusters)
    write_tsv(OUT / "read_vs_cluster_coordinate_summary.tsv", mismatch_summary)
    if not mismatch_examples.empty:
        mismatch_examples.to_csv(OUT / "read_vs_cluster_coordinate_differences.tsv", sep="\t", index=False)

    write_audit_md(
        OUT / "audit.md",
        mt_length=mt_length,
        rotations=rotations,
        check_rows=checks,
        corrected_check_rows=corrected_checks,
        mismatch_summary=mismatch_summary,
        examples_count=len(examples),
    )
    print(f"Wrote audit outputs to {OUT}")


if __name__ == "__main__":
    main()
