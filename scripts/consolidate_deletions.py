#!/usr/bin/env python3
from __future__ import annotations

import argparse
from statistics import median

from circular_deletions import (
    breakpoint_pair_id,
    circular_position_distance,
    deletion_id,
    directed_breakpoints,
    normalize_pos,
)
from common import read_tsv, write_tsv


def close(row: dict, cluster: dict, slop: int, mt_length: int) -> bool:
    return (
        row.get("species") == cluster.get("species")
        and circular_position_distance(int(row["left_breakpoint"]), int(cluster["left_breakpoint"]), mt_length) <= slop
        and circular_position_distance(int(row["right_breakpoint"]), int(cluster["right_breakpoint"]), mt_length) <= slop
    )


def circular_median(values: list[int], mt_length: int) -> int:
    anchor = int(values[0])
    unwrapped = []
    for value in values:
        delta = ((int(value) - anchor + mt_length // 2) % mt_length) - mt_length // 2
        unwrapped.append(anchor + delta)
    return normalize_pos(int(round(median(unwrapped))), mt_length)


def split_direction_conflicts(
    rows: list[dict],
    mt_length: int | None = None,
    slop: int = 0,
) -> tuple[list[dict], list[dict]]:
    grouped: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(rows):
        grouped.setdefault((row.get("sample", ""), row.get("read_id", "")), []).append(index)
    ambiguous_indices = set()
    for indices in grouped.values():
        for offset, first_index in enumerate(indices):
            first = rows[first_index]
            for second_index in indices[offset + 1 :]:
                second = rows[second_index]
                if mt_length:
                    reciprocal = (
                        circular_position_distance(first["left_breakpoint"], second["right_breakpoint"], mt_length) <= slop
                        and circular_position_distance(first["right_breakpoint"], second["left_breakpoint"], mt_length) <= slop
                    )
                else:
                    reciprocal = (
                        int(first["left_breakpoint"]) == int(second["right_breakpoint"])
                        and int(first["right_breakpoint"]) == int(second["left_breakpoint"])
                    )
                if reciprocal:
                    ambiguous_indices.update([first_index, second_index])
    accepted = []
    ambiguous = []
    for index, row in enumerate(rows):
        out = dict(row)
        pair_id = out.get("breakpoint_pair_id") or breakpoint_pair_id(out["left_breakpoint"], out["right_breakpoint"])
        out["breakpoint_pair_id"] = pair_id
        if index in ambiguous_indices:
            out["direction_status"] = "ambiguous_across_rotations"
            ambiguous.append(out)
        else:
            accepted.append(out)
    return accepted, ambiguous


def cluster_rows(
    rows: list[dict],
    slop: int,
    min_support: int,
    mt_length: int | None = None,
    result_schema_version: str = "2.0-alignment-directed-arcs",
) -> tuple[list[dict], list[dict], list[dict]]:
    if not mt_length:
        raise ValueError("mt_length is required for directed circular clustering")
    clusters: list[dict] = []
    bin_size = max(1, int(slop))
    bin_count = max(1, (int(mt_length) + bin_size - 1) // bin_size)
    cluster_bins: dict[tuple[str, int, int], list[int]] = {}

    def bin_index(pos: int) -> int:
        return ((int(pos) - 1) // bin_size) % bin_count

    for row in rows:
        left_bin = bin_index(int(row["left_breakpoint"]))
        right_bin = bin_index(int(row["right_breakpoint"]))
        # The final circular bin can be shorter than bin_size, so two neighboring
        # bins are needed to retain exact slop behavior across coordinate 1.
        radius = 0 if slop <= 0 else 2
        candidate_indices = set()
        for left_delta in range(-radius, radius + 1):
            for right_delta in range(-radius, radius + 1):
                key = (
                    row.get("species", ""),
                    (left_bin + left_delta) % bin_count,
                    (right_bin + right_delta) % bin_count,
                )
                candidate_indices.update(cluster_bins.get(key, []))
        cluster_index = next(
            (index for index in sorted(candidate_indices) if close(row, clusters[index], slop, int(mt_length))),
            None,
        )
        cluster = clusters[cluster_index] if cluster_index is not None else None
        if cluster is None:
            cluster = {
                "species": row["species"],
                "left_breakpoint": int(row["left_breakpoint"]),
                "right_breakpoint": int(row["right_breakpoint"]),
                "deleted_size": int(row["deleted_size"]),
                "rows": [],
                "support_keys": set(),
                "samples": set(),
                "rotations": set(),
                "arc_assignment_methods": set(),
                "direction_statuses": set(),
            }
            clusters.append(cluster)
            cluster_index = len(clusters) - 1
            cluster_bins.setdefault((row.get("species", ""), left_bin, right_bin), []).append(cluster_index)
        if row.get("rotation_name"):
            cluster["rotations"].add(row["rotation_name"])
        if row.get("arc_assignment_method"):
            cluster["arc_assignment_methods"].add(row["arc_assignment_method"])
        if row.get("direction_status"):
            cluster["direction_statuses"].add(row["direction_status"])
        support_key = (row["sample"], row["read_id"])
        if support_key in cluster["support_keys"]:
            continue
        cluster["rows"].append(row)
        cluster["support_keys"].add(support_key)
        cluster["samples"].add(row["sample"])

    kept = []
    id_rows = []
    all_rows = []
    for cluster in clusters:
        if len(cluster["support_keys"]) < min_support:
            continue
        left = circular_median([int(row["left_breakpoint"]) for row in cluster["rows"]], int(mt_length))
        right = circular_median([int(row["right_breakpoint"]) for row in cluster["rows"]], int(mt_length))
        directed = directed_breakpoints(left, right, int(mt_length))
        size = int(directed["deleted_size"])
        exact_id = deletion_id(left, right, size)
        pair_id = breakpoint_pair_id(left, right)
        wraps_origin = directed["wraps_origin"]
        direction_status = ";".join(sorted(cluster["direction_statuses"])) or "directed"
        arc_assignment_method = ";".join(sorted(cluster["arc_assignment_methods"])) or "alignment_directed"
        rotation_agreement = "multiple_rotations" if len(cluster["rotations"]) > 1 else "single_rotation"
        for row in cluster["rows"]:
            out = dict(row)
            out["read_left_breakpoint"] = row.get("left_breakpoint", "")
            out["read_right_breakpoint"] = row.get("right_breakpoint", "")
            out["read_deleted_size"] = row.get("deleted_size", "")
            out["exact_deletion_id"] = exact_id
            out["junction_id"] = exact_id
            out["deletion_id"] = exact_id
            out["left_breakpoint"] = left
            out["right_breakpoint"] = right
            out["deleted_size"] = size
            out["wraps_origin"] = wraps_origin
            out["deleted_interval"] = directed["deleted_interval"]
            out["complement_deleted_size"] = directed["complement_deleted_size"]
            out["complement_wraps_origin"] = directed["complement_wraps_origin"]
            out["breakpoint_pair_id"] = pair_id
            out["arc_assignment_method"] = arc_assignment_method
            out["direction_status"] = direction_status
            out["rotation_agreement"] = rotation_agreement
            out["result_schema_version"] = result_schema_version
            all_rows.append(out)
            id_rows.append(
                {
                    "sample": row["sample"],
                    "read_id": row["read_id"],
                    "exact_deletion_id": exact_id,
                    "junction_id": exact_id,
                    "breakpoint_pair_id": pair_id,
                    "left_breakpoint": left,
                    "right_breakpoint": right,
                    "reported_left_breakpoint": row.get("reported_left_breakpoint", ""),
                    "reported_right_breakpoint": row.get("reported_right_breakpoint", ""),
                    "direction_status": direction_status,
                    "result_schema_version": result_schema_version,
                }
            )
        kept.append(
            {
                "exact_deletion_id": exact_id,
                "junction_id": exact_id,
                "breakpoint_pair_id": pair_id,
                "species": cluster["species"],
                "left_breakpoint": left,
                "right_breakpoint": right,
                "deleted_size": size,
                "wraps_origin": wraps_origin,
                "deleted_interval": directed["deleted_interval"],
                "complement_deleted_size": directed["complement_deleted_size"],
                "complement_wraps_origin": directed["complement_wraps_origin"],
                "arc_assignment_method": arc_assignment_method,
                "direction_status": direction_status,
                "result_schema_version": result_schema_version,
                "total_supporting_reads": len(cluster["support_keys"]),
                "samples_with_signal": len(cluster["samples"]),
                "supporting_samples": ",".join(sorted(cluster["samples"])),
                "rotation_support": ",".join(sorted(cluster["rotations"])),
                "rotation_count": len(cluster["rotations"]),
                "rotation_agreement": rotation_agreement,
                "clustered_breakpoint_slop_bp": slop,
            }
        )
    return all_rows, kept, id_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--slop", type=int, required=True)
    parser.add_argument("--min-support", type=int, required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--all-reads", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--id-map", required=True)
    parser.add_argument("--ambiguous-reads", required=True)
    parser.add_argument("--ambiguous-direction-policy", choices=["exclude", "include"], default="exclude")
    parser.add_argument("--result-schema-version", default="2.0-alignment-directed-arcs")
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    rows = []
    for path in args.inputs:
        rows.extend(read_tsv(path))
    rows.sort(
        key=lambda row: (
            row.get("species", ""),
            int(row.get("left_breakpoint", 0) or 0),
            int(row.get("right_breakpoint", 0) or 0),
            row.get("sample", ""),
            row.get("read_id", ""),
            row.get("rotation_name", ""),
        )
    )
    accepted_rows, ambiguous_rows = split_direction_conflicts(rows, mt_length=args.mt_length, slop=args.slop)
    cluster_inputs = accepted_rows + ambiguous_rows if args.ambiguous_direction_policy == "include" else accepted_rows
    all_rows, clusters, id_rows = cluster_rows(
        cluster_inputs,
        args.slop,
        args.min_support,
        args.mt_length,
        result_schema_version=args.result_schema_version,
    )
    read_fields = [
        "exact_deletion_id",
        "junction_id",
        "sample",
        "species",
        "read_id",
        "deletion_id",
        "breakpoint_pair_id",
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "read_left_breakpoint",
        "read_right_breakpoint",
        "read_deleted_size",
        "wraps_origin",
        "deleted_interval",
        "complement_deleted_size",
        "complement_wraps_origin",
        "arc_assignment_method",
        "direction_status",
        "rotation_agreement",
        "result_schema_version",
        "reported_left_breakpoint",
        "reported_right_breakpoint",
        "reported_deleted_size",
        "canonical_orientation",
        "left_anchor_length",
        "right_anchor_length",
        "min_anchor_length",
        "strand",
        "source",
        "rotation_name",
        "rotation_start",
        "raw_left_breakpoint",
        "raw_right_breakpoint",
        "left_mapq",
        "right_mapq",
        "min_mapq",
        "query_overlap_bp",
        "query_gap_bp",
    ]
    provenance_fields = []
    for row in all_rows:
        for key in row:
            if key not in read_fields and key not in provenance_fields:
                provenance_fields.append(key)
    read_fields.extend(provenance_fields)
    for row in all_rows:
        if not row.get("deletion_id"):
            row["deletion_id"] = row.get("exact_deletion_id", row.get("junction_id", ""))
    write_tsv(args.all_reads, all_rows, read_fields)
    ambiguous_fields = list(ambiguous_rows[0]) if ambiguous_rows else read_fields
    write_tsv(args.ambiguous_reads, ambiguous_rows, ambiguous_fields)
    write_tsv(
        args.clusters,
        clusters,
        [
            "exact_deletion_id",
            "junction_id",
            "breakpoint_pair_id",
            "species",
            "left_breakpoint",
            "right_breakpoint",
            "deleted_size",
            "wraps_origin",
            "deleted_interval",
            "complement_deleted_size",
            "complement_wraps_origin",
            "arc_assignment_method",
            "direction_status",
            "result_schema_version",
            "total_supporting_reads",
            "samples_with_signal",
            "supporting_samples",
            "rotation_support",
            "rotation_count",
            "rotation_agreement",
            "clustered_breakpoint_slop_bp",
        ],
    )
    write_tsv(
        args.id_map,
        id_rows,
        [
            "sample",
            "read_id",
            "exact_deletion_id",
            "junction_id",
            "breakpoint_pair_id",
            "left_breakpoint",
            "right_breakpoint",
            "reported_left_breakpoint",
            "reported_right_breakpoint",
            "direction_status",
            "result_schema_version",
        ],
    )


if __name__ == "__main__":
    main()
