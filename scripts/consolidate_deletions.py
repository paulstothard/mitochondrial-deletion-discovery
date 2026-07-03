#!/usr/bin/env python3
from __future__ import annotations

import argparse
from statistics import median

from circular_deletions import circular_distance, deletion_id
from common import read_tsv, write_tsv


def close(row: dict, cluster: dict, slop: int) -> bool:
    return (
        row.get("species") == cluster.get("species")
        and abs(int(row["left_breakpoint"]) - int(cluster["left_breakpoint"])) <= slop
        and abs(int(row["right_breakpoint"]) - int(cluster["right_breakpoint"])) <= slop
    )


def cluster_rows(rows: list[dict], slop: int, min_support: int, mt_length: int | None = None) -> tuple[list[dict], list[dict], list[dict]]:
    clusters: list[dict] = []
    for row in rows:
        cluster = next((item for item in clusters if close(row, item, slop)), None)
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
            }
            clusters.append(cluster)
        support_key = (row["sample"], row["read_id"])
        if support_key in cluster["support_keys"]:
            continue
        cluster["rows"].append(row)
        cluster["support_keys"].add(support_key)
        cluster["samples"].add(row["sample"])
        if row.get("rotation_name"):
            cluster["rotations"].add(row["rotation_name"])

    kept = []
    id_rows = []
    all_rows = []
    for cluster in clusters:
        if len(cluster["support_keys"]) < min_support:
            continue
        left = int(round(median(int(row["left_breakpoint"]) for row in cluster["rows"])))
        right = int(round(median(int(row["right_breakpoint"]) for row in cluster["rows"])))
        size = circular_distance(left, right, int(mt_length)) if mt_length else int(round(median(int(row["deleted_size"]) for row in cluster["rows"])))
        exact_id = deletion_id(left, right, size)
        wraps_origin = "yes" if right <= left else "no"
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
            all_rows.append(out)
            id_rows.append(
                {
                    "sample": row["sample"],
                    "read_id": row["read_id"],
                    "exact_deletion_id": exact_id,
                    "junction_id": exact_id,
                    "left_breakpoint": left,
                    "right_breakpoint": right,
                    "reported_left_breakpoint": row.get("reported_left_breakpoint", ""),
                    "reported_right_breakpoint": row.get("reported_right_breakpoint", ""),
                }
            )
        kept.append(
            {
                "exact_deletion_id": exact_id,
                "junction_id": exact_id,
                "species": cluster["species"],
                "left_breakpoint": left,
                "right_breakpoint": right,
                "deleted_size": size,
                "wraps_origin": wraps_origin,
                "total_supporting_reads": len(cluster["support_keys"]),
                "samples_with_signal": len(cluster["samples"]),
                "supporting_samples": ",".join(sorted(cluster["samples"])),
                "rotation_support": ",".join(sorted(cluster["rotations"])),
                "rotation_count": len(cluster["rotations"]),
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
    all_rows, clusters, id_rows = cluster_rows(rows, args.slop, args.min_support, args.mt_length)
    read_fields = [
        "exact_deletion_id",
        "junction_id",
        "sample",
        "species",
        "read_id",
        "deletion_id",
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "read_left_breakpoint",
        "read_right_breakpoint",
        "read_deleted_size",
        "wraps_origin",
        "deleted_interval",
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
    for row in all_rows:
        if not row.get("deletion_id"):
            row["deletion_id"] = row.get("exact_deletion_id", row.get("junction_id", ""))
    write_tsv(args.all_reads, all_rows, read_fields)
    write_tsv(
        args.clusters,
        clusters,
        [
            "exact_deletion_id",
            "junction_id",
            "species",
            "left_breakpoint",
            "right_breakpoint",
            "deleted_size",
            "wraps_origin",
            "total_supporting_reads",
            "samples_with_signal",
            "supporting_samples",
            "rotation_support",
            "rotation_count",
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
            "left_breakpoint",
            "right_breakpoint",
            "reported_left_breakpoint",
            "reported_right_breakpoint",
        ],
    )


if __name__ == "__main__":
    main()
