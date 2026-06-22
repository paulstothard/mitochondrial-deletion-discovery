#!/usr/bin/env python3
from __future__ import annotations

import argparse
from common import read_tsv, write_tsv


def circular_deletion_size(left: int, right: int, mt_length: int) -> int:
    if right > left:
        return right - left - 1
    return mt_length - left + right - 1


def canonical_junction(row: dict, mt_length: int) -> dict:
    left = int(row["left_breakpoint"])
    right = int(row["right_breakpoint"])
    forward_size = circular_deletion_size(left, right, mt_length)
    reverse_size = circular_deletion_size(right, left, mt_length)
    canonical = dict(row)
    canonical["reported_left_breakpoint"] = left
    canonical["reported_right_breakpoint"] = right
    canonical["reported_deleted_size"] = int(row["deleted_size"])
    if reverse_size < forward_size or (reverse_size == forward_size and (right, left) < (left, right)):
        canonical["left_breakpoint"] = right
        canonical["right_breakpoint"] = left
        canonical["deleted_size"] = reverse_size
        canonical["canonical_orientation"] = "reversed_to_shorter_interval"
    else:
        canonical["left_breakpoint"] = left
        canonical["right_breakpoint"] = right
        canonical["deleted_size"] = forward_size
        canonical["canonical_orientation"] = "reported"
    return canonical


def close(row: dict, cluster: dict, slop: int) -> bool:
    return (
        row["species"] == cluster["species"]
        and abs(int(row["left_breakpoint"]) - int(cluster["left_breakpoint"])) <= slop
        and abs(int(row["right_breakpoint"]) - int(cluster["right_breakpoint"])) <= slop
    )


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

    reads = []
    for path in args.inputs:
        reads.extend(canonical_junction(row, args.mt_length) for row in read_tsv(path))
    reads.sort(key=lambda row: (row.get("species", ""), int(row["left_breakpoint"]), int(row["right_breakpoint"]), row["sample"], row["read_id"]))
    clusters = []
    id_rows = []
    id_seen = set()
    for row in reads:
        cluster = next((item for item in clusters if close(row, item, args.slop)), None)
        if cluster is None:
            cluster = {
                "junction_id": f"mtDelJunc_{len(clusters) + 1:06d}",
                "species": row["species"],
                "left_breakpoint": int(row["left_breakpoint"]),
                "right_breakpoint": int(row["right_breakpoint"]),
                "deleted_size": int(row["deleted_size"]),
                "supporting_read_keys": set(),
                "supporting_samples": set(),
                "rotations": set(),
                "reported_orientations": set(),
            }
            clusters.append(cluster)
        support_key = (row["sample"], row["read_id"])
        cluster["supporting_read_keys"].add(support_key)
        cluster["supporting_samples"].add(row["sample"])
        if row.get("rotation_name"):
            cluster["rotations"].add(row["rotation_name"])
        if row.get("canonical_orientation"):
            cluster["reported_orientations"].add(row["canonical_orientation"])
        id_key = (row["sample"], row["read_id"], cluster["junction_id"])
        if id_key not in id_seen:
            id_rows.append(
                {
                    "sample": row["sample"],
                    "read_id": row["read_id"],
                    "junction_id": cluster["junction_id"],
                    "left_breakpoint": row["left_breakpoint"],
                    "right_breakpoint": row["right_breakpoint"],
                    "reported_left_breakpoint": row["reported_left_breakpoint"],
                    "reported_right_breakpoint": row["reported_right_breakpoint"],
                }
            )
            id_seen.add(id_key)
        row["junction_id"] = cluster["junction_id"]

    cluster_rows = []
    for cluster in clusters:
        supporting_reads = len(cluster["supporting_read_keys"])
        if supporting_reads < args.min_support:
            continue
        cluster_rows.append(
            {
                "junction_id": cluster["junction_id"],
                "species": cluster["species"],
                "left_breakpoint": cluster["left_breakpoint"],
                "right_breakpoint": cluster["right_breakpoint"],
                "deleted_size": cluster["deleted_size"],
                "total_supporting_reads": supporting_reads,
                "samples_with_signal": len(cluster["supporting_samples"]),
                "supporting_samples": ",".join(sorted(cluster["supporting_samples"])),
                "rotation_support": ",".join(sorted(cluster["rotations"])),
                "rotation_count": len(cluster["rotations"]),
                "reported_orientation_support": ",".join(sorted(cluster["reported_orientations"])),
            }
        )
    read_fields = [
        "junction_id",
        "sample",
        "species",
        "read_id",
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "reported_left_breakpoint",
        "reported_right_breakpoint",
        "reported_deleted_size",
        "canonical_orientation",
        "left_anchor_length",
        "right_anchor_length",
        "strand",
        "source",
        "rotation_name",
        "rotation_start",
        "raw_left_breakpoint",
        "raw_right_breakpoint",
    ]
    write_tsv(args.all_reads, reads, fieldnames=read_fields)
    write_tsv(
        args.clusters,
        cluster_rows,
        fieldnames=[
            "junction_id",
            "species",
            "left_breakpoint",
            "right_breakpoint",
            "deleted_size",
            "total_supporting_reads",
            "samples_with_signal",
            "supporting_samples",
            "rotation_support",
            "rotation_count",
            "reported_orientation_support",
        ],
    )
    write_tsv(
        args.id_map,
        id_rows,
        fieldnames=[
            "sample",
            "read_id",
            "junction_id",
            "left_breakpoint",
            "right_breakpoint",
            "reported_left_breakpoint",
            "reported_right_breakpoint",
        ],
    )


if __name__ == "__main__":
    main()
