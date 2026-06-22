#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict

import pysam

from circular_deletions import canonical_breakpoints, deletion_id, normalize_pos
from common import ensure_parent


def read_key(read: pysam.AlignedSegment) -> str:
    return read.query_name.removesuffix("/1").removesuffix("/2")


def aligned_query_fraction(read: pysam.AlignedSegment) -> float:
    length = read.infer_query_length(always=True) or len(read.query_sequence or "") or 0
    if length <= 0:
        return 0.0
    return max(0, int(read.query_alignment_end or 0) - int(read.query_alignment_start or 0)) / length


def soft_clip_fraction(read: pysam.AlignedSegment) -> float:
    length = read.infer_query_length(always=True) or len(read.query_sequence or "") or 0
    if length <= 0 or not read.cigartuples:
        return 0.0
    soft = sum(size for op, size in read.cigartuples if op == 4)
    return soft / length


def segment_from_read(read: pysam.AlignedSegment, args: argparse.Namespace) -> dict | None:
    if read.is_unmapped:
        return None
    if read.mapping_quality < args.min_mapq:
        return None
    if read.is_secondary and not args.include_secondary:
        return None
    if read.is_supplementary and not args.include_supplementary:
        return None
    q_start = int(read.query_alignment_start or 0)
    q_end = int(read.query_alignment_end or 0)
    anchor = max(0, q_end - q_start)
    if anchor < args.min_anchor_length:
        return None
    if aligned_query_fraction(read) < args.min_segment_aligned_fraction:
        return None
    if soft_clip_fraction(read) > args.max_soft_clip_fraction:
        return None
    ref_start = int(read.reference_start) + 1
    ref_end = int(read.reference_end)
    return {
        "read_id": read_key(read),
        "query_name": read.query_name,
        "query_start": q_start,
        "query_end": q_end,
        "anchor_length": anchor,
        "ref_start_raw": ref_start,
        "ref_end_raw": ref_end,
        "ref_start": normalize_pos(ref_start, args.mt_length, args.rotation_start),
        "ref_end": normalize_pos(ref_end, args.mt_length, args.rotation_start),
        "strand": "-" if read.is_reverse else "+",
        "mapq": int(read.mapping_quality),
        "is_secondary": "yes" if read.is_secondary else "no",
        "is_supplementary": "yes" if read.is_supplementary else "no",
        "aligned_fraction": aligned_query_fraction(read),
        "soft_clip_fraction": soft_clip_fraction(read),
    }


def query_compatible(a: dict, b: dict, max_overlap: int, max_gap: int) -> bool:
    overlap = max(0, min(a["query_end"], b["query_end"]) - max(a["query_start"], b["query_start"]))
    gap = max(0, b["query_start"] - a["query_end"])
    return overlap <= max_overlap and gap <= max_gap


def deletion_from_segments(a: dict, b: dict, args: argparse.Namespace) -> dict | None:
    if a["strand"] != b["strand"]:
        return None
    if a["query_start"] > b["query_start"]:
        a, b = b, a
    if not query_compatible(a, b, args.max_query_overlap_bp, args.max_query_gap_bp):
        return None
    if a["strand"] == "+":
        reported_left = a["ref_end"]
        reported_right = b["ref_start"]
        raw_left = a["ref_end_raw"]
        raw_right = b["ref_start_raw"]
    else:
        reported_left = b["ref_end"]
        reported_right = a["ref_start"]
        raw_left = b["ref_end_raw"]
        raw_right = a["ref_start_raw"]
    canonical = canonical_breakpoints(reported_left, reported_right, args.mt_length)
    size = int(canonical["deleted_size"])
    if size < args.min_deletion_size or size > args.max_deletion_size:
        return None
    left = int(canonical["canonical_left_breakpoint"])
    right = int(canonical["canonical_right_breakpoint"])
    return {
        "sample": args.sample,
        "species": args.species,
        "read_id": a["read_id"],
        "deletion_id": deletion_id(left, right, size),
        "left_breakpoint": left,
        "right_breakpoint": right,
        "deleted_size": size,
        "wraps_origin": canonical["wraps_origin"],
        "deleted_interval": canonical["deleted_interval"],
        "reported_left_breakpoint": reported_left,
        "reported_right_breakpoint": reported_right,
        "reported_deleted_size": canonical.get("reported_deleted_size", ""),
        "canonical_orientation": canonical["canonical_orientation"],
        "left_anchor_length": a["anchor_length"],
        "right_anchor_length": b["anchor_length"],
        "min_anchor_length": min(a["anchor_length"], b["anchor_length"]),
        "strand": a["strand"],
        "source": "minimap2_split_alignment",
        "rotation_name": args.rotation_name,
        "rotation_start": args.rotation_start,
        "raw_left_breakpoint": raw_left,
        "raw_right_breakpoint": raw_right,
        "left_mapq": a["mapq"],
        "right_mapq": b["mapq"],
        "min_mapq": min(a["mapq"], b["mapq"]),
        "query_overlap_bp": max(0, min(a["query_end"], b["query_end"]) - max(a["query_start"], b["query_start"])),
        "query_gap_bp": max(0, b["query_start"] - a["query_end"]),
    }


def write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


FIELDS = [
    "sample",
    "species",
    "read_id",
    "deletion_id",
    "left_breakpoint",
    "right_breakpoint",
    "deleted_size",
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--rotation-start", type=int, required=True)
    parser.add_argument("--rotation-name", required=True)
    parser.add_argument("--min-anchor-length", type=int, required=True)
    parser.add_argument("--min-deletion-size", type=int, required=True)
    parser.add_argument("--max-deletion-size", type=int, required=True)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--min-segment-aligned-fraction", type=float, default=0.15)
    parser.add_argument("--max-soft-clip-fraction", type=float, default=0.9)
    parser.add_argument("--max-query-overlap-bp", type=int, default=10)
    parser.add_argument("--max-query-gap-bp", type=int, default=20)
    parser.add_argument("--include-secondary", action="store_true")
    parser.add_argument("--include-supplementary", action="store_true")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--filtered", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    by_read: dict[str, list[dict]] = defaultdict(list)
    counts = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            counts["alignment_records"] += 1
            segment = segment_from_read(read, args)
            if segment is None:
                counts["segments_rejected"] += 1
                continue
            counts["segments_used"] += 1
            by_read[segment["read_id"]].append(segment)

    candidate_rows = []
    seen = set()
    for read_id, segments in by_read.items():
        if len(segments) < 2:
            continue
        segments = sorted(segments, key=lambda item: (item["query_start"], item["query_end"], item["ref_start_raw"]))
        for i, a in enumerate(segments):
            for b in segments[i + 1 :]:
                row = deletion_from_segments(a, b, args)
                if row is None:
                    continue
                key = (read_id, row["deletion_id"], row["rotation_name"])
                if key in seen:
                    continue
                seen.add(key)
                candidate_rows.append(row)

    write_rows(args.candidates, candidate_rows, FIELDS)
    write_rows(args.filtered, candidate_rows, FIELDS)
    write_rows(
        args.summary,
        [
            {"metric": "alignment_records", "value": counts["alignment_records"]},
            {"metric": "segments_used", "value": counts["segments_used"]},
            {"metric": "segments_rejected", "value": counts["segments_rejected"]},
            {"metric": "reads_with_usable_segments", "value": len(by_read)},
            {"metric": "deletion_supporting_records", "value": len(candidate_rows)},
        ],
        ["metric", "value"],
    )


if __name__ == "__main__":
    main()
