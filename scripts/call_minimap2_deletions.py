#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from itertools import combinations

import pysam

from circular_deletions import (
    breakpoint_pair_id,
    canonical_breakpoints,
    deletion_id,
    directed_breakpoints,
    normalize_pos,
)
from common import ensure_parent


def read_key(read: pysam.AlignedSegment) -> str:
    return read.query_name.removesuffix("/1").removesuffix("/2")


def alignment_chain_key(read: pysam.AlignedSegment) -> str:
    """Identify one physical read sequence without merging paired-end mates."""
    fragment = read_key(read)
    if read.is_read1:
        return f"{fragment}/1"
    if read.is_read2:
        return f"{fragment}/2"
    if read.query_name.endswith(("/1", "/2")):
        return read.query_name
    return fragment


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


def alignment_tag(read: pysam.AlignedSegment, tag: str) -> object:
    try:
        return read.get_tag(tag)
    except KeyError:
        return ""


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
        "alignment_chain_id": alignment_chain_key(read),
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
        "cigar": read.cigarstring or "",
        "flag": int(read.flag),
        "is_primary": "no" if read.is_secondary or read.is_supplementary else "yes",
        "alignment_score": alignment_tag(read, "AS"),
        "edit_distance": alignment_tag(read, "NM"),
        "sa_tag": alignment_tag(read, "SA"),
        "minimap2_type": alignment_tag(read, "tp"),
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
    query_first = a
    query_second = b

    # SAM/BAM stores reverse-strand SEQ reverse-complemented, so CIGAR query
    # coordinates advance with the reference for same-strand split records.
    # Reversing the segment order again would select the reciprocal arc.
    left_segment = a
    right_segment = b
    directed_left = a["ref_end"]
    directed_right = b["ref_start"]
    raw_left = a["ref_end_raw"]
    raw_right = b["ref_start_raw"]

    directed = directed_breakpoints(directed_left, directed_right, args.mt_length)
    if args.arc_assignment == "legacy_shortest_arc":
        legacy = canonical_breakpoints(directed_left, directed_right, args.mt_length)
        left = int(legacy["canonical_left_breakpoint"])
        right = int(legacy["canonical_right_breakpoint"])
        arc = directed_breakpoints(left, right, args.mt_length)
        arc["arc_assignment_method"] = "legacy_shortest_arc"
        canonical_orientation = legacy["canonical_orientation"]
    else:
        left = int(directed["left_breakpoint"])
        right = int(directed["right_breakpoint"])
        arc = directed
        canonical_orientation = "alignment_directed"
    size = int(arc["deleted_size"])
    if size < args.min_deletion_size or size > args.max_deletion_size:
        return None
    return {
        "sample": args.sample,
        "species": args.species,
        "read_id": a["read_id"],
        "deletion_id": deletion_id(left, right, size),
        "breakpoint_pair_id": breakpoint_pair_id(left, right),
        "left_breakpoint": left,
        "right_breakpoint": right,
        "deleted_size": size,
        "wraps_origin": arc["wraps_origin"],
        "deleted_interval": arc["deleted_interval"],
        "complement_deleted_size": arc["complement_deleted_size"],
        "complement_wraps_origin": arc["complement_wraps_origin"],
        "arc_assignment_method": arc["arc_assignment_method"],
        "direction_status": "directed",
        "directed_left_breakpoint": directed_left,
        "directed_right_breakpoint": directed_right,
        "directed_deleted_size": directed["deleted_size"],
        "reported_left_breakpoint": directed_left,
        "reported_right_breakpoint": directed_right,
        "reported_deleted_size": directed["deleted_size"],
        "canonical_orientation": canonical_orientation,
        "left_anchor_length": left_segment["anchor_length"],
        "right_anchor_length": right_segment["anchor_length"],
        "min_anchor_length": min(left_segment["anchor_length"], right_segment["anchor_length"]),
        "strand": a["strand"],
        "source": "minimap2_split_alignment",
        "rotation_name": args.rotation_name,
        "rotation_start": args.rotation_start,
        "alignment_chain_id": f"{a.get('alignment_chain_id', a['read_id'])}:{args.rotation_name}",
        "raw_left_breakpoint": raw_left,
        "raw_right_breakpoint": raw_right,
        "left_mapq": left_segment["mapq"],
        "right_mapq": right_segment["mapq"],
        "min_mapq": min(left_segment["mapq"], right_segment["mapq"]),
        "query_overlap_bp": max(0, min(a["query_end"], b["query_end"]) - max(a["query_start"], b["query_start"])),
        "query_gap_bp": max(0, b["query_start"] - a["query_end"]),
        "query_first_start": query_first["query_start"],
        "query_first_end": query_first["query_end"],
        "query_second_start": query_second["query_start"],
        "query_second_end": query_second["query_end"],
        "left_segment_ref_start": left_segment["ref_start"],
        "left_segment_ref_end": left_segment["ref_end"],
        "right_segment_ref_start": right_segment["ref_start"],
        "right_segment_ref_end": right_segment["ref_end"],
        "left_segment_ref_start_raw": left_segment["ref_start_raw"],
        "left_segment_ref_end_raw": left_segment["ref_end_raw"],
        "right_segment_ref_start_raw": right_segment["ref_start_raw"],
        "right_segment_ref_end_raw": right_segment["ref_end_raw"],
        "left_cigar": left_segment["cigar"],
        "right_cigar": right_segment["cigar"],
        "left_flag": left_segment["flag"],
        "right_flag": right_segment["flag"],
        "left_is_primary": left_segment["is_primary"],
        "right_is_primary": right_segment["is_primary"],
        "left_is_secondary": left_segment["is_secondary"],
        "right_is_secondary": right_segment["is_secondary"],
        "left_is_supplementary": left_segment["is_supplementary"],
        "right_is_supplementary": right_segment["is_supplementary"],
        "left_alignment_score": left_segment["alignment_score"],
        "right_alignment_score": right_segment["alignment_score"],
        "left_edit_distance": left_segment["edit_distance"],
        "right_edit_distance": right_segment["edit_distance"],
        "left_sa_tag": left_segment["sa_tag"],
        "right_sa_tag": right_segment["sa_tag"],
        "left_minimap2_type": left_segment["minimap2_type"],
        "right_minimap2_type": right_segment["minimap2_type"],
    }


def write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def candidate_rows_from_chains(by_chain: dict[str, list[dict]], args: argparse.Namespace) -> list[dict]:
    candidate_rows = []
    seen = set()
    for segments in by_chain.values():
        if len(segments) < 2:
            continue
        segments = sorted(segments, key=lambda item: (item["query_start"], item["query_end"], item["ref_start_raw"]))
        pairs = zip(segments, segments[1:]) if args.pairing_mode == "adjacent" else combinations(segments, 2)
        for a, b in pairs:
            row = deletion_from_segments(a, b, args)
            if row is None:
                continue
            key = (row["read_id"], row["deletion_id"], row["rotation_name"])
            if key in seen:
                continue
            seen.add(key)
            candidate_rows.append(row)
    return candidate_rows


FIELDS = [
    "sample",
    "species",
    "read_id",
    "deletion_id",
    "breakpoint_pair_id",
    "left_breakpoint",
    "right_breakpoint",
    "deleted_size",
    "wraps_origin",
    "deleted_interval",
    "complement_deleted_size",
    "complement_wraps_origin",
    "arc_assignment_method",
    "direction_status",
    "directed_left_breakpoint",
    "directed_right_breakpoint",
    "directed_deleted_size",
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
    "alignment_chain_id",
    "raw_left_breakpoint",
    "raw_right_breakpoint",
    "left_mapq",
    "right_mapq",
    "min_mapq",
    "query_overlap_bp",
    "query_gap_bp",
    "query_first_start",
    "query_first_end",
    "query_second_start",
    "query_second_end",
    "left_segment_ref_start",
    "left_segment_ref_end",
    "right_segment_ref_start",
    "right_segment_ref_end",
    "left_segment_ref_start_raw",
    "left_segment_ref_end_raw",
    "right_segment_ref_start_raw",
    "right_segment_ref_end_raw",
    "left_cigar",
    "right_cigar",
    "left_flag",
    "right_flag",
    "left_is_primary",
    "right_is_primary",
    "left_is_secondary",
    "right_is_secondary",
    "left_is_supplementary",
    "right_is_supplementary",
    "left_alignment_score",
    "right_alignment_score",
    "left_edit_distance",
    "right_edit_distance",
    "left_sa_tag",
    "right_sa_tag",
    "left_minimap2_type",
    "right_minimap2_type",
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
    parser.add_argument("--arc-assignment", choices=["alignment_directed", "legacy_shortest_arc"], default="alignment_directed")
    parser.add_argument("--pairing-mode", choices=["adjacent", "all_compatible"], default="all_compatible")
    parser.add_argument("--ambiguous-direction-policy", choices=["exclude", "include"], default="exclude")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--filtered", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    by_chain: dict[str, list[dict]] = defaultdict(list)
    fragments_with_segments = set()
    counts = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            counts["alignment_records"] += 1
            segment = segment_from_read(read, args)
            if segment is None:
                counts["segments_rejected"] += 1
                continue
            counts["segments_used"] += 1
            by_chain[segment["alignment_chain_id"]].append(segment)
            fragments_with_segments.add(segment["read_id"])

    candidate_rows = candidate_rows_from_chains(by_chain, args)

    by_read_pair: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in candidate_rows:
        by_read_pair[(row["read_id"], row["breakpoint_pair_id"], row["rotation_name"])].append(row)
    ambiguous_keys = {
        key
        for key, rows in by_read_pair.items()
        if len({(int(row["directed_left_breakpoint"]), int(row["directed_right_breakpoint"])) for row in rows}) > 1
    }
    for row in candidate_rows:
        key = (row["read_id"], row["breakpoint_pair_id"], row["rotation_name"])
        if key in ambiguous_keys:
            row["direction_status"] = "ambiguous_reciprocal_alignments"
    filtered_rows = [
        row
        for row in candidate_rows
        if args.ambiguous_direction_policy == "include" or row["direction_status"] == "directed"
    ]

    write_rows(args.candidates, candidate_rows, FIELDS)
    write_rows(args.filtered, filtered_rows, FIELDS)
    write_rows(
        args.summary,
        [
            {"metric": "alignment_records", "value": counts["alignment_records"]},
            {"metric": "segments_used", "value": counts["segments_used"]},
            {"metric": "segments_rejected", "value": counts["segments_rejected"]},
            {"metric": "reads_with_usable_segments", "value": len(by_chain)},
            {"metric": "fragments_with_usable_segments", "value": len(fragments_with_segments)},
            {"metric": "candidate_deletion_records", "value": len(candidate_rows)},
            {"metric": "ambiguous_direction_records", "value": len(candidate_rows) - len(filtered_rows)},
            {"metric": "deletion_supporting_records", "value": len(filtered_rows)},
        ],
        ["metric", "value"],
    )


if __name__ == "__main__":
    main()
