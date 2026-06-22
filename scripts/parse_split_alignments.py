#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from collections import Counter

from common import write_tsv


FIELDS = [
    "sample",
    "species",
    "read_id",
    "left_breakpoint",
    "right_breakpoint",
    "deleted_size",
    "left_anchor_length",
    "right_anchor_length",
    "left_mapq",
    "right_mapq",
    "strand",
    "source",
    "rotation_name",
    "rotation_start",
    "raw_left_breakpoint",
    "raw_right_breakpoint",
    "filter_status",
    "filter_reason",
]

CIGAR_ALIGNED_RE = re.compile(r"(\d+)([M=X])")


def normalize_pos(pos: int, mt_length: int, padding: int = 0, rotation_start: int = 1) -> int:
    value = ((rotation_start - 1 + pos - padding - 1) % mt_length) + 1
    return value


def deletion_size(left: int, right: int, mt_length: int) -> int:
    if right > left:
        return right - left - 1
    return mt_length - left + right - 1


def aligned_bases_from_cigar(cigar: str) -> int:
    return sum(int(length) for length, op in CIGAR_ALIGNED_RE.findall(cigar))


def parse_star_junction_line(
    line: str,
    sample: str,
    species: str,
    mt_length: int,
    padding: int = 0,
    rotation_start: int = 1,
    rotation_name: str = "normal",
) -> dict | None:
    if not line.strip() or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 14:
        parts = line.split()
    if len(parts) < 14:
        return None
    try:
        raw_left_pos = int(parts[1])
        raw_right_pos = int(parts[4])
    except ValueError:
        return None
    left_pos = normalize_pos(raw_left_pos, mt_length, padding, rotation_start)
    right_pos = normalize_pos(raw_right_pos, mt_length, padding, rotation_start)
    read_id = parts[9]
    left_anchor = aligned_bases_from_cigar(parts[11])
    right_anchor = aligned_bases_from_cigar(parts[13])
    strand = f"{parts[2]}:{parts[5]}"
    return {
        "sample": sample,
        "species": species,
        "read_id": read_id,
        "left_breakpoint": left_pos,
        "right_breakpoint": right_pos,
        "deleted_size": deletion_size(left_pos, right_pos, mt_length),
        "left_anchor_length": left_anchor,
        "right_anchor_length": right_anchor,
        "left_mapq": "",
        "right_mapq": "",
        "strand": strand,
        "source": "STAR_Chimeric.out.junction",
        "rotation_name": rotation_name,
        "rotation_start": rotation_start,
        "raw_left_breakpoint": raw_left_pos,
        "raw_right_breakpoint": raw_right_pos,
    }


def filter_row(row: dict, min_anchor: int, min_deletion_size: int, max_deletion_size: int) -> dict:
    reasons = []
    if int(row["left_anchor_length"]) < min_anchor:
        reasons.append("left_anchor_short")
    if int(row["right_anchor_length"]) < min_anchor:
        reasons.append("right_anchor_short")
    if int(row["deleted_size"]) < min_deletion_size:
        reasons.append("deletion_too_small")
    if int(row["deleted_size"]) > max_deletion_size:
        reasons.append("deletion_too_large")
    row = dict(row)
    row["filter_status"] = "pass" if not reasons else "fail"
    row["filter_reason"] = ";".join(reasons)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--junction-file", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--padding", type=int, default=0)
    parser.add_argument("--rotation-start", type=int, default=1)
    parser.add_argument("--rotation-name", default="normal")
    parser.add_argument("--min-anchor-length", type=int, required=True)
    parser.add_argument("--min-deletion-size", type=int, required=True)
    parser.add_argument("--max-deletion-size", type=int, required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--filtered", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    candidates = []
    with open(args.junction_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            row = parse_star_junction_line(
                line,
                args.sample,
                args.species,
                args.mt_length,
                args.padding,
                args.rotation_start,
                args.rotation_name,
            )
            if row is not None:
                candidates.append(row)
    judged = [
        filter_row(row, args.min_anchor_length, args.min_deletion_size, args.max_deletion_size)
        for row in candidates
    ]
    passed = [row for row in judged if row["filter_status"] == "pass"]
    write_tsv(args.candidates, judged, fieldnames=FIELDS)
    write_tsv(args.filtered, passed, fieldnames=FIELDS)
    counts = Counter(row["filter_status"] for row in judged)
    write_tsv(
        args.summary,
        [{"sample": args.sample, "candidate_reads": len(judged), "passing_reads": counts["pass"], "failed_reads": counts["fail"]}],
        fieldnames=["sample", "candidate_reads", "passing_reads", "failed_reads"],
    )


if __name__ == "__main__":
    main()
