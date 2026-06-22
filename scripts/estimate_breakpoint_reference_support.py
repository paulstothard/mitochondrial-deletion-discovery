#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam

from circular_deletions import normalize_pos
from common import ensure_parent


def read_key(name: object) -> str:
    return str(name or "").strip().removesuffix("/1").removesuffix("/2")


def circular_window(pos: int, radius: int, mt_length: int) -> list[tuple[int, int]]:
    start = ((int(pos) - int(radius) - 1) % mt_length) + 1
    end = ((int(pos) + int(radius) - 1) % mt_length) + 1
    if start <= end:
        return [(start, end)]
    return [(start, mt_length), (1, end)]


def circular_block_pieces(start: int, end: int, mt_length: int, rotation_start: int) -> list[tuple[int, int]]:
    norm_start = normalize_pos(start, mt_length, rotation_start)
    norm_end = normalize_pos(end, mt_length, rotation_start)
    if norm_start <= norm_end:
        return [(norm_start, norm_end)]
    return [(norm_start, mt_length), (1, norm_end)]


def piece_covered(target: tuple[int, int], blocks: list[tuple[int, int]]) -> bool:
    target_start, target_end = target
    return any(block_start <= target_start and block_end >= target_end for block_start, block_end in blocks)


def window_covered(pos: int, radius: int, mt_length: int, blocks: list[tuple[int, int]]) -> bool:
    return all(piece_covered(piece, blocks) for piece in circular_window(pos, radius, mt_length))


def parse_bam_path(path: str) -> tuple[str, str]:
    p = Path(path)
    sample = p.name.removesuffix(".bam")
    rotation = p.parent.name
    return sample, rotation


def rotation_start(rotation: str, starts: dict[str, int]) -> int:
    return int(starts.get(rotation, 1))


def candidate_positions(clusters: pd.DataFrame) -> dict[str, set[int]]:
    by_sample: dict[str, set[int]] = defaultdict(set)
    if clusters.empty:
        return by_sample
    for _, row in clusters.iterrows():
        samples = [value for value in str(row.get("supporting_samples", "")).split(",") if value]
        if not samples and "sample" in row:
            samples = [str(row["sample"])]
        for sample in samples:
            by_sample[sample].add(int(row["left_breakpoint"]))
            by_sample[sample].add(int(row["right_breakpoint"]))
    return by_sample


def add_depth_range(diff: list[int], start: int, end: int) -> None:
    if start > end:
        return
    diff[start] += 1
    diff[end + 1] -= 1


def bam_spanning_depths(
    bam_path: str,
    rotation_starts: dict[str, int],
    mt_length: int,
    window_bp: int,
    min_mapq: int,
) -> tuple[str, dict[int, int]]:
    sample, rotation = parse_bam_path(bam_path)
    rot_start = rotation_start(rotation, rotation_starts)
    diff = [0] * (mt_length + 3)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if int(read.mapping_quality) < min_mapq:
                continue
            raw_blocks = read.get_blocks()
            if not raw_blocks:
                continue
            for raw_start0, raw_end0 in raw_blocks:
                raw_start = int(raw_start0) + 1
                raw_end = int(raw_end0)
                for block_start, block_end in circular_block_pieces(raw_start, raw_end, mt_length, rot_start):
                    add_depth_range(diff, block_start + window_bp, block_end - window_bp)
    depths = {}
    current = 0
    for pos in range(1, mt_length + 1):
        current += diff[pos]
        depths[pos] = current
    return sample, depths


def count_reference_spanners(
    clusters: pd.DataFrame,
    bam_paths: list[str],
    rotation_starts: dict[str, int],
    mt_length: int,
    window_bp: int,
    min_mapq: int,
) -> pd.DataFrame:
    positions = candidate_positions(clusters)
    reference_counts: dict[tuple[str, int], int] = defaultdict(int)

    for bam_path in bam_paths:
        sample, rotation = parse_bam_path(bam_path)
        if sample not in positions or not Path(bam_path).exists():
            continue
        depth_sample, depths = bam_spanning_depths(bam_path, rotation_starts, mt_length, window_bp, min_mapq)
        for breakpoint in positions[depth_sample]:
            key = (depth_sample, breakpoint)
            reference_counts[key] = max(reference_counts[key], int(depths.get(breakpoint, 0)))

    rows = []
    for _, row in clusters.iterrows():
        exact_id = str(row["exact_deletion_id"])
        samples = [value for value in str(row.get("supporting_samples", "")).split(",") if value]
        left_bp = int(row["left_breakpoint"])
        right_bp = int(row["right_breakpoint"])
        left_count = sum(reference_counts.get((sample, left_bp), 0) for sample in samples)
        right_count = sum(reference_counts.get((sample, right_bp), 0) for sample in samples)
        ref_min = min(left_count, right_count)
        split_support = int(row.get("total_supporting_reads", 0) or 0)
        denominator = split_support + ref_min
        fraction = split_support / denominator if denominator else float("nan")
        rows.append(
            {
                "exact_deletion_id": exact_id,
                "split_supporting_reads": split_support,
                "left_reference_spanning_reads": left_count,
                "right_reference_spanning_reads": right_count,
                "reference_spanning_reads_min": ref_min,
                "local_split_support_fraction": fraction,
                "reference_support_window_bp": window_bp,
                "reference_support_method": "primary_alignment_depth_max_across_rotations",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--all-reads")
    parser.add_argument("--bam", nargs="+", required=True)
    parser.add_argument("--rotation-starts", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--window-bp", type=int, default=20)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-reference-support", required=True)
    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters, sep="\t")
    rotation_starts = {}
    for item in args.rotation_starts.split(","):
        if not item:
            continue
        name, value = item.split(":", 1)
        rotation_starts[name] = int(value)

    support = count_reference_spanners(
        clusters,
        args.bam,
        rotation_starts,
        args.mt_length,
        args.window_bp,
        args.min_mapq,
    )
    if not clusters.empty and not support.empty:
        merged = clusters.merge(support, on="exact_deletion_id", how="left")
    else:
        merged = clusters.copy()
    ensure_parent(args.out_clusters)
    ensure_parent(args.out_reference_support)
    merged.to_csv(args.out_clusters, sep="\t", index=False)
    support.to_csv(args.out_reference_support, sep="\t", index=False)


if __name__ == "__main__":
    main()
