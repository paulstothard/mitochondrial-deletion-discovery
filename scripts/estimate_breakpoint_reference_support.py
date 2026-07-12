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


def eligible_cluster_samples(clusters: pd.DataFrame, observations: pd.DataFrame) -> dict[str, set[str]]:
    eligible: dict[str, set[str]] = defaultdict(set)
    if not observations.empty and "exact_deletion_id" in observations.columns and "minimap2_support" in observations.columns:
        minimap = observations.loc[
            observations["minimap2_support"].fillna("").astype(str).str.lower().eq("yes")
        ]
        for _, row in minimap.iterrows():
            eligible[str(row["exact_deletion_id"])].add(str(row["sample"]))
        return eligible
    for _, row in clusters.iterrows():
        exact_id = str(row["exact_deletion_id"])
        samples = [value for value in str(row.get("supporting_samples", "")).split(",") if value]
        if not samples and "sample" in row:
            samples = [str(row["sample"])]
        eligible[exact_id].update(samples)
    return eligible


def candidate_positions(clusters: pd.DataFrame, eligible: dict[str, set[str]]) -> dict[str, set[int]]:
    by_sample: dict[str, set[int]] = defaultdict(set)
    if clusters.empty:
        return by_sample
    for _, row in clusters.iterrows():
        for sample in eligible.get(str(row["exact_deletion_id"]), set()):
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


def original_to_rotated_pos(pos: int, mt_length: int, rotation_start: int) -> int:
    return ((int(pos) - int(rotation_start)) % int(mt_length)) + 1


def fetch_pieces(pos: int, radius: int, mt_length: int, rotation_start: int) -> list[tuple[int, int]]:
    raw_pos = original_to_rotated_pos(pos, mt_length, rotation_start)
    return circular_window(raw_pos, radius, mt_length)


def bam_indexed_spanning_depths(
    bam_path: str,
    positions: set[int],
    rotation_starts: dict[str, int],
    mt_length: int,
    window_bp: int,
    min_mapq: int,
) -> tuple[str, dict[int, int]]:
    sample, rotation = parse_bam_path(bam_path)
    rot_start = rotation_start(rotation, rotation_starts)
    depths = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index() or not bam.references:
            depth_sample, all_depths = bam_spanning_depths(
                bam_path, rotation_starts, mt_length, window_bp, min_mapq
            )
            return depth_sample, {position: all_depths.get(position, 0) for position in positions}
        contig = bam.references[0]
        for position in positions:
            candidates = {}
            for start, end in fetch_pieces(position, window_bp, mt_length, rot_start):
                for read in bam.fetch(contig, max(0, start - 1), min(mt_length, end)):
                    candidates[(read.query_name, read.flag, read.reference_start, read.cigarstring)] = read
            spanning = set()
            for read in candidates.values():
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if int(read.mapping_quality) < min_mapq:
                    continue
                blocks = []
                for raw_start0, raw_end0 in read.get_blocks():
                    blocks.extend(
                        circular_block_pieces(int(raw_start0) + 1, int(raw_end0), mt_length, rot_start)
                    )
                if blocks and window_covered(position, window_bp, mt_length, blocks):
                    spanning.add(read.query_name)
            depths[position] = len(spanning)
    return sample, depths


def count_reference_spanners(
    clusters: pd.DataFrame,
    observations: pd.DataFrame,
    bam_paths: list[str],
    rotation_starts: dict[str, int],
    mt_length: int,
    window_bp: int,
    min_mapq: int,
) -> pd.DataFrame:
    eligible = eligible_cluster_samples(clusters, observations)
    positions = candidate_positions(clusters, eligible)
    reference_counts: dict[tuple[str, int], int] = defaultdict(int)

    for bam_path in bam_paths:
        sample, rotation = parse_bam_path(bam_path)
        if sample not in positions or not Path(bam_path).exists():
            continue
        depth_sample, depths = bam_indexed_spanning_depths(
            bam_path,
            positions[sample],
            rotation_starts,
            mt_length,
            window_bp,
            min_mapq,
        )
        for breakpoint in positions[depth_sample]:
            key = (depth_sample, breakpoint)
            reference_counts[key] = max(reference_counts[key], int(depths.get(breakpoint, 0)))

    rows = []
    for _, row in clusters.iterrows():
        exact_id = str(row["exact_deletion_id"])
        samples = sorted(eligible.get(exact_id, set()))
        left_bp = int(row["left_breakpoint"])
        right_bp = int(row["right_breakpoint"])
        left_count = sum(reference_counts.get((sample, left_bp), 0) for sample in samples)
        right_count = sum(reference_counts.get((sample, right_bp), 0) for sample in samples)
        ref_min = min(left_count, right_count)
        if not observations.empty and "minimap2_support" in observations.columns:
            split_support = int(
                (
                    observations["exact_deletion_id"].astype(str).eq(exact_id)
                    & observations["minimap2_support"].fillna("").astype(str).str.lower().eq("yes")
                ).sum()
            )
        else:
            split_support = int(row.get("total_supporting_reads", 0) or 0)
        denominator = split_support + ref_min
        fraction = split_support / denominator if denominator and samples else float("nan")
        rows.append(
            {
                "exact_deletion_id": exact_id,
                "split_supporting_reads": split_support,
                "left_reference_spanning_reads": left_count if samples else float("nan"),
                "right_reference_spanning_reads": right_count if samples else float("nan"),
                "reference_spanning_reads_min": ref_min if samples else float("nan"),
                "local_split_support_fraction": fraction,
                "reference_support_window_bp": window_bp,
                "reference_support_method": (
                    "indexed_primary_alignment_depth_max_across_rotations"
                    if samples
                    else "not_available_without_minimap2_remap_evidence"
                ),
            }
        )
    return pd.DataFrame(rows)


def reuse_reference_support(
    clusters: pd.DataFrame,
    observations: pd.DataFrame,
    existing: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    eligible = eligible_cluster_samples(clusters, observations)
    cached = {
        str(row["exact_deletion_id"]): row.to_dict()
        for _, row in existing.iterrows()
        if str(row.get("exact_deletion_id", ""))
    }
    rows = []
    for _, cluster in clusters.iterrows():
        exact_id = str(cluster["exact_deletion_id"])
        if exact_id in cached:
            row = dict(cached[exact_id])
            row["reference_support_method"] = "reused_existing_" + str(
                row.get("reference_support_method", "remap_reference_support")
            )
            rows.append(row)
            continue
        has_minimap = bool(eligible.get(exact_id))
        rows.append(
            {
                "exact_deletion_id": exact_id,
                "split_supporting_reads": float("nan"),
                "left_reference_spanning_reads": float("nan"),
                "right_reference_spanning_reads": float("nan"),
                "reference_spanning_reads_min": float("nan"),
                "local_split_support_fraction": float("nan"),
                "reference_support_window_bp": window_bp,
                "reference_support_method": (
                    "not_available_in_existing_remap_support"
                    if has_minimap
                    else "not_available_without_minimap2_remap_evidence"
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--all-reads")
    parser.add_argument("--existing-reference-support", default="")
    parser.add_argument("--bam", nargs="+", required=True)
    parser.add_argument("--rotation-starts", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--window-bp", type=int, default=20)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-reference-support", required=True)
    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters, sep="\t")
    observations = pd.read_csv(args.all_reads, sep="\t") if args.all_reads else pd.DataFrame()
    rotation_starts = {}
    for item in args.rotation_starts.split(","):
        if not item:
            continue
        name, value = item.split(":", 1)
        rotation_starts[name] = int(value)

    existing = (
        pd.read_csv(args.existing_reference_support, sep="\t")
        if args.existing_reference_support and Path(args.existing_reference_support).exists()
        else pd.DataFrame()
    )
    existing_windows = set(
        pd.to_numeric(existing.get("reference_support_window_bp", pd.Series(dtype=float)), errors="coerce")
        .dropna()
        .astype(int)
    )
    if not existing.empty and (not existing_windows or existing_windows == {args.window_bp}):
        support = reuse_reference_support(clusters, observations, existing, args.window_bp)
    else:
        support = count_reference_spanners(
            clusters,
            observations,
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
