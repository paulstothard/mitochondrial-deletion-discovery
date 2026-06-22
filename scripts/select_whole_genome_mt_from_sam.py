#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter

import pysam

from classify_mt_reads import fastq_entry, read_key
from common import empty_gzip, ensure_parent, write_json


def aligned_fraction(read: pysam.AlignedSegment) -> float:
    query_length = read.query_length or len(read.query_sequence or "")
    if not query_length:
        return 0.0
    return float(read.query_alignment_length or 0) / float(query_length)


def is_mapped_to_mt(read: pysam.AlignedSegment, mt_names: set[str]) -> bool:
    return (not read.is_unmapped) and read.reference_name in mt_names


def is_strong_competitor(read: pysam.AlignedSegment, mt_names: set[str], min_aligned_fraction: float) -> bool:
    return (not read.is_unmapped) and read.reference_name not in mt_names and aligned_fraction(read) >= min_aligned_fraction


def write_classification(path: str, counts: Counter) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["table_type", "name", "count"])
        for name, count in sorted(counts.items()):
            writer.writerow(["alignment_category", name, count])
        writer.writerow(["selection_reason", "mt_primary_best", counts["mt_primary_best_written"]])
        writer.writerow(["selection_reason", "mt_primary_ambiguous", counts["mt_primary_ambiguous_written"]])


def classify_group(
    records: list[pysam.AlignedSegment],
    mt_names: set[str],
    min_mt_mapq: int,
    min_mt_aligned_fraction: float,
    ambiguous_mapq_below: int,
    competing_nuclear_aligned_fraction: float,
) -> tuple[str, pysam.AlignedSegment | None]:
    primary = next((read for read in records if not read.is_secondary and not read.is_supplementary), None)
    if primary is None:
        return "no_primary_alignment", None
    if primary.is_unmapped:
        return "primary_unmapped", primary
    if not is_mapped_to_mt(primary, mt_names):
        if any(is_mapped_to_mt(read, mt_names) for read in records if read.is_secondary or read.is_supplementary):
            return "nuclear_primary_with_mt_competitor", primary
        return "nuclear_primary", primary
    if primary.mapping_quality < min_mt_mapq or aligned_fraction(primary) < min_mt_aligned_fraction:
        return "weak_mt_primary", primary
    has_nuclear_competitor = any(
        is_strong_competitor(read, mt_names, competing_nuclear_aligned_fraction)
        for read in records
        if read.is_secondary or read.is_supplementary
    )
    if primary.mapping_quality < ambiguous_mapq_below or has_nuclear_competitor:
        return "mt_primary_ambiguous", primary
    return "mt_primary_best", primary


def process_group(records: list[pysam.AlignedSegment], args: argparse.Namespace, mt_names: set[str], out_handle, written: set[tuple[str, str]], counts: Counter) -> None:
    if not records:
        return
    category, primary = classify_group(
        records,
        mt_names,
        args.min_mt_mapq,
        args.min_mt_aligned_fraction,
        args.ambiguous_mapq_below,
        args.competing_nuclear_aligned_fraction,
    )
    counts[category] += 1
    if primary is None:
        return
    selected = category == "mt_primary_best" or (category == "mt_primary_ambiguous" and args.keep_ambiguous_mt_nuclear)
    if not selected:
        return
    key = read_key(primary)
    if key in written:
        return
    out_handle.write(fastq_entry(primary))
    written.add(key)
    counts["retained_for_mitochondrial_remap"] += 1
    counts[f"{category}_written"] += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--mt-evidence-fastq", required=True)
    parser.add_argument("--high-confidence-fastq", required=True)
    parser.add_argument("--ambiguous-fastq", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-mt-mapq", type=int, default=0)
    parser.add_argument("--min-mt-aligned-fraction", type=float, default=0.5)
    parser.add_argument("--ambiguous-mapq-below", type=int, default=10)
    parser.add_argument("--competing-nuclear-aligned-fraction", type=float, default=0.5)
    parser.add_argument("--keep-ambiguous-mt-nuclear", action="store_true")
    args = parser.parse_args()

    mt_names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    counts = Counter()
    written: set[tuple[str, str]] = set()
    current_key: tuple[str, str] | None = None
    current_records: list[pysam.AlignedSegment] = []
    ensure_parent(args.mt_evidence_fastq)
    with pysam.AlignmentFile("-", "r") as sam, gzip.open(args.mt_evidence_fastq, "wt", encoding="utf-8") as out_handle:
        for read in sam.fetch(until_eof=True):
            counts["alignment_records_examined"] += 1
            key = read_key(read)
            if current_key is None:
                current_key = key
            if key != current_key:
                process_group(current_records, args, mt_names, out_handle, written, counts)
                current_key = key
                current_records = []
            current_records.append(read)
        process_group(current_records, args, mt_names, out_handle, written, counts)

    empty_gzip(args.high_confidence_fastq)
    empty_gzip(args.ambiguous_fastq)
    write_classification(args.classification, counts)
    write_json(
        args.summary,
        {
            "sample": args.sample,
            "selection_source": "whole_genome_competitive_alignment",
            "selection_strategy": "whole_genome_mt_best",
            "mt_contig_names": sorted(mt_names),
            "min_mt_mapq": args.min_mt_mapq,
            "min_mt_aligned_fraction": args.min_mt_aligned_fraction,
            "ambiguous_mapq_below": args.ambiguous_mapq_below,
            "competing_nuclear_aligned_fraction": args.competing_nuclear_aligned_fraction,
            "keep_ambiguous_mt_nuclear": args.keep_ambiguous_mt_nuclear,
            "total_alignments_examined": counts["alignment_records_examined"],
            "mt_primary_best_reads": counts["mt_primary_best"],
            "mt_primary_ambiguous_reads": counts["mt_primary_ambiguous"],
            "nuclear_primary_reads": counts["nuclear_primary"],
            "nuclear_primary_with_mt_competitor_reads": counts["nuclear_primary_with_mt_competitor"],
            "weak_mt_primary_reads": counts["weak_mt_primary"],
            "primary_unmapped_reads": counts["primary_unmapped"],
            "mt_evidence_reads_selected": counts["retained_for_mitochondrial_remap"],
            "mt_evidence_fastq_records_written": counts["retained_for_mitochondrial_remap"],
            "high_confidence_mt": 0,
            "ambiguous_mt": 0,
        },
    )


if __name__ == "__main__":
    main()
