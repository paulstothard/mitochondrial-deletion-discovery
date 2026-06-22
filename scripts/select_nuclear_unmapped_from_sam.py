#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter

import pysam

from classify_mt_reads import fastq_entry, read_key
from common import empty_gzip, ensure_parent, write_json


def write_classification(path: str, counts: Counter) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["table_type", "name", "count"])
        for name, count in sorted(counts.items()):
            writer.writerow(["alignment_category", name, count])
        writer.writerow(["selection_reason", "retained_for_mitochondrial_remap", counts["retained_for_mitochondrial_remap"]])


def query_aligned_fraction(read: pysam.AlignedSegment) -> float:
    query_length = read.query_length or len(read.query_sequence or "")
    if not query_length:
        return 0.0
    return float(read.query_alignment_length or 0) / float(query_length)


def is_strong_nuclear_alignment(read: pysam.AlignedSegment, min_mapq: int, min_aligned_fraction: float) -> bool:
    return read.mapping_quality >= min_mapq and query_aligned_fraction(read) >= min_aligned_fraction


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--mt-evidence-fastq", required=True)
    parser.add_argument("--high-confidence-fastq", required=True)
    parser.add_argument("--ambiguous-fastq", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-nuclear-mapq", type=int, default=20)
    parser.add_argument("--min-nuclear-aligned-fraction", type=float, default=0.8)
    args = parser.parse_args()

    counts = Counter()
    written = set()
    ensure_parent(args.mt_evidence_fastq)
    with pysam.AlignmentFile("-", "r") as sam, gzip.open(args.mt_evidence_fastq, "wt", encoding="utf-8") as out_handle:
        for read in sam.fetch(until_eof=True):
            counts["alignment_records_examined"] += 1
            if read.is_secondary:
                counts["secondary_records"] += 1
                continue
            if read.is_supplementary:
                counts["supplementary_records"] += 1
                continue
            if read.is_unmapped:
                counts["primary_unmapped_records"] += 1
                retention_reason = "primary_unmapped"
            elif is_strong_nuclear_alignment(read, args.min_nuclear_mapq, args.min_nuclear_aligned_fraction):
                counts["strong_primary_nuclear_mapped_records"] += 1
                continue
            else:
                counts["weak_or_partial_primary_nuclear_mapped_records"] += 1
                retention_reason = "weak_or_partial_nuclear_alignment"

            if retention_reason:
                key = read_key(read)
                if key not in written:
                    out_handle.write(fastq_entry(read))
                    written.add(key)
                    counts["retained_for_mitochondrial_remap"] += 1
                    counts[f"{retention_reason}_written"] += 1

    empty_gzip(args.high_confidence_fastq)
    empty_gzip(args.ambiguous_fastq)
    write_classification(args.classification, counts)
    write_json(
        args.summary,
        {
            "sample": args.sample,
            "selection_source": args.source,
            "selection_strategy": "nuclear_unmapped_reads",
            "strong_nuclear_alignment_min_mapq": args.min_nuclear_mapq,
            "strong_nuclear_alignment_min_aligned_fraction": args.min_nuclear_aligned_fraction,
            "total_alignments_examined": counts["alignment_records_examined"],
            "strong_primary_nuclear_mapped_records": counts["strong_primary_nuclear_mapped_records"],
            "weak_or_partial_primary_nuclear_mapped_records": counts["weak_or_partial_primary_nuclear_mapped_records"],
            "nuclear_unmapped_records_written": counts["primary_unmapped_written"],
            "retained_for_mitochondrial_remap": counts["retained_for_mitochondrial_remap"],
            "mt_evidence_reads_selected": counts["retained_for_mitochondrial_remap"],
            "mt_evidence_fastq_records_written": counts["retained_for_mitochondrial_remap"],
            "high_confidence_mt": 0,
            "ambiguous_mt": 0,
        },
    )


if __name__ == "__main__":
    main()
