#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter

import pysam

from common import ensure_parent, write_json


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def fastq_entry(read: pysam.AlignedSegment) -> str:
    seq = read.query_sequence or ""
    qual = read.qual or "I" * len(seq)
    if read.is_reverse:
        seq = reverse_complement(seq)
        qual = qual[::-1]
    suffix = "/2" if read.is_read2 else "/1" if read.is_read1 else ""
    return f"@{read.query_name}{suffix}\n{seq}\n+\n{qual}\n"


def nh_tag(read: pysam.AlignedSegment) -> int:
    return int(read.get_tag("NH")) if read.has_tag("NH") else 1


def classify_alignment(read: pysam.AlignedSegment, mt_names: set[str], min_mapq: int) -> str:
    if getattr(read, "is_unmapped", False):
        return "unmapped"
    if getattr(read, "is_supplementary", False):
        return "supplementary_mt" if read.reference_name in mt_names else "supplementary_non_mt"
    if getattr(read, "is_secondary", False):
        return "secondary_mt" if read.reference_name in mt_names else "secondary_non_mt"
    if read.reference_name not in mt_names:
        return "non_mt_primary"
    if read.mapping_quality < min_mapq:
        return "low_quality_mt"
    if nh_tag(read) > 1:
        return "ambiguous_mt"
    return "high_confidence_mt"


def read_key(read: pysam.AlignedSegment) -> tuple[str, str]:
    mate = "2" if getattr(read, "is_read2", False) else "1" if getattr(read, "is_read1", False) else ""
    return read.query_name, mate


def clean_read_name(name: str) -> str:
    return name.removesuffix("/1").removesuffix("/2")


def read_names_from_chimeric_file(path: str, mt_names: set[str]) -> set[str]:
    names = set()
    if not path:
        return names
    try:
        handle = open(path, "r", encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return names
    with handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                parts = line.split()
            if len(parts) < 10:
                continue
            if parts[0] in mt_names or parts[3] in mt_names:
                names.add(clean_read_name(parts[9]))
    return names


def category_is_selected(category: str, args: argparse.Namespace) -> bool:
    if category == "high_confidence_mt":
        return True
    if category == "low_quality_mt":
        return args.include_low_mapq
    if category == "ambiguous_mt":
        return args.include_multimappers
    if category == "supplementary_mt":
        return args.include_supplementary
    if category == "secondary_mt":
        return args.include_secondary
    return False


def write_classification_row(
    writer: csv.DictWriter,
    read: pysam.AlignedSegment,
    category: str,
    selected: bool,
    reasons: set[str],
) -> None:
    writer.writerow(
        {
            "read_id": read.query_name,
            "mate": "2" if read.is_read2 else "1" if read.is_read1 else "",
            "category": category,
            "reference": "" if read.is_unmapped else read.reference_name,
            "mapq": read.mapping_quality,
            "nh": "" if read.is_unmapped else nh_tag(read),
            "selected_for_mt_remap": "yes" if selected else "no",
            "selection_reasons": ";".join(sorted(reasons)),
        }
    )


def selected_for_read(
    read: pysam.AlignedSegment,
    category: str,
    chimeric_names: set[str],
    args: argparse.Namespace,
    mate_names: set[str] | None = None,
) -> tuple[bool, set[str]]:
    reasons = set()
    if category_is_selected(category, args):
        reasons.add(category)
    if clean_read_name(read.query_name) in chimeric_names:
        reasons.add("full_genome_chimeric_mt_evidence")
    if mate_names is not None and read.query_name in mate_names:
        reasons.add("mate_of_mt_evidence_read")
    return bool(reasons), reasons


def write_summary_table(path: str, counts: Counter, selected_counts: Counter, reason_counts: Counter) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["table_type", "name", "count"])
        for name, count in sorted(counts.items()):
            writer.writerow(["alignment_category", name, count])
        for name, count in sorted(selected_counts.items()):
            writer.writerow(["selected_alignment_category", name, count])
        for name, count in sorted(reason_counts.items()):
            writer.writerow(["selection_reason", name, count])


def find_selected_names(
    bam_path: str,
    mt_names: set[str],
    min_mapq: int,
    chimeric_names: set[str],
    args: argparse.Namespace,
) -> set[str]:
    selected_names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            category = classify_alignment(read, mt_names, min_mapq)
            selected, _ = selected_for_read(read, category, chimeric_names, args)
            if selected:
                selected_names.add(read.query_name)
    return selected_names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--chimeric-junction", default="")
    parser.add_argument("--high-confidence-fastq", required=True)
    parser.add_argument("--ambiguous-fastq", required=True)
    parser.add_argument("--mt-evidence-fastq", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--include-low-mapq", action="store_true")
    parser.add_argument("--include-multimappers", action="store_true")
    parser.add_argument("--include-supplementary", action="store_true")
    parser.add_argument("--include-secondary", action="store_true")
    parser.add_argument("--include-chimeric-mt-reads", action="store_true")
    parser.add_argument("--include-mates-of-mt-evidence-reads", action="store_true")
    parser.add_argument("--write-read-classification-tsv", action="store_true")
    args = parser.parse_args()

    mt_names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    chimeric_names = read_names_from_chimeric_file(args.chimeric_junction, mt_names) if args.include_chimeric_mt_reads else set()
    counts = Counter()
    selected_counts = Counter()
    reason_counts = Counter()
    selected_names = (
        find_selected_names(args.bam, mt_names, args.min_mapq, chimeric_names, args)
        if args.include_mates_of_mt_evidence_reads
        else None
    )

    ensure_parent(args.high_confidence_fastq)
    ensure_parent(args.classification)
    fieldnames = [
        "read_id",
        "mate",
        "category",
        "reference",
        "mapq",
        "nh",
        "selected_for_mt_remap",
        "selection_reasons",
    ]
    written_evidence = set()
    written_high = set()
    written_ambiguous = set()
    selected_keys = set()
    with pysam.AlignmentFile(args.bam, "rb") as bam, gzip.open(
        args.high_confidence_fastq, "wt", encoding="utf-8"
    ) as high, gzip.open(args.ambiguous_fastq, "wt", encoding="utf-8") as ambiguous, gzip.open(
        args.mt_evidence_fastq, "wt", encoding="utf-8"
    ) as evidence:
        class_handle = open(args.classification, "w", encoding="utf-8", newline="") if args.write_read_classification_tsv else None
        writer = csv.DictWriter(class_handle, delimiter="\t", fieldnames=fieldnames) if class_handle else None
        if writer:
            writer.writeheader()
        for read in bam.fetch(until_eof=True):
            key = read_key(read)
            category = classify_alignment(read, mt_names, args.min_mapq)
            counts[category] += 1
            selected, reasons = selected_for_read(read, category, chimeric_names, args, selected_names)
            if selected:
                selected_keys.add(key)
                selected_counts[category] += 1
                for reason in reasons:
                    reason_counts[reason] += 1
            is_primary_record = not read.is_secondary and not read.is_supplementary
            if category == "high_confidence_mt" and key not in written_high and is_primary_record:
                high.write(fastq_entry(read))
                written_high.add(key)
            if category in {"ambiguous_mt", "low_quality_mt"} and key not in written_ambiguous and is_primary_record:
                ambiguous.write(fastq_entry(read))
                written_ambiguous.add(key)
            if selected and key not in written_evidence and is_primary_record:
                evidence.write(fastq_entry(read))
                written_evidence.add(key)
            if writer and (selected or category.endswith("_mt") or category in {"ambiguous_mt", "low_quality_mt", "high_confidence_mt"}):
                write_classification_row(writer, read, category, selected, reasons)
        if class_handle:
            class_handle.close()

    if not args.write_read_classification_tsv:
        write_summary_table(args.classification, counts, selected_counts, reason_counts)

    write_json(
        args.summary,
        {
            "total_alignments_examined": sum(counts.values()),
            "mt_evidence_reads_selected": len(selected_keys),
            "mt_evidence_fastq_records_written": len(written_evidence),
            "chimeric_mt_read_names": len(chimeric_names),
            "selection_includes_low_mapq": args.include_low_mapq,
            "selection_includes_multimappers": args.include_multimappers,
            "selection_includes_supplementary": args.include_supplementary,
            "selection_includes_secondary": args.include_secondary,
            "selection_includes_full_genome_chimeric_mt_reads": args.include_chimeric_mt_reads,
            "selection_includes_mates": args.include_mates_of_mt_evidence_reads,
            "read_classification_tsv_mode": "per_read" if args.write_read_classification_tsv else "summary",
            **dict(counts),
        },
    )


if __name__ == "__main__":
    main()
