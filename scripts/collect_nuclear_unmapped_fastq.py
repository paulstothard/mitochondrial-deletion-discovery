#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
from pathlib import Path

from common import empty_gzip, ensure_parent, write_json


def read_fastq_records(path: str | Path):
    if not str(path):
        return
    path = Path(path)
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                break
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def existing_fastqs(paths: list[str]) -> list[Path]:
    out = []
    for item in paths:
        if not item:
            continue
        path = Path(item)
        if path.exists() and path.is_file() and path.stat().st_size > 0:
            out.append(path)
    return out


def count_fastq_records(path: Path) -> int:
    lines = int(subprocess.check_output(["wc", "-l", str(path)], text=True).split()[0])
    return lines // 4


def gzip_concat_fastqs(paths: list[Path], output: str, threads: int) -> None:
    ensure_parent(output)
    with open(output, "wb") as out_handle:
        pigz = subprocess.Popen(["pigz", "-p", str(max(1, threads)), "-c"], stdin=subprocess.PIPE, stdout=out_handle)
        assert pigz.stdin is not None
        with pigz.stdin:
            for path in paths:
                with open(path, "rb") as in_handle:
                    shutil.copyfileobj(in_handle, pigz.stdin, length=16 * 1024 * 1024)
        status = pigz.wait()
    if status != 0:
        raise subprocess.CalledProcessError(status, ["pigz", "-c"])


def write_classification(path: str, counts: dict[str, int]) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["table_type", "name", "count"])
        for key, value in counts.items():
            writer.writerow(["selection_reason", key, value])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--mate1", required=True)
    parser.add_argument("--mate2", default="")
    parser.add_argument("--mt-evidence-fastq", required=True)
    parser.add_argument("--high-confidence-fastq", required=True)
    parser.add_argument("--ambiguous-fastq", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    mate1 = Path(args.mate1) if args.mate1 else None
    mate2 = Path(args.mate2) if args.mate2 else None
    mate1_written = count_fastq_records(mate1) if mate1 and mate1.exists() and mate1.is_file() else 0
    mate2_written = count_fastq_records(mate2) if mate2 and mate2.exists() and mate2.is_file() else 0
    records_written = mate1_written + mate2_written
    fastqs = existing_fastqs([args.mate1, args.mate2])
    if fastqs:
        gzip_concat_fastqs(fastqs, args.mt_evidence_fastq, args.threads)
    else:
        empty_gzip(args.mt_evidence_fastq)

    empty_gzip(args.high_confidence_fastq)
    empty_gzip(args.ambiguous_fastq)
    write_classification(args.classification, {"nuclear_unmapped_read": records_written})
    write_json(
        args.summary,
        {
            "sample": args.sample,
            "selection_source": args.source,
            "selection_strategy": "nuclear_unmapped_reads",
            "mt_evidence_reads_selected": records_written,
            "mt_evidence_fastq_records_written": records_written,
            "nuclear_unmapped_records_written": records_written,
            "nuclear_unmapped_mate1_records_written": mate1_written,
            "nuclear_unmapped_mate2_records_written": mate2_written,
            "high_confidence_mt": 0,
            "ambiguous_mt": 0,
        },
    )


if __name__ == "__main__":
    main()
