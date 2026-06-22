#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

from Bio import SeqIO

from common import ensure_parent, write_json


def open_text_maybe_gzip(path: str):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b" or Path(path).suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--expected-length", type=int, required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args()

    names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    chosen = None
    with open_text_maybe_gzip(args.genome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_names = {record.id, record.name, record.description.split()[0]}
            if names & record_names:
                chosen = record
                break
    if chosen is None:
        raise SystemExit(f"No mitochondrial contig matching {sorted(names)} in {args.genome}")
    ensure_parent(args.out_fasta)
    SeqIO.write(chosen, args.out_fasta, "fasta")
    write_json(
        args.out_json,
        {
            "contig": chosen.id,
            "length": len(chosen.seq),
            "expected_length": args.expected_length,
            "length_matches_expected": len(chosen.seq) == args.expected_length,
        },
    )


if __name__ == "__main__":
    main()
