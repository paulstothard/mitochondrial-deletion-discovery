#!/usr/bin/env python3
from __future__ import annotations

import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from common import ensure_parent, write_json


def rotate_sequence(sequence: str, start: int) -> str:
    if start < 1 or start > len(sequence):
        raise ValueError(f"rotation start must be within 1..{len(sequence)}")
    offset = start - 1
    return sequence[offset:] + sequence[:offset]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))
    if len(records) != 1:
        raise ValueError(f"expected one mitochondrial FASTA record, found {len(records)}")
    source = records[0]
    rotated = rotate_sequence(str(source.seq), args.start)
    ensure_parent(args.output)
    SeqIO.write(
        [SeqRecord(Seq(rotated), id=f"{source.id}_{args.name}", description=f"rotated_start_{args.start}")],
        args.output,
        "fasta",
    )
    write_json(
        args.metadata,
        {
            "source_record": source.id,
            "rotation_name": args.name,
            "rotation_start": args.start,
            "mt_length": len(source.seq),
            "duplicated_sequence": False,
        },
    )


if __name__ == "__main__":
    main()
