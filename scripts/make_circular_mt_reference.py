#!/usr/bin/env python3
from __future__ import annotations

import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from common import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--padding", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    record = next(SeqIO.parse(args.input, "fasta"))
    seq = record.seq
    padding = min(args.padding, len(seq))
    circular = seq[-padding:] + seq + seq[:padding]
    out = SeqRecord(circular, id=f"{record.id}_circular_padded", description=f"padding={padding};source={record.id}")
    ensure_parent(args.output)
    SeqIO.write(out, args.output, "fasta")


if __name__ == "__main__":
    main()
