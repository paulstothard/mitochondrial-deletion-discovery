#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import SeqIO

from common import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args()

    mt_names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    kept = []
    removed = []

    ensure_parent(args.out_fasta)
    with open(args.out_fasta, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(args.genome, "fasta"):
            row = {"id": record.id, "length": len(record.seq)}
            if record.id in mt_names or record.name in mt_names:
                removed.append(row)
                continue
            kept.append(row)
            SeqIO.write(record, out_handle, "fasta")

    ensure_parent(args.out_json)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "source_genome": str(Path(args.genome)),
                "mt_contig_names": sorted(mt_names),
                "nuclear_contigs_written": len(kept),
                "nuclear_bases_written": sum(item["length"] for item in kept),
                "excluded_contigs": removed,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    if not kept:
        raise SystemExit("No nuclear contigs were written; check mt_contig_names and genome FASTA.")


if __name__ == "__main__":
    main()
