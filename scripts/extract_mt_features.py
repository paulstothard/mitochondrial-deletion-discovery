#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from common import write_tsv


def parse_attributes(text: str) -> dict[str, str]:
    attrs = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text):
        attrs[key] = value
    return attrs


def open_text_maybe_gzip(path: str):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b" or Path(path).suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mt_names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    rows = []
    with open_text_maybe_gzip(args.gtf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[0] not in mt_names:
                continue
            feature_type = parts[2]
            if feature_type not in {"gene", "transcript", "exon", "tRNA", "rRNA", "CDS"}:
                continue
            attrs = parse_attributes(parts[8])
            gene_name = attrs.get("gene_name") or attrs.get("gene_id") or attrs.get("transcript_name") or attrs.get("transcript_id") or ""
            rows.append(
                {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "feature_type": feature_type,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": gene_name,
                    "transcript_id": attrs.get("transcript_id", ""),
                    "product": attrs.get("product", ""),
                }
            )
    rows.sort(key=lambda row: (row["start"], row["end"], row["feature_type"]))
    write_tsv(
        args.output,
        rows,
        fieldnames=["contig", "start", "end", "strand", "feature_type", "gene_id", "gene_name", "transcript_id", "product"],
    )


if __name__ == "__main__":
    main()
