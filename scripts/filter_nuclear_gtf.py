#!/usr/bin/env python3
from __future__ import annotations

import argparse

from common import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--out-gtf", required=True)
    args = parser.parse_args()

    mt_names = {name.strip() for name in args.mt_contig_names.split(",") if name.strip()}
    kept = 0
    skipped = 0
    ensure_parent(args.out_gtf)
    with open(args.gtf, "r", encoding="utf-8", errors="replace") as in_handle, open(args.out_gtf, "w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if line.startswith("#"):
                out_handle.write(line)
                continue
            contig = line.split("\t", 1)[0]
            if contig in mt_names:
                skipped += 1
                continue
            out_handle.write(line)
            kept += 1
    if kept == 0:
        raise SystemExit(f"No non-mitochondrial GTF rows were written from {args.gtf}; check mt_contig_names.")


if __name__ == "__main__":
    main()
