#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from common import copy_or_link, empty_gzip, ensure_parent, fastq_record_count_gz, gzip_is_nonempty, write_json, write_tsv


def fragment_count_from_fastp_json(path: str, paired: bool) -> int:
    with open(path, "r", encoding="utf-8") as handle:
        report = json.load(handle)
    if paired and "read1_after_filtering" in report:
        return int(report["read1_after_filtering"]["total_reads"])
    return int(report["summary"]["after_filtering"]["total_reads"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--in-r1", required=True)
    parser.add_argument("--in-r2", required=True)
    parser.add_argument("--out-r1", required=True)
    parser.add_argument("--out-r2", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--html", required=True)
    parser.add_argument("--decision", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--extra", default="")
    parser.add_argument("--skip", action="store_true")
    args = parser.parse_args()

    paired = gzip_is_nonempty(args.in_r2)
    if args.skip:
        copy_or_link(args.in_r1, args.out_r1)
        if paired:
            copy_or_link(args.in_r2, args.out_r2)
        else:
            empty_gzip(args.out_r2)
        fragments = fastq_record_count_gz(args.out_r1)
        write_json(
            args.json,
            {
                "summary": {
                    "before_filtering": {"total_reads": fragments},
                    "after_filtering": {"total_reads": fragments},
                },
                "read1_after_filtering": {"total_reads": fragments},
            },
        )
        ensure_parent(args.html)
        Path(args.html).write_text(
            f"<html><body><h1>{args.sample}</h1><p>Trimming was disabled; input reads were passed through.</p></body></html>\n",
            encoding="utf-8",
        )
        write_json(args.decision, {"sample": args.sample, "trimmed": False, "paired": paired, "minimum_length": args.min_length})
        write_tsv(
            args.counts,
            [{"sample": args.sample, "fragments": fragments, "million_fragments": fragments / 1_000_000}],
            fieldnames=["sample", "fragments", "million_fragments"],
        )
        return
    cmd = [
        "fastp",
        "-i",
        args.in_r1,
        "-o",
        args.out_r1,
        "-j",
        args.json,
        "-h",
        args.html,
        "--length_required",
        str(args.min_length),
        "--thread",
        str(args.threads),
    ]
    if paired:
        cmd.extend(["-I", args.in_r2, "-O", args.out_r2])
    for part in shlex.split(args.extra):
        cmd.append(part)
    subprocess.run(cmd, check=True)
    if not paired:
        empty_gzip(args.out_r2)
    fragments = fragment_count_from_fastp_json(args.json, paired)
    write_json(args.decision, {"sample": args.sample, "trimmed": True, "paired": paired, "minimum_length": args.min_length})
    write_tsv(
        args.counts,
        [{"sample": args.sample, "fragments": fragments, "million_fragments": fragments / 1_000_000}],
        fieldnames=["sample", "fragments", "million_fragments"],
    )


if __name__ == "__main__":
    main()
