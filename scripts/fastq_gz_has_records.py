#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys

from common import gzip_is_nonempty


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_gz")
    args = parser.parse_args()
    if gzip_is_nonempty(args.fastq_gz):
        return
    sys.exit(1)


if __name__ == "__main__":
    main()
