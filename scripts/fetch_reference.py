#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
from pathlib import Path

from common import ensure_parent


def open_maybe_gzip(path: Path):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ensure_parent(args.output)
    output = Path(args.output)
    tmp = output.with_suffix(output.suffix + ".tmp")
    source = args.source
    downloaded = None
    if source.startswith(("http://", "https://", "ftp://")):
        downloaded = output.with_suffix(output.suffix + ".download")
        urllib.request.urlretrieve(source, downloaded)
        source_path = downloaded
    else:
        source_path = Path(source)
    with open_maybe_gzip(source_path) as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp.replace(output)
    if downloaded and downloaded.exists():
        downloaded.unlink()


if __name__ == "__main__":
    main()
