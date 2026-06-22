from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path



def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    ensure_parent(path)
    if fieldnames is None:
        fieldnames = sorted({key for row in rows for key in row})
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_yaml(path: str | Path) -> dict:
    import yaml

    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def write_yaml(path: str | Path, data: dict) -> None:
    import yaml

    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)


def read_json(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: str | Path, data: dict) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def deep_update(base: dict, override: dict) -> dict:
    merged = dict(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_update(merged[key], value)
        else:
            merged[key] = value
    return merged


def run(cmd: list[str], cwd: str | Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def gzip_is_nonempty(path: str | Path) -> bool:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            return bool(handle.readline())
    except gzip.BadGzipFile:
        return False


def empty_gzip(path: str | Path) -> None:
    ensure_parent(path)
    with gzip.open(path, "wb"):
        pass


def copy_or_link(src: str | Path, dst: str | Path) -> None:
    ensure_parent(dst)
    dst = Path(dst)
    if dst.exists():
        dst.unlink()
    try:
        os.symlink(Path(src).resolve(), dst)
    except OSError:
        shutil.copyfile(src, dst)


def fastq_record_count_gz(path: str | Path) -> int:
    if not gzip_is_nonempty(path):
        return 0
    lines = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for lines, _ in enumerate(handle, 1):
            pass
    return lines // 4
