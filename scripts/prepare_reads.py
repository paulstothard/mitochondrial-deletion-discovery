#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

from common import deep_update, empty_gzip, ensure_parent, gzip_is_nonempty, read_tsv, read_yaml, write_json


def find_row(sample_table: str, sample: str) -> dict[str, str]:
    for row in read_tsv(sample_table):
        if row["sample"] == sample:
            return row
    raise SystemExit(f"Sample {sample!r} not found in {sample_table}")


def log_line(log_handle, message: str) -> None:
    print(message, file=log_handle, flush=True)


def run_logged(cmd: list[str], log_handle, cwd: str | Path | None = None, env: dict | None = None) -> None:
    log_line(log_handle, "$ " + " ".join(str(part) for part in cmd))
    result = subprocess.run(cmd, cwd=cwd, env=env, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")


def is_gzip_file(path: str | Path) -> bool:
    try:
        with open(path, "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def copy_or_link(src: str | Path, dst: str | Path) -> None:
    ensure_parent(dst)
    dst = Path(dst)
    if dst.exists():
        dst.unlink()
    try:
        dst.symlink_to(Path(src).resolve())
    except OSError:
        shutil.copyfile(src, dst)


def gzip_copy(src: str | Path, dst: str | Path, log_handle) -> None:
    ensure_parent(dst)
    tmp = Path(str(dst) + ".tmp")
    pigz = shutil.which("pigz")
    if pigz:
        log_line(log_handle, f"$ {pigz} -c {src} > {tmp}")
        with open(tmp, "wb") as out_handle:
            result = subprocess.run([pigz, "-c", str(src)], stdout=out_handle, stderr=log_handle, text=False)
        if result.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"pigz failed for {src}")
    else:
        log_line(log_handle, f"Compressing local FASTQ with Python gzip: {src}")
        with open(src, "rb") as in_handle, gzip.open(tmp, "wb") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)
    tmp.replace(dst)


def stage_fastq_source(src: str, dst: str, log_handle) -> None:
    if src.startswith(("http://", "https://", "ftp://")):
        stage_url(src, dst, log_handle)
        return
    src_path = Path(src)
    if is_gzip_file(src_path):
        copy_or_link(src_path, dst)
    else:
        gzip_copy(src_path, dst, log_handle)


def stage_url(src: str, dst: str, log_handle) -> None:
    ensure_parent(dst)
    tmp = str(dst) + ".download"
    curl = shutil.which("curl")
    if curl:
        run_logged(
            [
                curl,
                "-L",
                "--fail",
                "--retry",
                "5",
                "--retry-delay",
                "10",
                "-C",
                "-",
                "-o",
                tmp,
                src,
            ],
            log_handle,
        )
    else:
        log_line(log_handle, f"Downloading with urllib: {src}")
        urllib.request.urlretrieve(src, tmp)
    if is_gzip_file(tmp):
        shutil.move(tmp, dst)
    else:
        gzip_copy(tmp, dst, log_handle)
        Path(tmp).unlink(missing_ok=True)


def validate_fastq_gz(path: str, allow_empty: bool = False) -> None:
    p = Path(path)
    if allow_empty and p.exists():
        try:
            with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as handle:
                if handle.readline() == "":
                    return
        except gzip.BadGzipFile:
            pass
    if not gzip_is_nonempty(path):
        raise RuntimeError(f"FASTQ gzip validation failed for {path}")
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        first = handle.readline()
    if not first.startswith("@"):
        raise RuntimeError(f"FASTQ validation failed for {path}: first record does not start with @")


def existing_outputs_are_usable(out_r1: str, out_r2: str, validate: bool) -> bool:
    r1 = Path(out_r1)
    r2 = Path(out_r2)
    if not r1.exists() or not r2.exists():
        return False
    if validate:
        validate_fastq_gz(out_r1)
        validate_fastq_gz(out_r2, allow_empty=True)
    return True


def fastq_has_sequence(path: str | Path, max_records: int = 1000) -> bool:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for _ in range(max_records):
            header = handle.readline()
            if not header:
                return False
            sequence = handle.readline().strip()
            handle.readline()
            handle.readline()
            if sequence:
                return True
    return False


def stage_fastq_url(row: dict[str, str], out_r1: str, out_r2: str, log_handle) -> dict:
    fastq_1 = row.get("fastq_1", "")
    fastq_2 = row.get("fastq_2", "")
    declared_layout = normalized_layout(row)
    if not fastq_1:
        raise RuntimeError("No fastq_1 URL/path is available for ena_fastq download")
    stage_fastq_source(fastq_1, out_r1, log_handle)
    if fastq_2:
        stage_fastq_source(fastq_2, out_r2, log_handle)
        paired = True
    else:
        empty_gzip(out_r2)
        paired = False
    if declared_layout == "paired" and not paired:
        log_line(log_handle, "metadata declares paired layout but no fastq_2 is available; treating as single-end")
    return {"method": "ena_fastq", "paired": paired, "declared_layout": declared_layout, "fastq_1": fastq_1, "fastq_2": fastq_2}


def normalized_layout(row: dict[str, str]) -> str:
    layout = (row.get("layout") or row.get("LibraryLayout") or "").strip().lower()
    if layout in {"paired", "paired-end", "paired_end", "pe"}:
        return "paired"
    return "single"


def stage_fasterq_dump(
    row: dict[str, str],
    out_r1: str,
    out_r2: str,
    threads: int,
    prefetch: bool,
    sra_settings_dir: str,
    log_handle,
) -> dict:
    accession = row.get("run_accession", "")
    if not accession:
        raise RuntimeError("No run_accession is available for fasterq_dump")
    declared_layout = normalized_layout(row)
    outdir = Path(out_r1).parent
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = outdir / f".{accession}.fasterq"
    tempdir = workdir / "tmp"
    workdir.mkdir(parents=True, exist_ok=True)
    tempdir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if sra_settings_dir:
        env["NCBI_HOME"] = str(Path(sra_settings_dir).resolve())
        log_line(log_handle, f"NCBI_HOME={env['NCBI_HOME']}")
    if prefetch:
        run_logged(["prefetch", accession, "--output-directory", str(workdir)], log_handle, env=env)
    cmd = [
        "fasterq-dump",
        accession,
        "--outdir",
        str(workdir),
        "--threads",
        str(threads),
        "--temp",
        str(tempdir),
        "--progress",
    ]
    if declared_layout == "paired":
        cmd.insert(2, "--split-files")
    run_logged(cmd, log_handle, env=env)
    r1 = workdir / f"{accession}_1.fastq"
    r2 = workdir / f"{accession}_2.fastq"
    single = workdir / f"{accession}.fastq"
    if declared_layout == "single" and single.exists():
        run_logged(["pigz", "-f", str(single)], log_handle)
        shutil.move(str(single) + ".gz", out_r1)
        empty_gzip(out_r2)
        paired = False
    elif r1.exists():
        run_logged(["pigz", "-f", str(r1)], log_handle)
        shutil.move(str(r1) + ".gz", out_r1)
        if r2.exists():
            if fastq_has_sequence(r2):
                run_logged(["pigz", "-f", str(r2)], log_handle)
                shutil.move(str(r2) + ".gz", out_r2)
                paired = True
            else:
                log_line(log_handle, f"{r2} contains only zero-length reads; treating run as single-end")
                r2.unlink()
                empty_gzip(out_r2)
                paired = False
        else:
            empty_gzip(out_r2)
            paired = False
    elif single.exists():
        run_logged(["pigz", "-f", str(single)], log_handle)
        shutil.move(str(single) + ".gz", out_r1)
        empty_gzip(out_r2)
        paired = False
    else:
        raise RuntimeError(f"fasterq-dump produced no FASTQ for {accession} in {workdir}")
    shutil.rmtree(workdir, ignore_errors=True)
    return {"method": "fasterq_dump", "paired": paired, "declared_layout": declared_layout, "run_accession": accession}


def method_order(requested: str, config: dict) -> list[str]:
    if requested == "auto":
        return list(config.get("downloads", {}).get("auto_order", ["fasterq_dump", "ena_fastq"]))
    if requested in {"fasterq_dump", "ena_fastq"}:
        return [requested]
    raise SystemExit(f"Unsupported download method {requested!r}; use auto, fasterq_dump, or ena_fastq")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sample-table", required=True)
    parser.add_argument("--defaults", default="")
    parser.add_argument("--config", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--fasterq-threads", type=int, default=2)
    parser.add_argument("--prefetch", default="False")
    parser.add_argument("--out-r1", required=True)
    parser.add_argument("--out-r2", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    row = find_row(args.sample_table, args.sample)
    config = deep_update(read_yaml(args.defaults), read_yaml(args.config)) if args.defaults else read_yaml(args.config)
    validate = bool(config.get("downloads", {}).get("validate_fastq", True))
    sra_settings_dir = config.get("downloads", {}).get("sra_settings_dir", "")
    prefetch = str(args.prefetch).lower() in {"1", "true", "yes"}
    ensure_parent(args.log)
    attempts = []
    with open(args.log, "w", encoding="utf-8") as log_handle:
        log_line(log_handle, f"sample={args.sample}")
        log_line(log_handle, f"requested_method={args.method}")
        if existing_outputs_are_usable(args.out_r1, args.out_r2, validate):
            paired = gzip_is_nonempty(args.out_r2)
            summary = {
                "sample": args.sample,
                "method": "existing_fastq",
                "requested_method": args.method,
                "paired": paired,
                "declared_layout": normalized_layout(row),
                "reused_existing_outputs": True,
                "attempts": [{"method": "existing_fastq", "status": "success"}],
            }
            run_accession = row.get("run_accession", "")
            if run_accession:
                summary["run_accession"] = run_accession
            write_json(args.summary, summary)
            log_line(log_handle, "success_method=existing_fastq")
            return
        for method in method_order(args.method, config):
            try:
                log_line(log_handle, f"attempt_method={method}")
                if method == "fasterq_dump":
                    summary = stage_fasterq_dump(
                        row,
                        args.out_r1,
                        args.out_r2,
                        args.fasterq_threads,
                        prefetch,
                        sra_settings_dir,
                        log_handle,
                    )
                elif method == "ena_fastq":
                    summary = stage_fastq_url(row, args.out_r1, args.out_r2, log_handle)
                else:
                    raise RuntimeError(f"Unknown method {method}")
                if validate:
                    validate_fastq_gz(args.out_r1)
                    validate_fastq_gz(args.out_r2, allow_empty=True)
                summary.update({"sample": args.sample, "requested_method": args.method, "attempts": attempts + [{"method": method, "status": "success"}]})
                write_json(args.summary, summary)
                log_line(log_handle, f"success_method={method}")
                return
            except Exception as exc:
                attempts.append({"method": method, "status": "failed", "error": str(exc)})
                log_line(log_handle, f"failed_method={method}")
                log_line(log_handle, f"error={exc}")
                for path in (args.out_r1, args.out_r2):
                    tmp = Path(str(path) + ".download")
                    if tmp.exists():
                        tmp.unlink()
                if args.method != "auto":
                    break
    write_json(args.summary, {"sample": args.sample, "requested_method": args.method, "attempts": attempts, "status": "failed"})
    raise SystemExit(f"All download methods failed for sample {args.sample}; see {args.log}")


if __name__ == "__main__":
    main()
