#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml

from common import deep_update, ensure_parent, write_tsv


HIT_FIELDS = [
    "sample",
    "deletion_id",
    "deletion_name",
    "search_strategy",
    "read_id",
    "mate",
    "matched_sequence_ids",
    "matched_orientation",
]

SUMMARY_FIELDS = [
    "sample",
    "deletion_id",
    "deletion_name",
    "search_strategy",
    "reads_examined",
    "matching_reads",
    "matching_reads_per_million_examined",
]


def load_config(defaults: str, dataset_config: str) -> dict:
    with open(defaults, "r", encoding="utf-8") as handle:
        base = yaml.safe_load(handle) or {}
    with open(dataset_config, "r", encoding="utf-8") as handle:
        override = yaml.safe_load(handle) or {}
    return deep_update(base, override)


def read_samples(path: str) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_from_fastq(path: str, suffix: str) -> str:
    name = Path(path).name
    if name.endswith(suffix):
        return name[: -len(suffix)]
    if name.endswith(".mt_evidence.fastq.gz"):
        return name[: -len(".mt_evidence.fastq.gz")]
    return name.split("_R", 1)[0]


def fastq_records(path: str):
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            yield header[1:].strip().split()[0], seq


def read_fragment_counts(paths: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for path in paths:
        if not path or not Path(path).exists():
            continue
        with open(path, "r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row.get("sample", "")
                try:
                    counts[sample] = int(float(row.get("fragments", "0")))
                except ValueError:
                    counts[sample] = 0
    return counts


def read_mt_evidence_counts(paths: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for path in paths:
        if not path or not Path(path).exists():
            continue
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
        except (OSError, json.JSONDecodeError):
            continue
        sample = Path(path).name.replace(".mt_read_summary.json", "")
        value = data.get("mt_evidence_fastq_records_written", data.get("mt_evidence_reads_selected", 0))
        try:
            counts[sample] = int(float(value))
        except (TypeError, ValueError):
            counts[sample] = 0
    return counts


def sequence_variants(item: dict) -> list[tuple[str, str, str]]:
    seq_id = str(item.get("id", "sequence"))
    variants = []
    seq = str(item.get("sequence", "")).strip().upper()
    rc = str(item.get("reverse_complement", "")).strip().upper()
    if seq:
        variants.append((seq_id, "forward", seq))
    if rc and rc != seq:
        variants.append((seq_id, "reverse_complement", rc))
    return variants


def match_single(read_text: str, deletion: dict) -> tuple[bool, list[str], list[str]]:
    matched_ids = []
    orientations = []
    for search_seq in deletion.get("search_sequences", []) or []:
        for seq_id, orientation, seq in sequence_variants(search_seq):
            if seq in read_text:
                matched_ids.append(seq_id)
                orientations.append(orientation)
                return True, matched_ids, orientations
    return False, matched_ids, orientations


def match_multi_required(read_text: str, deletion: dict) -> tuple[bool, list[str], list[str]]:
    matched_ids = []
    orientations = []
    for search_seq in deletion.get("search_sequences", []) or []:
        found = False
        for seq_id, orientation, seq in sequence_variants(search_seq):
            if seq in read_text:
                matched_ids.append(seq_id)
                orientations.append(orientation)
                found = True
                break
        if not found:
            return False, matched_ids, orientations
    return True, matched_ids, orientations


def run_count_command(commands: list[list[str]]) -> int:
    previous = None
    processes = []
    try:
        for command in commands:
            proc = subprocess.Popen(
                command,
                stdin=previous.stdout if previous is not None else None,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=False,
            )
            if previous is not None and previous.stdout is not None:
                previous.stdout.close()
            previous = proc
            processes.append(proc)
        assert previous is not None and previous.stdout is not None
        output = previous.stdout.read().decode("utf-8", errors="ignore").strip()
        return_codes = [proc.wait() for proc in processes]
        if return_codes[-1] not in {0, 1}:
            return 0
        return int(output or "0")
    except (OSError, ValueError):
        for proc in processes:
            if proc.poll() is None:
                proc.kill()
        return 0


def rg_count_fastq_records(path: str, pattern_groups: list[list[str]]) -> int:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return 0
    if not pattern_groups:
        return 0
    if not shutil.which("rg") or not shutil.which("gzip") or not shutil.which("paste") or not shutil.which("wc"):
        return -1
    commands: list[list[str]] = [["gzip", "-cd", path], ["paste", "-", "-", "-", "-"]]
    for patterns in pattern_groups:
        rg_command = ["rg", "-F"]
        for pattern in patterns:
            rg_command.extend(["-e", pattern])
        commands.append(rg_command)
    commands.append(["wc", "-l"])
    return run_count_command(commands)


def seqkit_count_fastq_records(path: str, pattern_groups: list[list[str]]) -> int:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return 0
    if not pattern_groups:
        return 0
    if not shutil.which("seqkit") or not shutil.which("wc"):
        return -1
    commands: list[list[str]] = []
    for i, patterns in enumerate(pattern_groups):
        command = ["seqkit", "grep", "-s", "-j", "4"]
        for pattern in patterns:
            command.extend(["-p", pattern])
        if i == 0:
            command.append(path)
        commands.append(command)
    commands.extend([["seqkit", "seq", "-n"], ["wc", "-l"]])
    return run_count_command(commands)


def seqkit_locate_single_counts(path: str, deletions: list[dict]) -> dict[str, int] | None:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return {}
    if not shutil.which("seqkit"):
        return None
    records = []
    for deletion in deletions:
        strategy = (deletion.get("search_strategy", {}) or {}).get("type", "single_sequence")
        if strategy != "single_sequence":
            continue
        deletion_id = str(deletion.get("id", ""))
        for search_seq in deletion.get("search_sequences", []) or []:
            for seq_id, orientation, seq in sequence_variants(search_seq):
                records.append((f"{deletion_id}|{seq_id}|{orientation}", seq))
    if not records:
        return {}
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".fa", delete=False) as handle:
        pattern_path = Path(handle.name)
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")
    counts: dict[str, set[str]] = {}
    try:
        proc = subprocess.run(
            ["seqkit", "locate", "-j", "4", "-f", str(pattern_path), path],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        if proc.returncode not in {0, 1}:
            return None
        for line in proc.stdout.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read_id = parts[0]
            deletion_id = parts[1].split("|", 1)[0]
            counts.setdefault(deletion_id, set()).add(read_id)
        return {deletion_id: len(read_ids) for deletion_id, read_ids in counts.items()}
    finally:
        pattern_path.unlink(missing_ok=True)


def configured_pattern_groups(deletion: dict) -> list[list[str]]:
    strategy = (deletion.get("search_strategy", {}) or {}).get("type", "single_sequence")
    groups = []
    if strategy == "multi_sequence_required":
        for search_seq in deletion.get("search_sequences", []) or []:
            variants = [seq for _, _, seq in sequence_variants(search_seq)]
            if variants:
                groups.append(variants)
    else:
        variants = []
        for search_seq in deletion.get("search_sequences", []) or []:
            variants.extend([seq for _, _, seq in sequence_variants(search_seq)])
        if variants:
            groups.append(variants)
    return groups


def read_pair_iterator(r1: str, r2: str):
    r1_iter = fastq_records(r1)
    r2_iter = fastq_records(r2) if r2 and Path(r2).exists() and Path(r2).stat().st_size > 0 else None
    if r1_iter is None:
        return
    if r2_iter is None:
        for read_id, seq1 in r1_iter:
            yield read_id, seq1, "", "R1"
    else:
        for (read_id, seq1), (_, seq2) in zip(r1_iter, r2_iter):
            yield read_id, seq1, seq2, "R1+R2"


def external_count_all_searches(path: str, searches: list[dict]) -> dict[str, int] | None:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return {}
    decompressor = shutil.which("pigz") or shutil.which("gzip")
    perl = shutil.which("perl")
    counter = Path(__file__).resolve().parent / "count_known_sequences.pl"
    if not decompressor or not perl or not counter.exists():
        return None
    rows = []
    for deletion in searches:
        deletion_id = str(deletion.get("id", ""))
        strategy = (deletion.get("search_strategy", {}) or {}).get("type", "single_sequence")
        if strategy == "multi_sequence_required":
            for group_index, search_seq in enumerate(deletion.get("search_sequences", []) or []):
                for _, _, seq in sequence_variants(search_seq):
                    rows.append((deletion_id, strategy, group_index, seq))
        else:
            for search_seq in deletion.get("search_sequences", []) or []:
                for _, _, seq in sequence_variants(search_seq):
                    rows.append((deletion_id, strategy, 0, seq))
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".tsv", delete=False) as handle:
        pattern_path = Path(handle.name)
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    try:
        decompress_command = [decompressor, "-dc", path] if Path(decompressor).name == "pigz" else [decompressor, "-cd", path]
        decompress = subprocess.Popen(decompress_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        count_proc = subprocess.run(
            [perl, str(counter), str(pattern_path)],
            stdin=decompress.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        )
        if decompress.stdout is not None:
            decompress.stdout.close()
        decompress.wait()
        if count_proc.returncode != 0:
            return None
        counts = {str(deletion.get("id", "")): 0 for deletion in searches}
        for line in count_proc.stdout.splitlines():
            deletion_id, count = line.split("\t", 1)
            counts[deletion_id] = int(count)
        return counts
    finally:
        pattern_path.unlink(missing_ok=True)


def compiled_searches(searches: list[dict]) -> list[dict]:
    compiled = []
    for deletion in searches:
        strategy = (deletion.get("search_strategy", {}) or {}).get("type", "single_sequence")
        groups = []
        if strategy == "multi_sequence_required":
            for search_seq in deletion.get("search_sequences", []) or []:
                variants = [
                    {"sequence_id": seq_id, "orientation": orientation, "pattern": seq.encode("ascii")}
                    for seq_id, orientation, seq in sequence_variants(search_seq)
                ]
                if variants:
                    groups.append(variants)
        else:
            variants = []
            for search_seq in deletion.get("search_sequences", []) or []:
                variants.extend(
                    {"sequence_id": seq_id, "orientation": orientation, "pattern": seq.encode("ascii")}
                    for seq_id, orientation, seq in sequence_variants(search_seq)
                )
            if variants:
                groups.append(variants)
        compiled.append(
            {
                "id": str(deletion.get("id", "")),
                "name": deletion.get("name", ""),
                "strategy": strategy,
                "groups": groups,
            }
        )
    return compiled


def match_compiled_search(seq: bytes, item: dict) -> tuple[bool, list[str], list[str]]:
    matched_ids = []
    orientations = []
    groups = item["groups"]
    if not groups:
        return False, matched_ids, orientations
    for group in groups:
        group_match = None
        for variant in group:
            if variant["pattern"] in seq:
                group_match = variant
                break
        if group_match is None:
            return False, matched_ids, orientations
        matched_ids.append(str(group_match["sequence_id"]))
        orientations.append(str(group_match["orientation"]))
    return True, matched_ids, orientations


def scan_fastq_for_searches(path: str, searches: list[dict], sample: str = "", mate: str = "R1") -> tuple[dict[str, int], int, list[dict[str, object]]]:
    counts = {item["id"]: 0 for item in searches}
    hit_rows: list[dict[str, object]] = []
    examined = 0
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return counts, examined, hit_rows
    with gzip.open(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            examined += 1
            read_id = header[1:].decode("utf-8", errors="ignore").strip()
            for item in searches:
                matched, matched_ids, orientations = match_compiled_search(seq, item)
                if matched:
                    counts[item["id"]] += 1
                    hit_rows.append(
                        {
                            "sample": sample,
                            "deletion_id": item["id"],
                            "deletion_name": item["name"],
                            "search_strategy": item["strategy"],
                            "read_id": read_id,
                            "mate": mate,
                            "matched_sequence_ids": ",".join(matched_ids),
                            "matched_orientation": ",".join(orientations),
                        }
                    )
    return counts, examined, hit_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--defaults", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--r1-files", nargs="*", default=[])
    parser.add_argument("--r2-files", nargs="*", default=[])
    parser.add_argument("--fragment-counts", nargs="*", default=[])
    parser.add_argument("--mt-summaries", nargs="*", default=[])
    parser.add_argument("--hits", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    config = load_config(args.defaults, args.config)
    searches = config.get("analysis", {}).get("known_sequence_searches", []) or []
    compiled = compiled_searches(searches)
    samples = read_samples(args.samples)
    r1_by_sample = {sample_from_fastq(path, "_R1.fastq.gz"): path for path in args.r1_files}
    r2_by_sample = {sample_from_fastq(path, "_R2.fastq.gz"): path for path in args.r2_files}
    fragment_counts = read_fragment_counts(args.fragment_counts)
    mt_evidence_counts = read_mt_evidence_counts(args.mt_summaries)

    hit_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    if searches:
        for sample_row in samples:
            sample = sample_row["sample"]
            r1 = r1_by_sample.get(sample, "")
            r2 = r2_by_sample.get(sample, "")
            counts = {item.get("id", ""): 0 for item in searches}
            examined = mt_evidence_counts.get(sample, fragment_counts.get(sample, 0))
            r1_counts, examined_from_scan, r1_hits = scan_fastq_for_searches(r1, compiled, sample=sample, mate="R1")
            if examined_from_scan and not examined:
                examined = examined_from_scan
            for deletion_id, count in r1_counts.items():
                counts[deletion_id] = counts.get(deletion_id, 0) + count
            hit_rows.extend(r1_hits)
            if r2:
                r2_counts, _, r2_hits = scan_fastq_for_searches(r2, compiled, sample=sample, mate="R2")
                for deletion_id, count in r2_counts.items():
                    counts[deletion_id] = counts.get(deletion_id, 0) + count
                hit_rows.extend(r2_hits)
            for deletion in searches:
                deletion_id = str(deletion.get("id", ""))
                matching = counts.get(deletion_id, 0)
                strategy = (deletion.get("search_strategy", {}) or {}).get("type", "single_sequence")
                summary_rows.append(
                    {
                        "sample": sample,
                        "deletion_id": deletion_id,
                        "deletion_name": deletion.get("name", ""),
                        "search_strategy": strategy,
                        "reads_examined": examined,
                        "matching_reads": matching,
                        "matching_reads_per_million_examined": (matching / examined * 1_000_000) if examined else 0.0,
                    }
                )

    ensure_parent(args.hits)
    ensure_parent(args.summary)
    write_tsv(args.hits, hit_rows, HIT_FIELDS)
    write_tsv(args.summary, summary_rows, SUMMARY_FIELDS)


if __name__ == "__main__":
    main()
