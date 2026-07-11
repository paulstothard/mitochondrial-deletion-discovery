#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
from pathlib import Path

from common import deep_update, ensure_parent, read_tsv, read_yaml, write_tsv, write_yaml


RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
BIOSAMPLE_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def fetch_runinfo(bioproject: str) -> list[dict[str, str]]:
    query = urllib.parse.urlencode({"acc": bioproject})
    with urllib.request.urlopen(f"{RUNINFO_URL}?{query}", timeout=120) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines()))
    if not rows:
        raise SystemExit(f"NCBI returned no SRA run rows for {bioproject}")
    return rows


def fetch_ena_fastq_map(runs: list[str]) -> dict[str, dict[str, str]]:
    runs = sorted({run for run in runs if run})
    if not runs:
        return {}
    result = {}

    def parse_text(text: str) -> None:
        rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        for row in rows:
            run = row.get("run_accession", "")
            fastq_ftp = row.get("fastq_ftp", "")
            if not run or not fastq_ftp:
                continue
            urls = [item if item.startswith(("http://", "https://")) else f"https://{item}" for item in fastq_ftp.split(";") if item]
            result[run] = {
                "fastq_1": urls[0] if urls else "",
                "fastq_2": urls[1] if len(urls) > 1 else "",
                "fastq_md5": row.get("fastq_md5", ""),
                "fastq_bytes": row.get("fastq_bytes", ""),
            }

    query = urllib.parse.urlencode(
        {
            "accession": ",".join(runs),
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
            "format": "tsv",
        }
    )
    try:
        parse_text(fetch_url(f"{ENA_FILEREPORT_URL}?{query}", attempts=2))
    except Exception:
        pass
    if len(result) == len(runs):
        return result

    for run in runs:
        if run in result:
            continue
        query = urllib.parse.urlencode(
            {
                "accession": run,
                "result": "read_run",
                "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                "format": "tsv",
            }
        )
        try:
            parse_text(fetch_url(f"{ENA_FILEREPORT_URL}?{query}", attempts=3))
        except Exception:
            pass
        time.sleep(0.1)
    return result


def parse_biosample_attributes(biosample) -> dict[str, str]:
    accession = biosample.attrib.get("accession", "")
    attrs = {"BioSample": accession}
    for id_node in biosample.findall(".//Ids/Id"):
        if id_node.attrib.get("db_label") == "Sample name" and id_node.text:
            attrs["Sample Name"] = id_node.text
    for node in biosample.findall(".//Attributes/Attribute"):
        name = node.attrib.get("harmonized_name") or node.attrib.get("attribute_name") or node.attrib.get("display_name")
        if name and node.text:
            attrs[name] = node.text
            attrs[name.lower()] = node.text
    return attrs


def fetch_url(url: str, timeout: int = 120, attempts: int = 5) -> str:
    last_error = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except HTTPError as exc:
            last_error = exc
            if exc.code not in {400, 429, 500, 502, 503, 504}:
                raise
        except URLError as exc:
            last_error = exc
        time.sleep(min(30, 2**attempt))
    raise last_error


def fetch_biosample_attribute_map_batch(accessions: list[str]) -> dict[str, dict[str, str]]:
    accessions = sorted({acc for acc in accessions if acc})
    if not accessions:
        return {}
    query = urllib.parse.urlencode({"db": "biosample", "id": ",".join(accessions), "retmode": "xml"})
    text = fetch_url(f"{BIOSAMPLE_EFETCH_URL}?{query}", attempts=2)
    root = ET.fromstring(text)
    parsed = {}
    for biosample in root.findall(".//BioSample"):
        attrs = parse_biosample_attributes(biosample)
        if attrs.get("BioSample"):
            parsed[attrs["BioSample"]] = attrs
    return parsed


def fetch_biosample_attribute_map(accessions: list[str]) -> dict[str, dict[str, str]]:
    accessions = sorted({acc for acc in accessions if acc})
    if not accessions:
        return {}
    try:
        parsed = fetch_biosample_attribute_map_batch(accessions)
        if parsed:
            return parsed
    except HTTPError as exc:
        if exc.code != 400:
            raise

    parsed = {}
    for accession in accessions:
        query = urllib.parse.urlencode({"db": "biosample", "id": accession, "retmode": "xml"})
        try:
            text = fetch_url(f"{BIOSAMPLE_EFETCH_URL}?{query}", attempts=6)
        except Exception:
            continue
        root = ET.fromstring(text)
        biosample = root.find(".//BioSample")
        if biosample is not None:
            attrs = parse_biosample_attributes(biosample)
            if attrs.get("BioSample"):
                parsed[attrs["BioSample"]] = attrs
        time.sleep(0.35)
    return parsed


def enrich_with_biosamples(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    enriched = []
    accessions = [row.get("BioSample") or row.get("biosample") for row in rows]
    cache = fetch_biosample_attribute_map(accessions)
    for row in rows:
        merged = dict(row)
        accession = row.get("BioSample") or row.get("biosample")
        if accession and accession in cache:
            for key, value in cache[accession].items():
                merged.setdefault(key, value)
                if not merged.get(key):
                    merged[key] = value
        enriched.append(merged)
    return enriched


def load_run_table(path: str) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        delimiter = "\t" if first_line.count("\t") > first_line.count(",") else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def existing_sample_cache_paths(dataset: str) -> list[Path]:
    return [
        Path("metadata/generated") / f"{dataset}.samples.tsv",
        Path("results") / dataset / "analysis" / "deletion_burden.tsv",
        Path("results") / dataset / f"{dataset}_deliverables" / "tables" / "deletion_burden.tsv",
    ]


def existing_run_table_cache_paths(dataset: str) -> list[Path]:
    return [
        Path("metadata/cache") / f"{dataset}.sra_run_table.csv",
        Path("metadata/generated") / f"{dataset}.sra_run_table.csv",
    ]


def cached_samples_from_existing_outputs(dataset_cfg: dict) -> list[dict[str, str]]:
    dataset = dataset_cfg["dataset"]["name"]
    wanted = {
        "sample",
        "dataset",
        "species",
        "run_accession",
        "fastq_1",
        "fastq_2",
        "fastq_md5",
        "fastq_bytes",
        "layout",
        "age",
        "treatment",
        "condition",
        "tissue",
        "bioproject",
        "biosample",
        "sra_study",
        "sample_name",
        "biological_replicate",
    }
    for path in existing_sample_cache_paths(dataset):
        if not path.exists():
            continue
        rows = read_tsv(str(path))
        if rows and {"sample", "dataset", "species"}.issubset(rows[0].keys()):
            return [{key: value for key, value in row.items() if key in wanted} for row in rows]

    qc_root = Path("results") / dataset / "qc"
    rows = []
    if qc_root.exists():
        for read_json in sorted(qc_root.glob("*/read_input.json")):
            try:
                data = json.loads(read_json.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            sample = str(data.get("sample") or read_json.parent.name)
            parts = sample.split("_")
            age = parts[1] if len(parts) > 1 else derive_age({"Run": data.get("run_accession", sample)})
            treatment = parts[2] if len(parts) > 2 else derive_treatment({"Run": data.get("run_accession", sample)})
            rows.append(
                {
                    "sample": sample,
                    "dataset": dataset,
                    "species": dataset_cfg["dataset"]["species"],
                    "run_accession": data.get("run_accession", ""),
                    "fastq_1": "",
                    "fastq_2": "",
                    "layout": "paired" if data.get("paired") else "single",
                    "age": age,
                    "treatment": treatment,
                    "condition": sanitize(f"{age}_{treatment}"),
                    "bioproject": dataset_cfg["dataset"].get("bioproject", ""),
                }
            )
    return rows


def sanitize(value: str) -> str:
    value = str(value or "NA").strip()
    value = re.sub(r"[^A-Za-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "NA"


def derive_age(row: dict[str, str]) -> str:
    for key in ("Age", "age"):
        if row.get(key):
            return sanitize(row[key])
    biological = row.get("Biological_Replicate", "")
    match = re.search(r"(\d+)[A-Za-z]*", biological)
    return f"{match.group(1)}mo" if match else "NA"


def derive_treatment(row: dict[str, str]) -> str:
    for key in ("treatment", "Treatment"):
        if row.get(key):
            return sanitize(row[key])
    biological = row.get("Biological_Replicate", "")
    if re.match(r"\d+G", biological):
        return "GPA"
    if re.match(r"\d+C", biological):
        return "Control"
    return "NA"


def derive_replicate(row: dict[str, str]) -> str:
    biological = row.get("Biological_Replicate") or row.get("Sample Name") or row.get("Run", "")
    match = re.search(r"replicate[_ -]*(\d+)", biological, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"[-_](\d+)$", biological)
    if match:
        return match.group(1)
    return sanitize(row.get("Run", "1"))


def normalize_layout(value: str) -> str:
    value = (value or "").strip().lower()
    if value == "paired":
        return "paired"
    return "single"


def validate_dataset_inputs(dataset_cfg: dict, samples: list[dict[str, str]]) -> None:
    dataset = dataset_cfg.get("dataset", {}) or {}
    dataset_name = str(dataset.get("name", "")).strip()
    species = str(dataset.get("species", "")).strip()
    strategy = str(dataset.get("library_strategy", "unknown")).strip().lower()
    group_columns = dataset.get("group_columns", []) or []

    sample_ids = [str(row.get("sample", "")).strip() for row in samples]
    if any(not sample for sample in sample_ids):
        raise SystemExit("Resolved sample table contains an empty sample identifier")
    duplicates = sorted({sample for sample in sample_ids if sample_ids.count(sample) > 1})
    if duplicates:
        raise SystemExit(f"Resolved sample table contains duplicate sample identifiers: {', '.join(duplicates)}")

    expected_layout = ""
    if strategy == "paired_end_short_read":
        expected_layout = "paired"
    elif strategy == "single_end_short_read":
        expected_layout = "single"

    for row in samples:
        sample = str(row.get("sample", "")).strip()
        row_dataset = str(row.get("dataset", "")).strip()
        row_species = str(row.get("species", "")).strip()
        if dataset_name and row_dataset != dataset_name:
            raise SystemExit(f"Sample {sample} belongs to dataset {row_dataset!r}, expected {dataset_name!r}")
        if species and row_species != species:
            raise SystemExit(f"Sample {sample} has species {row_species!r}, expected {species!r}")
        raw_layout = str(row.get("layout", "")).strip().lower()
        if raw_layout not in {"single", "single-end", "single_end", "se", "paired", "paired-end", "paired_end", "pe"}:
            raise SystemExit(f"Sample {sample} has missing or unsupported layout {raw_layout!r}")
        layout = normalize_layout(raw_layout)
        if expected_layout and layout != expected_layout:
            raise SystemExit(
                f"Sample {sample} has {layout}-end layout but dataset.library_strategy is {strategy}"
            )
        if layout == "paired" and not str(row.get("fastq_2", "")).strip():
            raise SystemExit(f"Paired-end sample {sample} does not define fastq_2")
        if layout == "single" and str(row.get("fastq_2", "")).strip():
            raise SystemExit(f"Single-end sample {sample} unexpectedly defines fastq_2")
        if not str(row.get("fastq_1", "")).strip() and not str(row.get("run_accession", "")).strip():
            raise SystemExit(f"Sample {sample} defines neither fastq_1 nor run_accession")
        empty_groups = [column for column in group_columns if not str(row.get(column, "")).strip()]
        if empty_groups:
            raise SystemExit(f"Sample {sample} has empty grouping values: {', '.join(empty_groups)}")

    reference = (dataset_cfg.get("references", {}) or {}).get(species, {}) or {}
    mt_length = int(reference.get("mt_length", 0) or 0)
    for target in (dataset_cfg.get("analysis", {}) or {}).get("known_deletions", []) or []:
        values = [target.get(key) for key in ("left_breakpoint", "right_breakpoint", "deleted_size")]
        if mt_length <= 0 or any(value in {None, ""} for value in values):
            continue
        left, right, configured_size = map(int, values)
        implied_size = right - left - 1 if right > left else mt_length - left + right - 1
        if configured_size != implied_size:
            name = str(target.get("name", "unnamed configured deletion"))
            raise SystemExit(
                f"Configured deletion {name!r} has deleted_size {configured_size}, but retained "
                f"breakpoints {left}->{right} imply {implied_size} bases on mtDNA length {mt_length}"
            )


def make_sample_id(dataset: str, row: dict[str, str], template: str | None) -> str:
    age = derive_age(row)
    treatment = derive_treatment(row)
    replicate = derive_replicate(row)
    values = {
        "dataset": sanitize(dataset),
        "run": sanitize(row.get("Run", "")),
        "age": age,
        "treatment": treatment,
        "condition": sanitize(f"{age}_{treatment}"),
        "replicate": replicate,
        "sample_name": sanitize(row.get("Sample Name", "")),
    }
    if template:
        try:
            return sanitize(template.format(**values))
        except KeyError as exc:
            raise SystemExit(f"Unknown sample_naming template field: {exc}") from exc
    return sanitize(f"{dataset}_{row.get('Run', '')}")


def resolve_samples(dataset_cfg: dict, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    dataset = dataset_cfg["dataset"]["name"]
    species = dataset_cfg["dataset"]["species"]
    template = dataset_cfg["dataset"].get("sample_naming", {}).get("template")
    ena_fastqs = fetch_ena_fastq_map([row.get("Run") or row.get("run_accession") for row in rows])
    resolved = []
    seen = set()
    for row in rows:
        run = row.get("Run") or row.get("run_accession")
        if not run:
            continue
        age = derive_age(row)
        treatment = derive_treatment(row)
        condition = sanitize(f"{age}_{treatment}")
        sample = make_sample_id(dataset, row, template)
        base = sample
        index = 2
        while sample in seen:
            sample = f"{base}_{index}"
            index += 1
        seen.add(sample)
        resolved.append(
            {
                "sample": sample,
                "dataset": dataset,
                "species": species,
                "run_accession": run,
                "fastq_1": ena_fastqs.get(run, {}).get("fastq_1", ""),
                "fastq_2": ena_fastqs.get(run, {}).get("fastq_2", ""),
                "layout": normalize_layout(row.get("LibraryLayout", row.get("layout", ""))),
                "age": age,
                "treatment": treatment,
                "condition": condition,
                "tissue": row.get("tissue") or row.get("Tissue") or row.get("source_name", ""),
                "bioproject": row.get("BioProject", dataset_cfg["dataset"].get("bioproject", "")),
                "biosample": row.get("BioSample", ""),
                "sra_study": row.get("SRA Study", row.get("SRAStudy", "")),
                "sample_name": row.get("Sample Name", ""),
                "biological_replicate": row.get("Biological_Replicate") or row.get("biological_replicate", ""),
            }
        )
    if not resolved:
        raise SystemExit("No usable samples were resolved from SRA metadata")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--defaults", required=True)
    parser.add_argument("--dataset-config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--bioproject", default="")
    parser.add_argument("--sra-run-table", default="")
    parser.add_argument("--sample-table", default="")
    parser.add_argument("--out-run-table", required=True)
    parser.add_argument("--out-samples", required=True)
    parser.add_argument("--out-config", required=True)
    args = parser.parse_args()

    defaults = read_yaml(args.defaults)
    dataset_cfg = deep_update(defaults, read_yaml(args.dataset_config))
    if args.sample_table:
        samples = read_tsv(args.sample_table)
        rows = samples
        ensure_parent(args.out_run_table)
        fieldnames = list(samples[0].keys()) if samples else []
        with open(args.out_run_table, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(samples)
    else:
        if args.sra_run_table:
            rows = load_run_table(args.sra_run_table)
        elif args.bioproject:
            cached_run_table = next((path for path in existing_run_table_cache_paths(args.dataset) if path.exists()), None)
            if cached_run_table is not None:
                rows = load_run_table(str(cached_run_table))
                print(f"Using cached SRA run table at {cached_run_table}")
            else:
                try:
                    rows = enrich_with_biosamples(fetch_runinfo(args.bioproject))
                except (HTTPError, URLError) as exc:
                    samples = cached_samples_from_existing_outputs(dataset_cfg)
                    if not samples:
                        raise
                    rows = samples
                    print(f"Warning: using cached resolved sample metadata after metadata fetch failed: {exc}")
        else:
            raise SystemExit("Provide dataset.bioproject, dataset.sra_run_table, or dataset.sample_table")
        ensure_parent(args.out_run_table)
        with open(args.out_run_table, "w", encoding="utf-8", newline="") as handle:
            fieldnames = list(rows[0].keys())
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        if {"sample", "dataset", "species"}.issubset(rows[0].keys()):
            samples = rows
        else:
            samples = resolve_samples(dataset_cfg, rows)
        if args.bioproject:
            cache_path = Path("metadata/cache") / f"{args.dataset}.sra_run_table.csv"
            ensure_parent(str(cache_path))
            with cache_path.open("w", encoding="utf-8", newline="") as handle:
                fieldnames = list(rows[0].keys())
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

    required = ["sample", "dataset", "species"]
    missing = [col for col in required if col not in samples[0]]
    if missing:
        raise SystemExit(f"Resolved sample table missing columns: {', '.join(missing)}")
    group_columns = dataset_cfg["dataset"].get("group_columns", [])
    missing_groups = [col for col in group_columns if col not in samples[0]]
    if missing_groups:
        raise SystemExit(f"Resolved sample table missing group columns: {', '.join(missing_groups)}")
    validate_dataset_inputs(dataset_cfg, samples)

    fieldnames = [
        "sample",
        "dataset",
        "species",
        "run_accession",
        "fastq_1",
        "fastq_2",
        "fastq_md5",
        "fastq_bytes",
        "layout",
        "age",
        "treatment",
        "condition",
        "tissue",
        "bioproject",
        "biosample",
        "sra_study",
        "sample_name",
        "biological_replicate",
    ]
    fieldnames = [name for name in fieldnames if any(name in row for row in samples)]
    for row in samples:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    write_tsv(args.out_samples, samples, fieldnames=fieldnames)
    write_yaml(args.out_config, dataset_cfg)


if __name__ == "__main__":
    main()
