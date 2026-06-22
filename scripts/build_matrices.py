#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import ensure_parent, read_json, read_tsv, read_yaml


def load_fragment_counts(paths: list[str]) -> dict[str, float]:
    values = {}
    for path in paths:
        rows = read_tsv(path)
        for row in rows:
            values[row["sample"]] = float(row.get("million_fragments") or 0)
    return values


def load_mt_counts(paths: list[str]) -> dict[str, float]:
    values = {}
    for path in paths:
        data = read_json(path)
        sample = Path(path).name.split(".")[0]
        mt_reads = data.get("mt_evidence_fastq_records_written", data.get("high_confidence_mt", 0))
        values[sample] = float(mt_reads or 0) / 1_000_000
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--id-map", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--ffpm", required=True)
    parser.add_argument("--mtpm", required=True)
    parser.add_argument("--presence", required=True)
    parser.add_argument("--config", default="")
    parser.add_argument("--include-expected-transcript-junctions", action="store_true")
    parser.add_argument("--fragment-counts", nargs="*", default=[])
    parser.add_argument("--mt-summaries", nargs="*", default=[])
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t")
    clusters = pd.read_csv(args.clusters, sep="\t")
    config = read_yaml(args.config) if args.config else {}
    exclude_expected = bool(config.get("junctions", {}).get("exclude_expected_transcript_junctions", True))
    if args.include_expected_transcript_junctions:
        exclude_expected = False
    if exclude_expected and "junction_interpretation" in clusters.columns:
        clusters = clusters[clusters["junction_interpretation"] != "expected_transcript_junction"].copy()
    id_map = pd.read_csv(args.id_map, sep="\t")
    junction_ids = list(clusters["junction_id"]) if not clusters.empty else []
    raw_counts = pd.DataFrame(0, index=samples["sample"], columns=junction_ids, dtype=int)
    if not id_map.empty and junction_ids:
        counts = (
            id_map[id_map["junction_id"].isin(junction_ids)]
            .groupby(["sample", "junction_id"])
            .size()
            .unstack(fill_value=0)
        )
        raw_counts.update(counts)
    raw_counts.index.name = "sample"
    raw = raw_counts.reset_index()
    raw = samples.merge(raw, on="sample", how="left")
    for path in (args.raw, args.ffpm, args.mtpm, args.presence):
        ensure_parent(path)
    raw.to_csv(args.raw, sep="\t", index=False)

    fragment_millions = load_fragment_counts(args.fragment_counts)
    mt_millions = load_mt_counts(args.mt_summaries)
    metadata = raw.drop(columns=junction_ids)
    if junction_ids:
        raw_junctions = raw[junction_ids].astype(float)
        frag_denominator = raw["sample"].map(fragment_millions).replace(0, np.nan).astype(float)
        mt_denominator = raw["sample"].map(mt_millions).replace(0, np.nan).astype(float)
        ffpm = pd.concat([metadata, raw_junctions.div(frag_denominator, axis=0)], axis=1)
        mtpm = pd.concat([metadata, raw_junctions.div(mt_denominator, axis=0)], axis=1)
        presence = pd.concat([metadata, (raw_junctions > 0).astype(int)], axis=1)
    else:
        ffpm = raw.copy()
        mtpm = raw.copy()
        presence = raw.copy()
    ffpm.to_csv(args.ffpm, sep="\t", index=False)
    mtpm.to_csv(args.mtpm, sep="\t", index=False)
    presence.to_csv(args.presence, sep="\t", index=False)


if __name__ == "__main__":
    main()
