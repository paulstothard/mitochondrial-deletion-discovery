#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pandas as pd

from circular_deletions import (
    affected_feature_impact,
    features_at,
    known_deletion_match,
    nearest_feature,
    replication_arc_annotation,
    size_class,
)
from common import ensure_parent, read_yaml


FEATURE_PRIORITY = {
    "gene": 0,
    "tRNA": 1,
    "rRNA": 1,
    "CDS": 2,
    "transcript": 3,
    "exon": 4,
    "region": 5,
}


def inferred_feature_type(name: object, feature_type: object) -> str:
    text = str(name or "").lower()
    raw_type = str(feature_type or "").lower()
    if raw_type == "region":
        return "region"
    if "trna" in raw_type or text.startswith(("mt-t", "trn")):
        return "tRNA"
    if "rrna" in raw_type or text.startswith(("mt-r", "rrn")):
        return "rRNA"
    if "protein_coding" in raw_type or "cds" in raw_type or text.startswith(("mt-co", "mt-cy", "mt-nd", "mt-atp")):
        return "protein_coding"
    return str(feature_type or "")


def append_configured_regions(features: pd.DataFrame, config: dict, mt_length: int) -> pd.DataFrame:
    rows = []
    for item in config.get("analysis", {}).get("mt_regions", []) or []:
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        start = int(item.get("start", 0) or 0)
        end = int(item.get("end", 0) or 0)
        if start <= 0 or end <= 0:
            continue
        pieces = [(start, end)] if start <= end else [(start, mt_length), (1, end)]
        for piece_start, piece_end in pieces:
            rows.append(
                {
                    "contig": "",
                    "start": piece_start,
                    "end": piece_end,
                    "strand": ".",
                    "feature_type": "region",
                    "gene_id": name,
                    "gene_name": name,
                    "transcript_id": "",
                    "product": item.get("reason", ""),
                }
            )
    if not rows:
        return features
    return pd.concat([features, pd.DataFrame(rows)], ignore_index=True)


def biological_features(features: pd.DataFrame) -> pd.DataFrame:
    if features.empty:
        return features.copy()
    out = features.copy()
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out = out.dropna(subset=["start", "end"])
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out["_name"] = out.apply(lambda row: str(row.get("gene_name") or row.get("gene_id") or row.get("product") or "").strip(), axis=1)
    out = out[out["_name"] != ""].copy()
    out["_priority"] = out["feature_type"].map(FEATURE_PRIORITY).fillna(99).astype(int)
    out["_length"] = out["end"] - out["start"] + 1
    out["_feature_key"] = out.apply(
        lambda row: f"{row['_name']}:{row['start']}-{row['end']}" if row.get("feature_type") == "region" else row["_name"],
        axis=1,
    )
    out = out.sort_values(["_feature_key", "_priority", "_length", "start", "end"])
    out = out.drop_duplicates(subset=["_feature_key"], keep="first").copy()
    out["gene_name"] = out["_name"]
    return out.drop(columns=["_name", "_feature_key", "_priority", "_length"], errors="ignore").sort_values(["start", "end", "gene_name"]).reset_index(drop=True)


def apply_feature_aliases(features: pd.DataFrame, config: dict) -> pd.DataFrame:
    if features.empty:
        return features.copy()
    out = features.copy()
    out["raw_gene_name"] = out.get("gene_name", "").astype(str)
    out["display_name"] = out.get("gene_name", "").astype(str)
    for item in config.get("annotations", {}).get("feature_aliases", []) or []:
        display = str(item.get("display_name", "")).strip()
        if not display:
            continue
        mask = pd.Series(True, index=out.index)
        if item.get("raw_name"):
            mask &= out.get("gene_name", "").astype(str).eq(str(item["raw_name"]))
        if item.get("gene_id"):
            mask &= out.get("gene_id", "").astype(str).eq(str(item["gene_id"]))
        if item.get("start") is not None:
            mask &= pd.to_numeric(out.get("start"), errors="coerce").eq(int(item["start"]))
        if item.get("end") is not None:
            mask &= pd.to_numeric(out.get("end"), errors="coerce").eq(int(item["end"]))
        out.loc[mask, "display_name"] = display
    out["gene_name"] = out["display_name"]
    out["feature_type"] = out.apply(lambda row: inferred_feature_type(row.get("gene_name"), row.get("feature_type")), axis=1)
    return out


def feature_label(row: pd.Series) -> str:
    name = row.get("gene_name") or row.get("gene_id") or row.get("product") or "unknown"
    return f"{name}:{row.get('feature_type', '')}:{int(row['start'])}-{int(row['end'])}"


def simple_feature_name(annotation: str) -> str:
    if not annotation:
        return ""
    return str(annotation).split(";")[0].split(":")[0].lower()


def expected_transcript_annotation(left_feature: str, right_feature: str, expected_pairs: list[dict]) -> tuple[str, str]:
    left = simple_feature_name(left_feature)
    right = simple_feature_name(right_feature)
    for pair in expected_pairs:
        expected_left = str(pair.get("left_feature", "")).lower()
        expected_right = str(pair.get("right_feature", "")).lower()
        bidirectional = bool(pair.get("bidirectional", False))
        if left == expected_left and right == expected_right:
            return "expected_transcript_junction", pair.get("reason", "configured_expected_transcript")
        if bidirectional and left == expected_right and right == expected_left:
            return "expected_transcript_junction", pair.get("reason", "configured_expected_transcript")
    return "candidate_deletion_junction", "not_configured_expected_transcript_pair"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--config", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters, sep="\t")
    features = pd.read_csv(args.features, sep="\t")
    config = read_yaml(args.config) if args.config else {}
    features = apply_feature_aliases(biological_features(append_configured_regions(features, config, args.mt_length)), config)
    expected_pairs = config.get("annotations", {}).get("expected_adjacent_transcripts", []) or []
    if clusters.empty:
        ensure_parent(args.output)
        clusters.to_csv(args.output, sep="\t", index=False)
        return
    annotations = []
    for _, row in clusters.iterrows():
        left = int(row["left_breakpoint"])
        right = int(row["right_breakpoint"])
        left_overlap = features_at(features, left)
        right_overlap = features_at(features, right)
        impact = affected_feature_impact(features, left, right, args.mt_length)
        arc_annotation = replication_arc_annotation(config, left, right, args.mt_length)
        known_label, known_reason = known_deletion_match(left, right, int(row["deleted_size"]), config, args.mt_length)
        annotations.append(
            {
                "left_feature_overlap": left_overlap,
                "right_feature_overlap": right_overlap,
                "nearest_left_feature": nearest_feature(features, left),
                "nearest_right_feature": nearest_feature(features, right),
                "deleted_interval_features": impact.affected_features,
                "affected_feature_label": impact.affected_feature_label,
                "affected_features": impact.affected_features,
                "fully_removed_features": impact.fully_removed_features,
                "partially_overlapped_features": impact.partially_overlapped_features,
                "feature_impact_class": impact.feature_impact_class,
                "per_feature_overlap_details": impact.per_feature_hits,
                "size_class": size_class(int(row["deleted_size"])),
                **arc_annotation,
                "known_deletion_label": known_label,
                "known_deletion_match_reason": known_reason,
            }
        )
    annotation_df = pd.DataFrame(annotations)
    statuses = [
        expected_transcript_annotation(row["nearest_left_feature"], row["nearest_right_feature"], expected_pairs)
        for _, row in annotation_df.iterrows()
    ]
    annotation_df["junction_interpretation"] = [status for status, _ in statuses]
    annotation_df["junction_interpretation_reason"] = [reason for _, reason in statuses]
    annotated = pd.concat([clusters.reset_index(drop=True), annotation_df], axis=1)
    ensure_parent(args.output)
    annotated.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
