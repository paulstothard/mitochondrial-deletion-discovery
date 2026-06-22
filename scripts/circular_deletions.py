#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable

import pandas as pd


def normalize_pos(pos: int, mt_length: int, rotation_start: int = 1) -> int:
    """Convert a 1-based position on a rotated mtDNA reference to original coordinates."""
    if mt_length <= 0:
        raise ValueError("mt_length must be positive")
    return ((int(pos) + int(rotation_start) - 2) % int(mt_length)) + 1


def circular_distance(start: int, end: int, mt_length: int) -> int:
    """Number of bases in the circular interval from start-exclusive to end-exclusive."""
    start = int(start)
    end = int(end)
    if end > start:
        return end - start - 1
    return mt_length - start + end - 1


def interval_pieces(left: int, right: int, mt_length: int) -> list[tuple[int, int]]:
    """Return 1-based closed intervals deleted between breakpoints left -> right."""
    left = int(left)
    right = int(right)
    if right > left:
        return [(left + 1, right - 1)] if right - left > 1 else []
    pieces = []
    if left < mt_length:
        pieces.append((left + 1, mt_length))
    if right > 1:
        pieces.append((1, right - 1))
    return pieces


def interval_length(pieces: Iterable[tuple[int, int]]) -> int:
    return sum(max(0, int(end) - int(start) + 1) for start, end in pieces)


def canonical_breakpoints(left: int, right: int, mt_length: int) -> dict:
    """Choose the shorter circular deleted interval for a breakpoint pair."""
    left = int(left)
    right = int(right)
    forward_size = circular_distance(left, right, mt_length)
    reverse_size = circular_distance(right, left, mt_length)
    if reverse_size < forward_size or (reverse_size == forward_size and (right, left) < (left, right)):
        can_left, can_right, size = right, left, reverse_size
        orientation = "reversed_to_shorter_interval"
    else:
        can_left, can_right, size = left, right, forward_size
        orientation = "reported"
    pieces = interval_pieces(can_left, can_right, mt_length)
    return {
        "canonical_left_breakpoint": can_left,
        "canonical_right_breakpoint": can_right,
        "deleted_size": size,
        "wraps_origin": "yes" if can_right <= can_left else "no",
        "canonical_orientation": orientation,
        "deleted_interval": ";".join(f"{start}-{end}" for start, end in pieces),
    }


def deletion_id(left: int, right: int, deleted_size: int) -> str:
    return f"mtDel_{int(left):05d}_{int(right):05d}_{int(deleted_size):05d}"


def feature_name(row: pd.Series) -> str:
    for col in ("gene_name", "gene_id", "product", "feature_id"):
        value = row.get(col, "")
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return "unknown"


def feature_type(row: pd.Series) -> str:
    value = row.get("feature_type", "")
    return "" if pd.isna(value) else str(value)


def ordered_features(features: pd.DataFrame) -> pd.DataFrame:
    if features.empty:
        return features.copy()
    out = features.copy()
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")
    return out.sort_values(["start", "end", "gene_name" if "gene_name" in out.columns else "feature_type"]).reset_index(drop=True)


def overlaps_interval(row: pd.Series, start: int, end: int) -> bool:
    return int(row["end"]) >= int(start) and int(row["start"]) <= int(end)


def overlap_len(row: pd.Series, start: int, end: int) -> int:
    return max(0, min(int(row["end"]), int(end)) - max(int(row["start"]), int(start)) + 1)


@dataclass(frozen=True)
class FeatureImpact:
    affected_feature_label: str
    affected_features: str
    fully_removed_features: str
    partially_overlapped_features: str
    feature_impact_class: str
    per_feature_hits: str


def affected_feature_impact(features: pd.DataFrame, left: int, right: int, mt_length: int) -> FeatureImpact:
    pieces = interval_pieces(left, right, mt_length)
    ordered = ordered_features(features)
    affected = []
    full = []
    partial = []
    hit_rows = []
    for _, feat in ordered.iterrows():
        if pd.isna(feat.get("start")) or pd.isna(feat.get("end")):
            continue
        total_overlap = sum(overlap_len(feat, start, end) for start, end in pieces)
        if total_overlap <= 0:
            continue
        name = feature_name(feat)
        if name not in affected:
            affected.append(name)
        feat_len = int(feat["end"]) - int(feat["start"]) + 1
        if total_overlap >= feat_len:
            full.append(name)
            overlap_type = "full"
        else:
            partial.append(name)
            overlap_type = "partial"
        hit_rows.append(f"{name}:{overlap_type}:{total_overlap}")
    if not affected:
        label = "intergenic"
        impact_class = "intergenic"
    else:
        label = "+".join(affected)
        if len(affected) == 1:
            impact_class = "single_feature"
        elif len(affected) == 2:
            impact_class = "two_features"
        else:
            coding_types = {
                feature_type(row).lower()
                for _, row in ordered.iterrows()
                if feature_name(row) in affected
            }
            if coding_types and coding_types <= {"cds", "gene", "protein_coding"}:
                impact_class = "multi_feature_protein_coding_only"
            else:
                impact_class = "multi_feature_mixed"
    lower_names = [name.lower() for name in affected]
    if any("d-loop" in name or "d_loop" in name or "control" in name for name in lower_names):
        impact_class = "d_loop_involved"
    elif any(
        "trn" in name
        or "rrn" in name
        or "trna" in name
        or "rrna" in name
        or name.startswith("mt-t")
        or name.startswith("mt-r")
        for name in lower_names
    ):
        impact_class = "rrna_trna_involved"
    return FeatureImpact(
        affected_feature_label=label,
        affected_features=";".join(affected),
        fully_removed_features=";".join(full),
        partially_overlapped_features=";".join(partial),
        feature_impact_class=impact_class,
        per_feature_hits=";".join(hit_rows),
    )


def nearest_feature(features: pd.DataFrame, pos: int) -> str:
    ordered = ordered_features(features)
    if ordered.empty:
        return ""
    distances = ordered.apply(
        lambda row: 0
        if int(row["start"]) <= int(pos) <= int(row["end"])
        else min(abs(int(pos) - int(row["start"])), abs(int(pos) - int(row["end"]))),
        axis=1,
    )
    row = ordered.loc[distances.idxmin()]
    return f"{feature_name(row)}:{feature_type(row)}:{int(row['start'])}-{int(row['end'])};distance={int(distances.min())}"


def features_at(features: pd.DataFrame, pos: int) -> str:
    ordered = ordered_features(features)
    hits = ordered[(ordered["start"] <= int(pos)) & (ordered["end"] >= int(pos))]
    return ";".join(f"{feature_name(row)}:{feature_type(row)}:{int(row['start'])}-{int(row['end'])}" for _, row in hits.iterrows())


def size_class(size: int) -> str:
    size = int(size)
    if size < 250:
        return "small_lt_250bp"
    if size < 1000:
        return "small_250bp_to_999bp"
    if size < 5000:
        return "medium_1kb_to_4999bp"
    if size < 10000:
        return "large_5kb_to_9999bp"
    return "very_large_ge_10kb"


def pos_within_circular_window(pos: int, center: int, tolerance: int, mt_length: int) -> bool:
    pos = int(pos)
    center = int(center)
    tolerance = int(tolerance)
    delta = abs(pos - center)
    return min(delta, mt_length - delta) <= tolerance


def configured_deletion_targets(config: dict, mt_length: int) -> list[dict[str, object]]:
    targets: list[dict[str, object]] = []
    analysis = config.get("analysis", {}) or {}
    for item in analysis.get("known_deletions", []) or []:
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        target: dict[str, object] = {
            "name": name,
            "source": "analysis.known_deletions",
            "match_reason": "configured_known_deletion_match",
            "breakpoint_tolerance_bp": int(item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", 50))),
            "size_tolerance_bp": int(item.get("size_tolerance_bp", item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", 50)))),
        }
        for key in ["left_breakpoint", "right_breakpoint", "deleted_size"]:
            if item.get(key) not in {None, ""}:
                target[key] = int(item[key])
        targets.append(target)

    for item in analysis.get("known_sequence_searches", []) or []:
        name = str(item.get("name") or item.get("id") or "").strip()
        if not name:
            continue
        try:
            if item.get("left_breakpoint") not in {None, ""} and item.get("right_breakpoint") not in {None, ""}:
                left = int(item["left_breakpoint"])
                right = int(item["right_breakpoint"])
            else:
                text = " ".join([str(item.get("id", "")), str(item.get("name", "")), str(item.get("description", ""))])
                match = re.search(r"(\d{2,6})\D+(\d{2,6})", text)
                if not match:
                    continue
                left = int(match.group(1))
                right = int(match.group(2))
            size = int(item.get("deleted_size", circular_distance(left, right, mt_length)))
        except (TypeError, ValueError):
            continue
        breakpoint_tolerance = int(item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", analysis.get("sequence_search_breakpoint_tolerance_bp", 100))))
        size_tolerance = int(item.get("size_tolerance_bp", analysis.get("sequence_search_size_tolerance_bp", 150)))
        duplicate = False
        for target in targets:
            target_left = target.get("left_breakpoint")
            target_right = target.get("right_breakpoint")
            target_size = target.get("deleted_size")
            if target_left is None or target_right is None:
                continue
            target_tol = int(target.get("breakpoint_tolerance_bp", breakpoint_tolerance))
            target_size_tol = int(target.get("size_tolerance_bp", size_tolerance))
            size_ok = target_size is None or abs(size - int(target_size)) <= target_size_tol
            if (
                pos_within_circular_window(left, int(target_left), target_tol, mt_length)
                and pos_within_circular_window(right, int(target_right), target_tol, mt_length)
                and size_ok
            ):
                duplicate = True
                break
        if duplicate:
            continue
        targets.append(
            {
                "name": name,
                "source": "analysis.known_sequence_searches",
                "match_reason": "configured_sequence_search_target_match",
                "left_breakpoint": left,
                "right_breakpoint": right,
                "deleted_size": size,
                "breakpoint_tolerance_bp": breakpoint_tolerance,
                "size_tolerance_bp": size_tolerance,
            }
        )
    return targets


def known_deletion_match(left: int, right: int, size: int, config: dict, mt_length: int) -> tuple[str, str]:
    matches = []
    reasons = []
    for item in configured_deletion_targets(config, mt_length):
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        tol = int(item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", 50)))
        size_tol = int(item.get("size_tolerance_bp", tol))
        left_ok = "left_breakpoint" not in item or pos_within_circular_window(left, int(item["left_breakpoint"]), tol, mt_length)
        right_ok = "right_breakpoint" not in item or pos_within_circular_window(right, int(item["right_breakpoint"]), tol, mt_length)
        size_ok = "deleted_size" not in item or abs(int(size) - int(item["deleted_size"])) <= size_tol
        if left_ok and right_ok and size_ok:
            matches.append(name)
            reasons.append(str(item.get("match_reason", "configured_deletion_target_match")))
    if matches:
        return ";".join(dict.fromkeys(matches)), ";".join(dict.fromkeys(reasons))
    return "", ""
