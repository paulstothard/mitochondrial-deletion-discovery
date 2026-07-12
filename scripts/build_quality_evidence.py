#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from annotate_junctions import append_configured_regions, apply_feature_aliases, biological_features
from circular_deletions import circular_position_distance, directed_breakpoints
from common import read_tsv, read_yaml, write_tsv
from consolidate_deletions import circular_median, cluster_rows, split_direction_conflicts


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ALIGNED_OPS = {"M", "=", "X"}
QUERY_OPS = {"M", "I", "S", "H", "=", "X"}


def clean_read_id(value: object) -> str:
    return str(value or "").removesuffix("/1").removesuffix("/2")


def yes(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def numeric(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def integer(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def cigar_query_metrics(cigar: str) -> dict[str, int]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(str(cigar or ""))]
    query_length = sum(length for length, op in operations if op in QUERY_OPS)
    leading_clip = 0
    for length, op in operations:
        if op not in {"S", "H"}:
            break
        leading_clip += length
    aligned_query = sum(length for length, op in operations if op in ALIGNED_OPS or op == "I")
    aligned_reference = sum(length for length, op in operations if op in ALIGNED_OPS or op in {"D", "N"})
    return {
        "query_length": query_length,
        "query_start": leading_clip,
        "query_end": leading_clip + aligned_query,
        "aligned_query_length": aligned_query,
        "aligned_reference_length": aligned_reference,
    }


def query_pair_metrics(first: dict, second: dict) -> dict[str, object]:
    ordered = sorted((first, second), key=lambda item: (item["query_start"], item["query_end"]))
    left, right = ordered
    overlap = max(0, int(left["query_end"]) - int(right["query_start"]))
    gap = max(0, int(right["query_start"]) - int(left["query_end"]))
    union = (
        int(left["query_end"])
        - int(left["query_start"])
        + int(right["query_end"])
        - int(right["query_start"])
        - overlap
    )
    query_length = max(int(first["query_length"]), int(second["query_length"]), 1)
    return {
        "query_first_start": int(left["query_start"]),
        "query_first_end": int(left["query_end"]),
        "query_second_start": int(right["query_start"]),
        "query_second_end": int(right["query_end"]),
        "query_overlap_bp": overlap,
        "query_gap_bp": gap,
        "query_union_coverage": min(1.0, max(0.0, union / query_length)),
        "query_union_aligned_length": union,
        "total_aligned_query_length": (
            int(left["query_end"])
            - int(left["query_start"])
            + int(right["query_end"])
            - int(right["query_start"])
        ),
        "read_length": query_length,
        "query_segments_adjacent": "yes" if gap == 0 else "no",
        "query_junction_boundary_distance_from_read_end_min": min(
            int(left["query_end"]),
            query_length - int(left["query_end"]),
            int(right["query_start"]),
            query_length - int(right["query_start"]),
        ),
    }


def inferred_library_layout(config: dict) -> str:
    dataset = config.get("dataset", {}) or {}
    explicit = str(dataset.get("library_layout", "")).strip().lower()
    if explicit in {"single", "paired"}:
        return explicit
    strategy = str(dataset.get("library_strategy", "")).lower()
    if "paired" in strategy:
        return "paired"
    if "single" in strategy:
        return "single"
    return "unknown"


def alignment_segment_count(row: dict) -> int:
    counts = [2]
    for key in ("left_sa_tag", "right_sa_tag"):
        entries = [value for value in str(row.get(key, "")).split(";") if value]
        if entries:
            counts.append(1 + len(entries))
    return max(counts)


def star_sample(path: str | Path) -> str:
    name = Path(path).name
    suffix = ".Chimeric.out.junction"
    return name[: -len(suffix)] if name.endswith(suffix) else name.split(".", 1)[0]


def load_gene_features(path: str, config: dict, mt_length: int) -> list[dict]:
    features = pd.read_csv(path, sep="\t")
    features = apply_feature_aliases(biological_features(append_configured_regions(features, config, mt_length)), config)
    rows = []
    for _, row in features.iterrows():
        if str(row.get("feature_type", "")) == "region":
            continue
        name = str(row.get("gene_name") or row.get("gene_id") or "").strip()
        if not name:
            continue
        rows.append({"name": name, "start": int(row["start"]), "end": int(row["end"])})
    return rows


def segment_gene_hits(start: int, cigar: str, features: list[dict]) -> list[dict]:
    reference_length = cigar_query_metrics(cigar)["aligned_reference_length"]
    end = int(start) + max(0, reference_length - 1)
    hits = []
    for feature in features:
        overlap = max(0, min(end, feature["end"]) - max(int(start), feature["start"]) + 1)
        if overlap:
            hits.append({**feature, "overlap": overlap})
    return sorted(hits, key=lambda item: (-item["overlap"], item["name"], item["start"]))


def genes_overlap_substantially(left: dict, right: dict, threshold: float = 0.5) -> bool:
    overlap = max(0, min(left["end"], right["end"]) - max(left["start"], right["start"]) + 1)
    shorter = min(left["end"] - left["start"] + 1, right["end"] - right["start"] + 1)
    return shorter > 0 and overlap / shorter > threshold


def configured_expected_pair(left: str, right: str, config: dict) -> bool:
    left_lower = left.lower()
    right_lower = right.lower()
    for pair in config.get("annotations", {}).get("expected_adjacent_transcripts", []) or []:
        expected_left = str(pair.get("left_feature", "")).lower()
        expected_right = str(pair.get("right_feature", "")).lower()
        if (left_lower, right_lower) == (expected_left, expected_right):
            return True
        if pair.get("bidirectional") and (left_lower, right_lower) == (expected_right, expected_left):
            return True
    return False


def parse_star_chimeric_line(
    line: str,
    sample: str,
    species: str,
    mt_names: set[str],
    mt_length: int,
    min_anchor: int,
    min_deletion_size: int,
    max_deletion_size: int,
    max_query_overlap: int,
    max_query_gap: int,
    require_same_orientation: bool = True,
    gene_features: list[dict] | None = None,
    config: dict | None = None,
    require_gene_anchors: bool = False,
    exclude_same_gene: bool = False,
) -> dict | None:
    if not line.strip() or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 14 or parts[0] in {"chr_donorA", "chromosome"}:
        return None
    if parts[0] not in mt_names or parts[3] not in mt_names:
        return None
    try:
        raw_left = int(parts[1])
        raw_right = int(parts[4])
    except ValueError:
        return None
    left_cigar = parts[11]
    right_cigar = parts[13]
    left_query = cigar_query_metrics(left_cigar)
    right_query = cigar_query_metrics(right_cigar)
    pair = query_pair_metrics(left_query, right_query)
    directed = directed_breakpoints(raw_left, raw_right, mt_length)
    left_anchor = int(left_query["aligned_query_length"])
    right_anchor = int(right_query["aligned_query_length"])
    reasons = []
    config = config or {}
    if require_same_orientation and parts[2] != parts[5]:
        reasons.append("incompatible_reference_orientation")
    if left_anchor < min_anchor:
        reasons.append("left_anchor_short")
    if right_anchor < min_anchor:
        reasons.append("right_anchor_short")
    if int(directed["deleted_size"]) < min_deletion_size:
        reasons.append("deletion_too_small")
    if int(directed["deleted_size"]) > max_deletion_size:
        reasons.append("deletion_too_large")
    if int(pair["query_overlap_bp"]) > max_query_overlap:
        reasons.append("query_overlap_too_large")
    if int(pair["query_gap_bp"]) > max_query_gap:
        reasons.append("query_gap_too_large")
    left_gene_hits = segment_gene_hits(integer(parts[10]), left_cigar, gene_features or [])
    right_gene_hits = segment_gene_hits(integer(parts[12]), right_cigar, gene_features or [])
    left_gene = left_gene_hits[0]["name"] if left_gene_hits else ""
    right_gene = right_gene_hits[0]["name"] if right_gene_hits else ""
    if require_gene_anchors and (not left_gene_hits or not right_gene_hits):
        reasons.append("missing_annotated_gene_anchor")
    common_genes = {hit["name"] for hit in left_gene_hits} & {hit["name"] for hit in right_gene_hits}
    if exclude_same_gene and common_genes:
        reasons.append("same_gene_chimeric_alignment")
    if left_gene_hits and right_gene_hits and genes_overlap_substantially(left_gene_hits[0], right_gene_hits[0]):
        reasons.append("substantially_overlapping_gene_anchors")
    if left_gene and right_gene and configured_expected_pair(left_gene, right_gene, config):
        reasons.append("configured_expected_transcript_pair")
    this_score = integer(parts[17]) if len(parts) > 17 else 0
    best_score = integer(parts[18]) if len(parts) > 18 else this_score
    row = {
        "sample": sample,
        "species": species,
        "read_id": clean_read_id(parts[9]),
        "physical_observation_id": clean_read_id(parts[9]),
        "left_breakpoint": raw_left,
        "right_breakpoint": raw_right,
        "deleted_size": directed["deleted_size"],
        "wraps_origin": directed["wraps_origin"],
        "deleted_interval": directed["deleted_interval"],
        "complement_deleted_size": directed["complement_deleted_size"],
        "complement_wraps_origin": directed["complement_wraps_origin"],
        "arc_assignment_method": "alignment_directed",
        "direction_status": "directed",
        "breakpoint_pair_id": "",
        "source": "STAR_Chimeric.out.junction",
        "evidence_source": "star_chimeric",
        "rotation_name": "full_genome",
        "rotation_start": 1,
        "raw_left_breakpoint": raw_left,
        "raw_right_breakpoint": raw_right,
        "reported_left_breakpoint": raw_left,
        "reported_right_breakpoint": raw_right,
        "left_anchor_length": left_anchor,
        "right_anchor_length": right_anchor,
        "min_anchor_length": min(left_anchor, right_anchor),
        "strand": f"{parts[2]}:{parts[5]}",
        "left_cigar": left_cigar,
        "right_cigar": right_cigar,
        "left_mapq": "",
        "right_mapq": "",
        "min_mapq": "",
        "left_alignment_score": this_score,
        "right_alignment_score": this_score,
        "left_edit_distance": "",
        "right_edit_distance": "",
        "star_num_chimeric_alignments": integer(parts[14]) if len(parts) > 14 else "",
        "star_max_possible_alignment_score": integer(parts[15]) if len(parts) > 15 else "",
        "star_non_chimeric_alignment_score": integer(parts[16]) if len(parts) > 16 else "",
        "star_chimeric_alignment_score": this_score,
        "star_best_chimeric_alignment_score": best_score,
        "star_chimeric_score_delta": best_score - this_score,
        "star_repeat_left_length": integer(parts[7]) if len(parts) > 7 else 0,
        "star_repeat_right_length": integer(parts[8]) if len(parts) > 8 else 0,
        "star_left_gene": left_gene,
        "star_right_gene": right_gene,
        "star_gene_pair_label": f"{left_gene}--{right_gene}" if left_gene and right_gene else "",
        "star_left_gene_hits": ";".join(hit["name"] for hit in left_gene_hits),
        "star_right_gene_hits": ";".join(hit["name"] for hit in right_gene_hits),
        "total_alignment_segments": 2,
        "alignment_pattern_id": f"star:{raw_left}:{raw_right}:{left_cigar}:{right_cigar}",
        "left_clipped_fraction": (
            1.0 - left_anchor / max(1, int(left_query["query_length"]))
        ),
        "right_clipped_fraction": (
            1.0 - right_anchor / max(1, int(right_query["query_length"]))
        ),
        "star_junction_repeat_length": min(
            integer(parts[7]) if len(parts) > 7 else 0,
            integer(parts[8]) if len(parts) > 8 else 0,
        ),
        "nuclear_competition_status": "not_available_in_chimeric_record",
        "left_is_primary": "unknown",
        "right_is_primary": "unknown",
        "left_is_secondary": "unknown",
        "right_is_secondary": "unknown",
        "left_is_supplementary": "unknown",
        "right_is_supplementary": "unknown",
        "primary_chain_evidence": "unknown",
        "secondary_only_evidence": "unknown",
        "filter_status": "pass" if not reasons else "fail",
        "filter_reason": ";".join(reasons),
        **pair,
    }
    return row


def minimap_observation_row(row: dict) -> dict:
    out = dict(row)
    out["read_id"] = clean_read_id(row.get("read_id"))
    out["physical_observation_id"] = out["read_id"]
    out["source"] = row.get("source") or "minimap2_split_alignment"
    out["evidence_source"] = "minimap2_remap"
    left_query = {
        "query_start": integer(row.get("query_first_start")),
        "query_end": integer(row.get("query_first_end")),
        "query_length": 0,
    }
    right_query = {
        "query_start": integer(row.get("query_second_start")),
        "query_end": integer(row.get("query_second_end")),
        "query_length": 0,
    }
    cigar_lengths = [cigar_query_metrics(row.get(name, ""))["query_length"] for name in ("left_cigar", "right_cigar")]
    read_length = max(cigar_lengths + [left_query["query_end"], right_query["query_end"], 1])
    left_query["query_length"] = read_length
    right_query["query_length"] = read_length
    pair = query_pair_metrics(left_query, right_query)
    for key, value in pair.items():
        out.setdefault(key, value)
    out["read_length"] = read_length
    left_aligned = integer(row.get("left_anchor_length"))
    right_aligned = integer(row.get("right_anchor_length"))
    edit_total = integer(row.get("left_edit_distance")) + integer(row.get("right_edit_distance"))
    aligned_total = max(1, left_aligned + right_aligned)
    out["estimated_alignment_error_rate"] = edit_total / aligned_total
    out["left_clipped_fraction"] = 1.0 - left_aligned / read_length
    out["right_clipped_fraction"] = 1.0 - right_aligned / read_length
    out["total_alignment_segments"] = alignment_segment_count(row)
    out["alignment_pattern_id"] = ":".join(
        [
            "minimap2",
            str(row.get("rotation_name", "")),
            str(row.get("raw_left_breakpoint", row.get("left_breakpoint", ""))),
            str(row.get("raw_right_breakpoint", row.get("right_breakpoint", ""))),
            str(row.get("left_cigar", "")),
            str(row.get("right_cigar", "")),
        ]
    )
    out["nuclear_competition_status"] = str(row.get("nuclear_competition_status", "not_propagated"))
    left_primary = yes(row.get("left_is_primary"))
    right_primary = yes(row.get("right_is_primary"))
    left_secondary = yes(row.get("left_is_secondary"))
    right_secondary = yes(row.get("right_is_secondary"))
    out["primary_chain_evidence"] = "yes" if left_primary or right_primary else "no"
    out["secondary_only_evidence"] = "yes" if left_secondary and right_secondary else "no"
    out["filter_status"] = "pass"
    out["filter_reason"] = ""
    return out


def rows_close(first: dict, second: dict, slop: int, mt_length: int) -> bool:
    return (
        circular_position_distance(first["left_breakpoint"], second["left_breakpoint"], mt_length) <= slop
        and circular_position_distance(first["right_breakpoint"], second["right_breakpoint"], mt_length) <= slop
    )


def row_quality_flags(row: dict) -> list[str]:
    flags = []
    if row.get("evidence_source") == "minimap2_remap" and integer(row.get("min_mapq")) == 0:
        flags.append("min_mapq_zero")
    if row.get("secondary_only_evidence") == "yes":
        flags.append("secondary_only")
    if numeric(row.get("estimated_alignment_error_rate"), 0.0) > 0.15:
        flags.append("high_alignment_error")
    if integer(row.get("star_num_chimeric_alignments"), 1) > 1:
        flags.append("star_multimapping")
    if integer(row.get("query_gap_bp")) > 0:
        flags.append("query_gap")
    if integer(row.get("query_overlap_bp")) > 0:
        flags.append("query_overlap")
    if integer(row.get("total_alignment_segments"), 2) > 2:
        flags.append("complex_alignment_chain")
    return flags


def aggregate_hypothesis(rows: list[dict], mt_length: int) -> dict:
    sources = sorted({row["evidence_source"] for row in rows})
    source_positions = []
    for source in sources:
        source_rows = [row for row in rows if row["evidence_source"] == source]
        source_positions.append(
            (
                circular_median([integer(row["left_breakpoint"]) for row in source_rows], mt_length),
                circular_median([integer(row["right_breakpoint"]) for row in source_rows], mt_length),
            )
        )
    left = circular_median([value[0] for value in source_positions], mt_length)
    right = circular_median([value[1] for value in source_positions], mt_length)
    directed = directed_breakpoints(left, right, mt_length)
    flags = sorted({flag for row in rows for flag in row_quality_flags(row)})
    rotations = sorted({str(row.get("rotation_name", "")) for row in rows if row.get("rotation_name")})
    min_mapq_values = [integer(row.get("min_mapq")) for row in rows if str(row.get("min_mapq", "")) != ""]
    error_values = [numeric(row.get("estimated_alignment_error_rate")) for row in rows if str(row.get("estimated_alignment_error_rate", "")) != ""]
    union_values = [numeric(row.get("query_union_coverage")) for row in rows if str(row.get("query_union_coverage", "")) != ""]
    anchors = [integer(row.get("min_anchor_length"), min(integer(row.get("left_anchor_length")), integer(row.get("right_anchor_length")))) for row in rows]
    star_gene_pairs = [
        str(row.get("star_gene_pair_label", ""))
        for row in rows
        if row.get("evidence_source") == "star_chimeric" and str(row.get("star_gene_pair_label", ""))
    ]
    modal_star_gene_pair = Counter(star_gene_pairs).most_common(1)[0][0] if star_gene_pairs else ""
    alignment_patterns = sorted({str(row.get("alignment_pattern_id", "")) for row in rows if row.get("alignment_pattern_id")})
    def joined(field: str) -> str:
        return ";".join(sorted({str(row.get(field, "")) for row in rows if str(row.get(field, "")).strip()}))

    def minimum_numeric(field: str) -> int | float | str:
        values = [numeric(row.get(field)) for row in rows if str(row.get(field, "")).strip()]
        return min(values) if values else ""

    def maximum_numeric(field: str) -> int | float | str:
        values = [numeric(row.get(field)) for row in rows if str(row.get(field, "")).strip()]
        return max(values) if values else ""

    def any_yes(field: str) -> str:
        values = {str(row.get(field, "")).lower() for row in rows if str(row.get(field, "")).strip()}
        if "yes" in values:
            return "yes"
        if "no" in values:
            return "no"
        return "unknown"

    sample = rows[0]["sample"]
    read_id = clean_read_id(rows[0]["read_id"])
    return {
        "sample": sample,
        "species": rows[0]["species"],
        "read_id": read_id,
        "physical_observation_id": read_id,
        "left_breakpoint": left,
        "right_breakpoint": right,
        "deleted_size": directed["deleted_size"],
        "wraps_origin": directed["wraps_origin"],
        "deleted_interval": directed["deleted_interval"],
        "complement_deleted_size": directed["complement_deleted_size"],
        "complement_wraps_origin": directed["complement_wraps_origin"],
        "arc_assignment_method": "alignment_directed",
        "direction_status": "directed",
        "evidence_sources": ";".join(sources),
        "source": ";".join(sources),
        "star_support": "yes" if "star_chimeric" in sources else "no",
        "minimap2_support": "yes" if "minimap2_remap" in sources else "no",
        "both_callers_support": "yes" if len(sources) > 1 else "no",
        "strand": joined("strand"),
        "source_record_count": len(rows),
        "star_source_record_count": sum(row["evidence_source"] == "star_chimeric" for row in rows),
        "minimap2_source_record_count": sum(row["evidence_source"] == "minimap2_remap" for row in rows),
        "star_gene_pair_label": modal_star_gene_pair,
        "star_gene_pair_labels": ";".join(sorted(set(star_gene_pairs))),
        "rotation_name": ";".join(rotations),
        "rotation_support": ";".join(rotations),
        "rotation_count": len(rotations),
        "rotation_agreement": "multiple_rotations" if len([value for value in rotations if value != "full_genome"]) > 1 else "single_rotation",
        "left_anchor_length": min(integer(row.get("left_anchor_length")) for row in rows),
        "right_anchor_length": min(integer(row.get("right_anchor_length")) for row in rows),
        "left_aligned_length": min(integer(row.get("left_anchor_length")) for row in rows),
        "right_aligned_length": min(integer(row.get("right_anchor_length")) for row in rows),
        "left_mapq": minimum_numeric("left_mapq"),
        "right_mapq": minimum_numeric("right_mapq"),
        "left_alignment_score": minimum_numeric("left_alignment_score"),
        "right_alignment_score": minimum_numeric("right_alignment_score"),
        "left_edit_distance": maximum_numeric("left_edit_distance"),
        "right_edit_distance": maximum_numeric("right_edit_distance"),
        "left_cigar": joined("left_cigar"),
        "right_cigar": joined("right_cigar"),
        "left_is_primary": any_yes("left_is_primary"),
        "right_is_primary": any_yes("right_is_primary"),
        "left_is_secondary": any_yes("left_is_secondary"),
        "right_is_secondary": any_yes("right_is_secondary"),
        "left_is_supplementary": any_yes("left_is_supplementary"),
        "right_is_supplementary": any_yes("right_is_supplementary"),
        "min_anchor_length": min(anchors) if anchors else 0,
        "read_length": max(integer(row.get("read_length")) for row in rows),
        "query_overlap_bp": max(integer(row.get("query_overlap_bp")) for row in rows),
        "query_gap_bp": max(integer(row.get("query_gap_bp")) for row in rows),
        "query_union_coverage": min(union_values) if union_values else "",
        "query_segments_adjacent": "yes" if all(row.get("query_segments_adjacent") == "yes" for row in rows) else "no",
        "query_union_aligned_length": min(
            [integer(row.get("query_union_aligned_length")) for row in rows if str(row.get("query_union_aligned_length", "")) != ""],
            default="",
        ),
        "total_aligned_query_length": min(
            [integer(row.get("total_aligned_query_length")) for row in rows if str(row.get("total_aligned_query_length", "")) != ""],
            default="",
        ),
        "query_junction_boundary_distance_from_read_end_min": min(
            [integer(row.get("query_junction_boundary_distance_from_read_end_min")) for row in rows if str(row.get("query_junction_boundary_distance_from_read_end_min", "")) != ""],
            default="",
        ),
        "left_clipped_fraction": max(
            [numeric(row.get("left_clipped_fraction")) for row in rows if str(row.get("left_clipped_fraction", "")) != ""],
            default="",
        ),
        "right_clipped_fraction": max(
            [numeric(row.get("right_clipped_fraction")) for row in rows if str(row.get("right_clipped_fraction", "")) != ""],
            default="",
        ),
        "total_alignment_segments": max(integer(row.get("total_alignment_segments"), 2) for row in rows),
        "alignment_pattern_id": ";".join(alignment_patterns),
        "alignment_pattern_count": len(alignment_patterns),
        "source_coordinate_pairs": ";".join(
            sorted(
                {
                    f"{row.get('evidence_source', '')}:{row.get('raw_left_breakpoint', row.get('left_breakpoint', ''))}:{row.get('raw_right_breakpoint', row.get('right_breakpoint', ''))}"
                    for row in rows
                }
            )
        ),
        "star_junction_repeat_length": max(integer(row.get("star_junction_repeat_length")) for row in rows),
        "nuclear_competition_status": ";".join(
            sorted({str(row.get("nuclear_competition_status", "")) for row in rows if row.get("nuclear_competition_status")})
        ),
        "min_mapq": min(min_mapq_values) if min_mapq_values else "",
        "estimated_alignment_error_rate": max(error_values) if error_values else "",
        "primary_chain_evidence": "yes" if any(row.get("primary_chain_evidence") == "yes" for row in rows) else "unknown" if all(row.get("primary_chain_evidence") == "unknown" for row in rows) else "no",
        "secondary_only_evidence": "yes" if rows and all(row.get("secondary_only_evidence") == "yes" for row in rows) else "no",
        "library_layout": joined("library_layout"),
        "observation_unit": joined("observation_unit"),
        "paired_end_collapsed": joined("paired_end_collapsed"),
        "mate_context_status": joined("mate_context_status"),
        "both_mates_support_same_cluster": joined("both_mates_support_same_cluster"),
        "non_supporting_mate_near_breakpoint": joined("non_supporting_mate_near_breakpoint"),
        "mate_mapping_class": joined("mate_mapping_class"),
        "mate_placement_discordant": joined("mate_placement_discordant"),
        "observation_quality_flags": ";".join(flags),
        "filter_status": "pass",
        "filter_reason": "",
    }


def collapse_physical_observations(rows: list[dict], slop: int, mt_length: int) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], clean_read_id(row["read_id"]))].append(row)
    observations = []
    for key in sorted(grouped):
        hypotheses: list[list[dict]] = []
        for row in sorted(grouped[key], key=lambda item: (integer(item["left_breakpoint"]), integer(item["right_breakpoint"]), item["evidence_source"])):
            match = next((group for group in hypotheses if rows_close(row, group[0], slop, mt_length)), None)
            if match is None:
                hypotheses.append([row])
            else:
                match.append(row)
        count = len(hypotheses)
        for group in hypotheses:
            observation = aggregate_hypothesis(group, mt_length)
            observation["deletion_hypotheses_from_read"] = count
            if count > 1:
                flags = {value for value in observation["observation_quality_flags"].split(";") if value}
                flags.add("multiple_deletion_hypotheses")
                observation["observation_quality_flags"] = ";".join(sorted(flags))
            observations.append(observation)
    return observations


def read_star_rows(paths: list[str], args: argparse.Namespace) -> list[dict]:
    rows = []
    mt_names = {value for value in args.mt_contig_names.split(",") if value}
    config = read_yaml(args.config)
    gene_features = load_gene_features(args.features, config, args.mt_length)
    for path in paths:
        sample = star_sample(path)
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                row = parse_star_chimeric_line(
                    line,
                    sample,
                    args.species,
                    mt_names,
                    args.mt_length,
                    args.star_min_anchor_length,
                    args.min_deletion_size,
                    args.max_deletion_size,
                    args.star_max_query_overlap_bp,
                    args.star_max_query_gap_bp,
                    args.require_same_orientation,
                    gene_features,
                    config,
                    args.star_require_gene_anchors,
                    args.star_exclude_same_gene,
                )
                if row is not None:
                    rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--mt-length", type=int, required=True)
    parser.add_argument("--mt-contig-names", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--breakpoint-slop-bp", type=int, required=True)
    parser.add_argument("--min-deletion-size", type=int, required=True)
    parser.add_argument("--max-deletion-size", type=int, required=True)
    parser.add_argument("--star-min-anchor-length", type=int, default=12)
    parser.add_argument("--star-max-query-overlap-bp", type=int, default=20)
    parser.add_argument("--star-max-query-gap-bp", type=int, default=20)
    parser.add_argument("--require-same-orientation", action="store_true")
    parser.add_argument("--star-require-gene-anchors", action="store_true")
    parser.add_argument("--star-exclude-same-gene", action="store_true")
    parser.add_argument("--ambiguous-direction-policy", choices=["exclude", "include"], default="exclude")
    parser.add_argument("--result-schema-version", default="3.0-quality-evidence")
    parser.add_argument("--minimap-reads", nargs="*", default=[])
    parser.add_argument("--star-junctions", nargs="*", default=[])
    parser.add_argument("--out-source-candidates", required=True)
    parser.add_argument("--out-observations", required=True)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-id-map", required=True)
    parser.add_argument("--out-ambiguous", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    layout = inferred_library_layout(config)

    minimap_rows = [minimap_observation_row(row) for path in args.minimap_reads for row in read_tsv(path)]
    star_rows = read_star_rows(args.star_junctions, args)
    source_candidates = minimap_rows + star_rows
    for row in source_candidates:
        row["library_layout"] = layout
        row["observation_unit"] = "sequenced_fragment" if layout == "paired" else "read"
        row["paired_end_collapsed"] = "yes" if layout == "paired" else "not_applicable"
        mate_status = "not_available_from_retained_intermediates" if layout == "paired" else "not_applicable"
        row["mate_context_status"] = mate_status
        row["both_mates_support_same_cluster"] = mate_status
        row["non_supporting_mate_near_breakpoint"] = mate_status
        row["mate_mapping_class"] = mate_status
        row["mate_placement_discordant"] = mate_status
    passing = [row for row in source_candidates if row.get("filter_status") == "pass"]
    passing, ambiguous = split_direction_conflicts(passing, args.mt_length, args.breakpoint_slop_bp)
    if args.ambiguous_direction_policy == "include":
        passing.extend(ambiguous)
    observations = collapse_physical_observations(passing, args.breakpoint_slop_bp, args.mt_length)
    clustered_observations, clusters, id_map = cluster_rows(
        observations,
        args.breakpoint_slop_bp,
        1,
        args.mt_length,
        result_schema_version=args.result_schema_version,
    )
    source_fields = [
        "sample", "species", "read_id", "physical_observation_id", "evidence_source", "source",
        "left_breakpoint", "right_breakpoint", "deleted_size", "wraps_origin", "strand",
        "left_anchor_length", "right_anchor_length", "min_anchor_length", "read_length",
        "total_aligned_query_length", "query_union_aligned_length", "query_overlap_bp", "query_gap_bp",
        "query_union_coverage", "query_segments_adjacent", "query_junction_boundary_distance_from_read_end_min",
        "left_clipped_fraction", "right_clipped_fraction", "total_alignment_segments",
        "min_mapq", "primary_chain_evidence", "secondary_only_evidence",
        "rotation_name", "library_layout", "observation_unit", "paired_end_collapsed", "mate_context_status",
        "both_mates_support_same_cluster", "non_supporting_mate_near_breakpoint", "mate_mapping_class", "mate_placement_discordant",
        "nuclear_competition_status", "filter_status", "filter_reason",
    ]
    source_fields.extend(sorted({key for row in source_candidates for key in row if key not in source_fields}))
    observation_fields = [
        "exact_deletion_id", "junction_id", "sample", "species", "read_id", "physical_observation_id",
        "left_breakpoint", "right_breakpoint", "deleted_size", "read_left_breakpoint", "read_right_breakpoint",
        "wraps_origin", "deleted_interval", "complement_deleted_size", "complement_wraps_origin",
        "evidence_sources", "star_support", "minimap2_support", "both_callers_support",
        "source_record_count", "star_source_record_count", "minimap2_source_record_count",
        "star_gene_pair_label", "star_gene_pair_labels",
        "strand", "left_anchor_length", "right_anchor_length", "left_aligned_length", "right_aligned_length", "min_anchor_length", "read_length",
        "left_mapq", "right_mapq", "left_alignment_score", "right_alignment_score",
        "left_edit_distance", "right_edit_distance", "left_cigar", "right_cigar",
        "left_is_primary", "right_is_primary", "left_is_secondary", "right_is_secondary",
        "left_is_supplementary", "right_is_supplementary",
        "total_aligned_query_length", "query_union_aligned_length", "query_overlap_bp", "query_gap_bp",
        "query_union_coverage", "query_segments_adjacent", "query_junction_boundary_distance_from_read_end_min",
        "left_clipped_fraction", "right_clipped_fraction", "total_alignment_segments",
        "min_mapq", "estimated_alignment_error_rate", "primary_chain_evidence", "secondary_only_evidence",
        "rotation_support", "rotation_count", "rotation_agreement", "deletion_hypotheses_from_read",
        "library_layout", "observation_unit", "paired_end_collapsed", "mate_context_status",
        "both_mates_support_same_cluster", "non_supporting_mate_near_breakpoint", "mate_mapping_class", "mate_placement_discordant",
        "nuclear_competition_status", "alignment_pattern_id", "alignment_pattern_count",
        "source_coordinate_pairs",
        "observation_quality_flags", "direction_status", "result_schema_version",
    ]
    observation_fields.extend(sorted({key for row in clustered_observations for key in row if key not in observation_fields}))
    write_tsv(args.out_source_candidates, source_candidates, source_fields)
    write_tsv(args.out_observations, clustered_observations, observation_fields)
    write_tsv(args.out_clusters, clusters)
    write_tsv(args.out_id_map, id_map)
    write_tsv(args.out_ambiguous, ambiguous, source_fields)
    counts = Counter(row.get("evidence_source", "unknown") for row in source_candidates)
    summary = [
        {"metric": "minimap2_source_candidates", "value": counts["minimap2_remap"]},
        {"metric": "star_source_candidates", "value": counts["star_chimeric"]},
        {"metric": "source_candidates_passing", "value": len(passing)},
        {"metric": "ambiguous_direction_records", "value": len(ambiguous)},
        {"metric": "canonical_distinct_observations", "value": len(clustered_observations)},
        {"metric": "canonical_deletion_clusters", "value": len(clusters)},
        {"metric": "both_caller_observations", "value": sum(row.get("both_callers_support") == "yes" for row in clustered_observations)},
    ]
    write_tsv(args.out_summary, summary, ["metric", "value"])


if __name__ == "__main__":
    main()
