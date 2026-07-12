#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from statistics import median

from circular_deletions import circular_position_distance

from common import read_tsv, read_yaml, write_tsv


DEFAULT_PROFILES = {
    "stringent": ["strong"],
    "standard": ["strong", "supported"],
    "exploratory": ["strong", "supported", "review"],
}


def integer(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def numeric(value: object) -> float | None:
    try:
        if str(value).strip() == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def configured_profiles(config: dict) -> dict[str, list[str]]:
    configured = config.get("quality", {}).get("report_profiles", {}) or {}
    profiles = {}
    for name, value in configured.items():
        tiers = value.get("include_tiers", []) if isinstance(value, dict) else value
        tiers = [str(tier).strip().lower() for tier in tiers or [] if str(tier).strip()]
        if tiers:
            profiles[str(name)] = tiers
    return profiles or dict(DEFAULT_PROFILES)


def feature_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "intergenic"
    return text.split(";", 1)[0].split(":", 1)[0].strip() or "intergenic"


def gene_pair_label(cluster: dict) -> str:
    left = feature_name(cluster.get("nearest_left_feature") or cluster.get("left_feature_overlap"))
    right = feature_name(cluster.get("nearest_right_feature") or cluster.get("right_feature_overlap"))
    return f"{left}--{right}"


def modal_value(values: list[str]) -> str:
    counts = Counter(value for value in values if value)
    if not counts:
        return ""
    return sorted(counts, key=lambda value: (-counts[value], value))[0]


def proportion(values: list[bool]) -> float | str:
    return sum(values) / len(values) if values else ""


def joined_flags(observations: list[dict], cluster: dict) -> list[str]:
    flags = {
        flag
        for row in observations
        for flag in str(row.get("observation_quality_flags", "")).split(";")
        if flag
    }
    if len(observations) == 1:
        flags.add("singleton_cluster")
    if cluster.get("junction_interpretation") == "expected_transcript_junction":
        flags.add("expected_transcript_junction")
    if any(str(row.get("direction_status", "directed")) != "directed" for row in observations):
        flags.add("direction_ambiguity")
    return sorted(flags)


def cluster_quality_tier(cluster: dict, observations: list[dict], config: dict) -> tuple[str, str]:
    quality = config.get("quality", {}) or {}
    minimum_supported = int(quality.get("minimum_supported_observations", 2))
    minimum_strong = int(quality.get("minimum_strong_observations", 2))
    exclude_expected = bool(config.get("junctions", {}).get("exclude_expected_transcript_junctions", True))
    if exclude_expected and cluster.get("junction_interpretation") == "expected_transcript_junction":
        return "rejected", "configured_expected_transcript_junction"
    if any(str(row.get("direction_status", "directed")) != "directed" for row in observations):
        return "rejected", "unresolved_direction_ambiguity"
    total = len(observations)
    if total < minimum_supported:
        return "review", "fewer_than_minimum_supported_observations"
    per_sample = Counter(row["sample"] for row in observations)
    both_callers = sum(row.get("both_callers_support") == "yes" for row in observations)
    if total >= minimum_strong and (max(per_sample.values(), default=0) >= 2 or both_callers > 0):
        return "strong", "replicated_within_sample_or_cross_caller_corroboration"
    return "supported", "multi_observation_support_without_strong_corroboration"


def summarize_cluster(cluster: dict, observations: list[dict], config: dict) -> dict:
    out = dict(cluster)
    per_sample = Counter(row["sample"] for row in observations)
    sources = {
        source
        for row in observations
        for source in str(row.get("evidence_sources", "")).split(";")
        if source
    }
    star_count = sum(row.get("star_support") == "yes" for row in observations)
    minimap_count = sum(row.get("minimap2_support") == "yes" for row in observations)
    both_count = sum(row.get("both_callers_support") == "yes" for row in observations)
    anchors = [integer(row.get("min_anchor_length")) for row in observations if integer(row.get("min_anchor_length")) > 0]
    errors = [value for row in observations if (value := numeric(row.get("estimated_alignment_error_rate"))) is not None]
    coverages = [value for row in observations if (value := numeric(row.get("query_union_coverage"))) is not None]
    segment_counts = [integer(row.get("total_alignment_segments"), 2) for row in observations]
    observed_left = [integer(row.get("read_left_breakpoint", row.get("left_breakpoint"))) for row in observations]
    observed_right = [integer(row.get("read_right_breakpoint", row.get("right_breakpoint"))) for row in observations]
    mt_length = int((config.get("references", {}).get(str(out.get("species", "")), {}) or {}).get("mt_length", 0) or 0)
    modal_pairs = Counter(zip(observed_left, observed_right))
    modal_pair, modal_pair_support = modal_pairs.most_common(1)[0] if modal_pairs else ((0, 0), 0)
    if mt_length:
        left_deviations = [circular_position_distance(value, integer(out.get("left_breakpoint")), mt_length) for value in observed_left]
        right_deviations = [circular_position_distance(value, integer(out.get("right_breakpoint")), mt_length) for value in observed_right]
    else:
        left_deviations = [abs(value - integer(out.get("left_breakpoint"))) for value in observed_left]
        right_deviations = [abs(value - integer(out.get("right_breakpoint"))) for value in observed_right]
    rotation_sets = [set(value for value in str(row.get("rotation_support", "")).split(";") if value and value != "full_genome") for row in observations]
    both_rotation_values = [len(value) > 1 for value in rotation_sets]
    star_gene_pairs = [str(row.get("star_gene_pair_label", "")) for row in observations if str(row.get("star_gene_pair_label", ""))]
    primary_values = [row.get("primary_chain_evidence") == "yes" for row in observations if row.get("primary_chain_evidence") in {"yes", "no"}]
    secondary_values = [row.get("secondary_only_evidence") == "yes" for row in observations if row.get("secondary_only_evidence") in {"yes", "no"}]
    tier, reason = cluster_quality_tier(out, observations, config)
    flags = joined_flags(observations, out)
    out.update(
        {
            "total_supporting_reads": len(observations),
            "distinct_observation_count": len(observations),
            "raw_source_record_count": sum(integer(row.get("source_record_count"), 1) for row in observations),
            "raw_supporting_read_count": len(observations),
            "paired_end_fragment_deduplicated_observation_count": len(observations),
            "distinct_read_identifier_count": len({str(row.get("read_id", "")) for row in observations}),
            "distinct_alignment_pattern_count": len(
                {
                    value
                    for row in observations
                    for value in str(row.get("alignment_pattern_id", "")).split(";")
                    if value
                }
            ),
            "support_per_sample": ";".join(f"{sample}:{count}" for sample, count in sorted(per_sample.items())),
            "maximum_support_within_sample": max(per_sample.values(), default=0),
            "samples_with_replicated_support": sum(count >= 2 for count in per_sample.values()),
            "samples_with_signal": len(per_sample),
            "supporting_samples": ",".join(sorted(per_sample)),
            "evidence_sources": ";".join(sorted(sources)),
            "star_supporting_observations": star_count,
            "minimap2_supporting_observations": minimap_count,
            "both_caller_supporting_observations": both_count,
            "both_caller_supporting_observation_fraction": both_count / len(observations) if observations else "",
            "evidence_status": "star_and_minimap2" if star_count and minimap_count else "star_only" if star_count else "minimap2_only",
            "primary_chain_evidence_fraction": proportion(primary_values),
            "secondary_only_evidence_fraction": proportion(secondary_values),
            "median_min_anchor_length": median(anchors) if anchors else "",
            "worst_min_anchor_length": min(anchors) if anchors else "",
            "median_alignment_error_rate": median(errors) if errors else "",
            "worst_alignment_error_rate": max(errors) if errors else "",
            "median_query_union_coverage": median(coverages) if coverages else "",
            "adjacent_query_segment_fraction": proportion([row.get("query_segments_adjacent") == "yes" for row in observations]),
            "complex_alignment_chain_observation_count": sum(value > 2 for value in segment_counts),
            "maximum_alignment_segment_count": max(segment_counts, default=""),
            "multiple_hypothesis_observation_count": sum(integer(row.get("deletion_hypotheses_from_read"), 1) > 1 for row in observations),
            "distinct_exact_breakpoint_pair_count": len(modal_pairs),
            "modal_observed_breakpoint_pair": f"{modal_pair[0]}:{modal_pair[1]}" if modal_pairs else "",
            "modal_observed_breakpoint_pair_support": modal_pair_support,
            "median_left_breakpoint_deviation_bp": median(left_deviations) if left_deviations else "",
            "maximum_left_breakpoint_deviation_bp": max(left_deviations, default=""),
            "median_right_breakpoint_deviation_bp": median(right_deviations) if right_deviations else "",
            "maximum_right_breakpoint_deviation_bp": max(right_deviations, default=""),
            "both_rotation_supporting_observation_count": sum(both_rotation_values),
            "both_rotation_supporting_observation_fraction": proportion(both_rotation_values),
            "star_gene_pair_label": modal_value(star_gene_pairs),
            "star_gene_pair_labels": ";".join(sorted(set(star_gene_pairs))),
            "breakpoint_flanking_gene_pair_label": gene_pair_label(out),
            "gene_pair_label": modal_value(star_gene_pairs) or gene_pair_label(out),
            "microhomology_or_repeat_length": max(
                [integer(row.get("star_junction_repeat_length")) for row in observations],
                default=0,
            ),
            "nuclear_competition_status": ";".join(
                sorted({str(row.get("nuclear_competition_status", "")) for row in observations if row.get("nuclear_competition_status")})
            ),
            "quality_tier": tier,
            "quality_tier_reason": reason,
            "quality_flags": ";".join(flags),
        }
    )
    return out


def finalize(clusters: list[dict], observations: list[dict], config: dict) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    by_cluster: dict[str, list[dict]] = defaultdict(list)
    for row in observations:
        cluster_id = row.get("exact_deletion_id") or row.get("junction_id")
        if cluster_id:
            by_cluster[cluster_id].append(row)
    final_clusters = []
    cluster_by_id = {}
    for cluster in clusters:
        cluster_id = cluster.get("exact_deletion_id") or cluster.get("junction_id")
        summarized = summarize_cluster(cluster, by_cluster.get(cluster_id, []), config)
        final_clusters.append(summarized)
        cluster_by_id[cluster_id] = summarized
    final_observations = []
    id_map = []
    for row in observations:
        out = dict(row)
        cluster_id = out.get("exact_deletion_id") or out.get("junction_id")
        cluster = cluster_by_id.get(cluster_id, {})
        out["quality_tier"] = cluster.get("quality_tier", "review")
        out["quality_tier_reason"] = cluster.get("quality_tier_reason", "")
        out["cluster_quality_flags"] = cluster.get("quality_flags", "")
        out["breakpoint_flanking_gene_pair_label"] = cluster.get("breakpoint_flanking_gene_pair_label", "")
        out["gene_pair_label"] = out.get("star_gene_pair_label") or cluster.get("breakpoint_flanking_gene_pair_label", "")
        final_observations.append(out)
        id_map.append(
            {
                "sample": out.get("sample", ""),
                "read_id": out.get("read_id", ""),
                "physical_observation_id": out.get("physical_observation_id", out.get("read_id", "")),
                "exact_deletion_id": cluster_id,
                "junction_id": cluster_id,
                "quality_tier": out["quality_tier"],
                "evidence_sources": out.get("evidence_sources", ""),
            }
        )
    profiles = configured_profiles(config)
    membership = []
    for cluster in final_clusters:
        cluster_id = cluster.get("exact_deletion_id") or cluster.get("junction_id")
        for profile, tiers in profiles.items():
            included = cluster.get("quality_tier") in tiers
            membership.append(
                {
                    "exact_deletion_id": cluster_id,
                    "report_profile": profile,
                    "included": "yes" if included else "no",
                    "quality_tier": cluster.get("quality_tier", ""),
                    "distinct_observation_count": cluster.get("distinct_observation_count", 0),
                    "reason": "included_tier" if included else cluster.get("quality_tier_reason", "tier_not_in_profile"),
                }
            )
    return final_clusters, final_observations, id_map, membership


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--observations", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-observations", required=True)
    parser.add_argument("--out-id-map", required=True)
    parser.add_argument("--out-membership", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()
    config = read_yaml(args.config)
    clusters, observations, id_map, membership = finalize(
        read_tsv(args.clusters),
        read_tsv(args.observations),
        config,
    )
    write_tsv(args.out_clusters, clusters)
    write_tsv(args.out_observations, observations)
    write_tsv(args.out_id_map, id_map)
    write_tsv(args.out_membership, membership)
    cluster_tiers = Counter(row.get("quality_tier", "unknown") for row in clusters)
    observation_tiers = Counter(row.get("quality_tier", "unknown") for row in observations)
    summary = []
    for tier in ["strong", "supported", "review", "rejected", "unknown"]:
        summary.append({"level": "cluster", "quality_tier": tier, "count": cluster_tiers[tier]})
        summary.append({"level": "observation", "quality_tier": tier, "count": observation_tiers[tier]})
    write_tsv(args.out_summary, summary, ["level", "quality_tier", "count"])


if __name__ == "__main__":
    main()
