#!/usr/bin/env python3
from __future__ import annotations

import argparse

from common import read_tsv, read_yaml, write_tsv
from finalize_quality_evidence import configured_profiles


def filter_profile(
    clusters: list[dict],
    observations: list[dict],
    id_map: list[dict],
    profile: str,
    config: dict,
) -> tuple[list[dict], list[dict], list[dict]]:
    profiles = configured_profiles(config)
    if profile not in profiles:
        raise ValueError(f"Unknown quality report profile: {profile}")
    included_tiers = set(profiles[profile])
    kept_clusters = [row for row in clusters if row.get("quality_tier") in included_tiers]
    kept_ids = {row.get("exact_deletion_id") or row.get("junction_id") for row in kept_clusters}
    kept_observations = [
        row
        for row in observations
        if (row.get("exact_deletion_id") or row.get("junction_id")) in kept_ids
    ]
    kept_map = [
        row
        for row in id_map
        if (row.get("exact_deletion_id") or row.get("junction_id")) in kept_ids
    ]
    return kept_clusters, kept_observations, kept_map


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--observations", required=True)
    parser.add_argument("--id-map", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-clusters", required=True)
    parser.add_argument("--out-observations", required=True)
    parser.add_argument("--out-id-map", required=True)
    args = parser.parse_args()
    clusters, observations, id_map = filter_profile(
        read_tsv(args.clusters),
        read_tsv(args.observations),
        read_tsv(args.id_map),
        args.profile,
        read_yaml(args.config),
    )
    write_tsv(args.out_clusters, clusters)
    write_tsv(args.out_observations, observations)
    write_tsv(args.out_id_map, id_map)


if __name__ == "__main__":
    main()
