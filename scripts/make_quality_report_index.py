#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
from pathlib import Path

from common import ensure_parent, read_tsv, read_yaml
from finalize_quality_evidence import configured_profiles


def profile_rows(config: dict, membership: list[dict]) -> list[dict]:
    profiles = configured_profiles(config)
    rows = []
    for name, tiers in profiles.items():
        included = {
            row.get("exact_deletion_id", "")
            for row in membership
            if row.get("report_profile") == name and row.get("included") == "yes"
        }
        observations = sum(
            int(float(row.get("distinct_observation_count", 0) or 0))
            for row in membership
            if row.get("report_profile") == name and row.get("included") == "yes"
        )
        rows.append(
            {
                "profile": name,
                "included_tiers": ", ".join(tiers),
                "deletion_clusters": len(included),
                "distinct_observations": observations,
                "role": "Primary interpretation" if name == config.get("quality", {}).get("primary_report_profile", "standard") else "Sensitivity view",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--title", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--reports", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = read_yaml(args.config)
    rows = profile_rows(config, read_tsv(args.membership))
    report_by_profile = {Path(path).parents[1].name: Path(path) for path in args.reports}
    cards = []
    output_parent = Path(args.output).resolve().parent
    for row in rows:
        profile = row["profile"]
        report = report_by_profile.get(profile)
        if report is None:
            continue
        relative = Path("..") / "profiles" / profile / ".report" / "index.html"
        cards.append(
            f'<article><h2>{html.escape(profile.capitalize())}</h2>'
            f'<p><strong>{row["deletion_clusters"]}</strong> exact deletion clusters</p>'
            f'<p><strong>{row["distinct_observations"]}</strong> distinct observations</p>'
            f'<p>Included tiers: {html.escape(row["included_tiers"])}</p>'
            f'<p>{html.escape(row["role"])}</p>'
            f'<a href="{html.escape(str(relative))}">Open report</a></article>'
        )
    primary = html.escape(str(config.get("quality", {}).get("primary_report_profile", "standard")))
    document = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(args.title)} quality reports</title>
  <style>
    body {{ margin: 0; font-family: system-ui, sans-serif; color: #1f2933; background: #f4f6f8; }}
    header {{ background: #243447; color: white; padding: 32px 40px; }}
    main {{ padding: 28px 40px; max-width: 1180px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 16px; }}
    article {{ background: white; border: 1px solid #d8dee8; border-radius: 7px; padding: 20px; }}
    article h2 {{ margin-top: 0; }}
    a {{ display: inline-block; background: #285f8f; color: white; padding: 8px 12px; border-radius: 5px; text-decoration: none; }}
  </style>
</head>
<body>
  <header><h1>{html.escape(args.title)}</h1><p>Evidence-quality report views generated from one shared canonical call set.</p></header>
  <main>
    <p>The predefined primary profile is <strong>{primary}</strong>. Each report rebuilds its matrices, PCA, plots, and summaries from its retained evidence. PCA axes are not assumed equivalent across profiles.</p>
    <div class="grid">{''.join(cards)}</div>
  </main>
</body>
</html>
"""
    ensure_parent(args.output)
    Path(args.output).write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
