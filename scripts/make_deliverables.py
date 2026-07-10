#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

from common import deep_update, ensure_parent, read_yaml, write_tsv, write_yaml


PLOTS = [
    "deletion_burden_by_sample.pdf",
    "unique_exact_deletions_by_sample.pdf",
    "deletion_burden_factorial_interaction.pdf",
    "unique_exact_deletions_factorial_interaction.pdf",
    "deletion_size_distribution_unweighted.pdf",
    "deletion_size_distribution_support_weighted.pdf",
    "deletion_size_distribution_support_weighted_log_y.pdf",
    "deletion_size_distribution_small.pdf",
    "deletion_size_distribution_medium.pdf",
    "deletion_size_distribution_large.pdf",
    "deletion_rainfall_left_breakpoint.pdf",
    "deletion_rainfall_right_breakpoint.pdf",
    "deletion_rainfall_midpoint.pdf",
    "breakpoint_pair_support_map.pdf",
    "pooled_breakpoint_support_density.pdf",
    "pooled_breakpoint_support_density_capped.pdf",
    "affected_feature_support.pdf",
    "affected_feature_counts.pdf",
    "affected_feature_proportions.pdf",
    "feature_impact_classes.pdf",
    "per_gene_affected_burden.pdf",
    "exact_deletion_recurrence.pdf",
    "exact_deletion_pca.pdf",
    "exact_deletion_bray_curtis_mds.pdf",
    "affected_feature_pca.pdf",
    "affected_feature_bray_curtis_mds.pdf",
]


TABLES = [
    ("junctions/junction_clusters.tsv", "tables/exact_deletions.tsv"),
    ("junctions/ambiguous_direction_reads.tsv", "tables/ambiguous_direction_reads.tsv"),
    ("analysis/deletion_burden.tsv", "tables/deletion_burden.tsv"),
    ("analysis/exact_deletion_comparison.tsv", "tables/exact_deletion_comparison.tsv"),
    ("analysis/affected_feature_comparison.tsv", "tables/affected_feature_comparison.tsv"),
    ("analysis/feature_impact_class_comparison.tsv", "tables/feature_impact_class_comparison.tsv"),
    ("analysis/deletion_size_distribution_tests.tsv", "tables/deletion_size_distribution_tests.tsv"),
    ("analysis/deletion_size_bin_summary.tsv", "tables/deletion_size_bin_summary.tsv"),
    ("analysis/breakpoint_reference_support.tsv", "tables/breakpoint_reference_support.tsv"),
    ("analysis/factorial_model_summary.tsv", "tables/factorial_model_summary.tsv"),
    ("analysis/deletion_metadata_associations.tsv", "tables/deletion_metadata_associations.tsv"),
    ("analysis/per_gene_affected_burden.tsv", "tables/per_gene_affected_burden.tsv"),
    ("analysis/qc_summary.tsv", "tables/qc_summary.tsv"),
    ("analysis/known_sequence_search_summary.tsv", "tables/known_sequence_search_summary.tsv"),
    ("analysis/known_sequence_search_hits.tsv", "tables/known_sequence_search_hits.tsv"),
]


COLUMN_DEFINITIONS = {
    "exact_deletion_id": "Stable identifier for one directed coordinate-level inferred deletion model.",
    "breakpoint_pair_id": "Unordered breakpoint-pair identifier for diagnostics only; reciprocal directions can have the same pair ID.",
    "left_breakpoint": "Retained flanking base before the alignment-directed deleted circular arc.",
    "right_breakpoint": "Retained flanking base after the alignment-directed deleted circular arc.",
    "deleted_size": "Number of reference bases in the directed circular interval, excluding both retained breakpoint bases.",
    "deleted_interval": "One or two 1-based closed reference intervals containing the inferred deleted bases.",
    "wraps_origin": "Whether the alignment-directed deleted interval crosses the configured coordinate origin.",
    "complement_deleted_size": "Size of the reciprocal circular arc, retained as a diagnostic and not used to assign the deletion.",
    "arc_assignment_method": "Method used to select the deleted circular arc; primary corrected results use alignment_directed.",
    "direction_status": "Whether directed evidence is accepted or conflicts with a reciprocal alignment hypothesis.",
    "rotation_agreement": "Whether evidence for the exact deletion was recorded from one or multiple reference rotations.",
    "affected_feature_label": "Deterministic genomic-order label of reference features overlapped by the directed deleted interval.",
    "total_supporting_reads": "Distinct sample/read evidence count after deduplication across reference rotations.",
    "local_split_support_fraction": "Split support divided by split support plus minimum local reference-spanning support; not automatically heteroplasmy.",
    "normalization_denominator": "Configured population of reads used for per-million normalization.",
    "normalization_reads": "Per-sample count used as the normalized-support denominator.",
}


def flatten_config(value: object, prefix: str = "") -> list[dict[str, str]]:
    rows = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(flatten_config(item, child))
    else:
        rendered = json.dumps(value, sort_keys=True) if isinstance(value, (list, bool)) or value is None else str(value)
        rows.append({"setting": prefix, "value": rendered})
    return rows


def column_category(column: str) -> str:
    if column in {"sample", "dataset", "species", "condition", "age", "treatment", "tissue"}:
        return "biological_metadata"
    if "breakpoint" in column or column in {"deleted_size", "deleted_interval", "wraps_origin"}:
        return "coordinates"
    if "support" in column or column.endswith("_reads") or column.endswith("_count"):
        return "evidence_or_count"
    if "normalization" in column or "per_million" in column:
        return "normalization"
    if column in {"affected_feature_label", "affected_features", "fully_removed_features", "partially_overlapped_features", "feature_impact_class"}:
        return "derived_annotation"
    return "workflow_or_analysis_field"


def data_dictionary_rows(table_dir: Path) -> list[dict[str, str]]:
    rows = []
    for path in sorted(table_dir.glob("*.tsv")):
        if path.name in {"data_dictionary.tsv", "run_methods.tsv"}:
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
        for column in header:
            rows.append(
                {
                    "table": path.name,
                    "column": column,
                    "category": column_category(column),
                    "definition": COLUMN_DEFINITIONS.get(column, "Workflow-generated field; see the report method and table guide for context."),
                }
            )
    return rows


MATRICES = [
    "exact_deletion_raw_counts.tsv",
    "exact_deletion_support_per_million_mt_reads.tsv",
    "affected_feature_raw_counts.tsv",
    "affected_feature_support_per_million_mt_reads.tsv",
    "feature_impact_class_raw_counts.tsv",
    "feature_impact_class_support_per_million_mt_reads.tsv",
]


PLOT_EXPLANATIONS = {
    "deletion_burden_by_sample.pdf": "Total post-remap deletion burden by sample and group, using the configured per-million normalization denominator.",
    "unique_exact_deletions_by_sample.pdf": "Number of distinct alignment-directed exact deletions per sample.",
    "deletion_burden_factorial_interaction.pdf": "Factorial age-by-treatment interaction view for normalized deletion burden when age and treatment metadata are present.",
    "unique_exact_deletions_factorial_interaction.pdf": "Factorial age-by-treatment interaction view for distinct exact deletions when age and treatment metadata are present.",
    "deletion_size_distribution_unweighted.pdf": "Deletion size distribution where each supporting read contributes one count.",
    "deletion_size_distribution_support_weighted.pdf": "Deletion size distribution weighted by normalized deletion support.",
    "deletion_size_distribution_support_weighted_log_y.pdf": "Support-weighted deletion size distribution with a log y-axis for low-abundance larger peaks.",
    "deletion_size_distribution_small.pdf": "Support-weighted deletion size distribution below 1 kb.",
    "deletion_size_distribution_medium.pdf": "Support-weighted deletion size distribution from 1 kb to 5 kb.",
    "deletion_size_distribution_large.pdf": "Support-weighted deletion size distribution at 5 kb and above.",
    "deletion_rainfall_left_breakpoint.pdf": "Location-size plot of post-remap deletion calls placed by alignment-directed left breakpoint.",
    "deletion_rainfall_right_breakpoint.pdf": "Location-size plot of post-remap deletion calls placed by alignment-directed right breakpoint.",
    "deletion_rainfall_midpoint.pdf": "Location-size plot of post-remap deletion calls placed by circular deleted-interval midpoint.",
    "breakpoint_pair_support_map.pdf": "Breakpoint-pair support map showing which deletion starts pair with which deletion ends.",
    "pooled_breakpoint_support_density.pdf": "Group-split pooled breakpoint endpoint support density using stacked binned left/right endpoint bars and a circular-smoothed pooled support curve.",
    "pooled_breakpoint_support_density_capped.pdf": "Group-split pooled breakpoint endpoint support density with a capped y-axis to reveal secondary endpoint hotspots.",
    "affected_feature_support.pdf": "Affected-feature categories compared by normalized support.",
    "affected_feature_counts.pdf": "Affected-feature categories compared by read support counts.",
    "affected_feature_proportions.pdf": "Within-group composition of affected-feature categories.",
    "feature_impact_classes.pdf": "Collapsed feature-impact classes for stable high-level interpretation.",
    "per_gene_affected_burden.pdf": "Per-gene affected burden; a deletion contributes to each feature it overlaps.",
    "exact_deletion_recurrence.pdf": "Top alignment-directed exact deletions ranked by supporting reads.",
    "exact_deletion_pca.pdf": "PCA using normalized exact-deletion support.",
    "exact_deletion_bray_curtis_mds.pdf": "Bray-Curtis MDS using normalized exact-deletion support.",
    "affected_feature_pca.pdf": "PCA using normalized affected-feature-category support.",
    "affected_feature_bray_curtis_mds.pdf": "Bray-Curtis MDS using normalized affected-feature-category support.",
}


TABLE_EXPLANATIONS = {
    "exact_deletions.tsv": "Alignment-directed post-remap exact deletion calls with coordinates, size, wrapping status, direction and rotation provenance, affected features, and configured deletion-target annotation.",
    "ambiguous_direction_reads.tsv": "Read-level reciprocal-direction conflicts retained for audit and excluded from primary summaries by default.",
    "run_methods.tsv": "Machine-readable resolved workflow configuration used to describe the run.",
    "data_dictionary.tsv": "Table and column inventory with field categories and definitions.",
    "deletion_burden.tsv": "Per-sample deletion burden, unique exact deletion count, and the configured normalization denominator.",
    "exact_deletion_comparison.tsv": "Group comparisons for recurrent exact deletions.",
    "affected_feature_comparison.tsv": "Group comparisons after exact deletions are summed by deterministic affected-feature label.",
    "feature_impact_class_comparison.tsv": "Group comparisons for broad collapsed feature-impact classes.",
    "deletion_size_distribution_tests.tsv": "Distribution-level tests comparing deletion sizes between groups.",
    "deletion_size_bin_summary.tsv": "Group summaries for deletion size bins, including small, medium, and large deletion support fractions.",
    "breakpoint_reference_support.tsv": "Local reference-spanning read counts at each deletion breakpoint and the resulting local split-support fraction.",
    "factorial_model_summary.tsv": "Age, treatment, and age-by-treatment model terms for sample-level deletion outcomes when a two-factor design is available.",
    "deletion_metadata_associations.tsv": "Associations between deletion outcomes and available metadata, including continuous variables when present.",
    "per_gene_affected_burden.tsv": "Per-gene affected burden; each deletion contributes to every feature overlapped by its deleted interval.",
    "qc_summary.tsv": "Brief processing QC and denominators for post-remap deletion analysis.",
    "known_sequence_search_summary.tsv": "Supplementary configured sequence-search counts for named breakpoint motifs.",
    "known_sequence_search_hits.tsv": "Read-level hits from configured sequence searches.",
}


def copy_if_exists(src: Path, dst: Path, copied: list[str]) -> None:
    if not src.exists() or not src.is_file():
        return
    ensure_parent(dst)
    shutil.copy2(src, dst)
    copied.append(str(dst))


def write_guide(path: Path, title: str, explanations: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [title, ""]
    for item in sorted(path.parent.iterdir()):
        if item.name == path.name or item.suffix not in {".pdf", ".tsv"}:
            continue
        lines.extend([item.name, f"  {explanations.get(item.name, 'Workflow-generated file.')}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--defaults", default="config/defaults.yaml")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--complete", required=True)
    args = parser.parse_args()

    root = Path(args.results_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    stale_star = out / "star_full_genome_split_read"
    if stale_star.exists():
        shutil.rmtree(stale_star)
    for subdir in ("plots", "tables", "matrices", "config"):
        path = out / subdir
        if path.exists():
            shutil.rmtree(path)
    copied: list[str] = []

    copy_if_exists(root / ".report" / "index.html", out / "index.html", copied)
    report_read_lists = root / ".report" / "read_lists"
    deliverable_read_lists = out / "read_lists"
    if report_read_lists.exists():
        if deliverable_read_lists.exists():
            shutil.rmtree(deliverable_read_lists)
        shutil.copytree(report_read_lists, deliverable_read_lists)
        copied.append(str(deliverable_read_lists))
    for plot in PLOTS:
        src = root / "plots" / plot
        copy_if_exists(src, out / "plots" / plot, copied)
        copy_if_exists(src.with_suffix(".svg"), out / "plots" / src.with_suffix(".svg").name, copied)
        if plot in {
            "deletion_rainfall_left_breakpoint.pdf",
            "deletion_rainfall_right_breakpoint.pdf",
            "deletion_rainfall_midpoint.pdf",
            "breakpoint_pair_support_map.pdf",
            "pooled_breakpoint_support_density.pdf",
            "pooled_breakpoint_support_density_capped.pdf",
        }:
            for sidecar in sorted((root / "plots").glob(f"{Path(plot).stem}__*.pdf")):
                copy_if_exists(sidecar, out / "plots" / sidecar.name, copied)
                copy_if_exists(sidecar.with_suffix(".svg"), out / "plots" / sidecar.with_suffix(".svg").name, copied)
    for src_rel, dst_rel in TABLES:
        copy_if_exists(root / src_rel, out / dst_rel, copied)
    for matrix in MATRICES:
        copy_if_exists(root / "matrices" / matrix, out / "matrices" / matrix, copied)
    merged_config = deep_update(read_yaml(args.defaults), read_yaml(args.config))
    write_yaml(out / "config" / "resolved_config.yaml", merged_config)
    copied.append(str(out / "config" / "resolved_config.yaml"))

    write_tsv(out / "tables" / "run_methods.tsv", flatten_config(merged_config), ["setting", "value"])
    write_tsv(out / "tables" / "data_dictionary.tsv", data_dictionary_rows(out / "tables"), ["table", "column", "category", "definition"])
    copied.extend([str(out / "tables" / "run_methods.tsv"), str(out / "tables" / "data_dictionary.tsv")])

    write_guide(out / "plots" / "README.txt", "Plot guide", PLOT_EXPLANATIONS)
    write_guide(out / "tables" / "README.txt", "Analysis table guide", TABLE_EXPLANATIONS)

    readme = out / "README.txt"
    readme.write_text(
        "\n".join(
            [
                f"Deliverables for dataset: {args.dataset}",
                "",
                "This folder excludes large FASTQ, BAM, reference index, and intermediate alignment files.",
                "",
                "Start here:",
                "- index.html",
                "",
                "Main analysis files:",
                "- tables/exact_deletions.tsv",
                "- tables/deletion_burden.tsv",
                "- tables/exact_deletion_comparison.tsv",
                "- tables/affected_feature_comparison.tsv",
                "- tables/deletion_metadata_associations.tsv",
                "- tables/feature_impact_class_comparison.tsv",
                "- matrices/exact_deletion_support_per_million_mt_reads.tsv",
                "- matrices/affected_feature_support_per_million_mt_reads.tsv",
                "- plots/*.pdf and plots/*.svg",
                "- read_lists/*.read_names.tsv and configured_vs_remap__*.tsv",
                "",
                "Guides:",
                "- plots/README.txt",
                "- tables/README.txt",
                "",
                "Copied files:",
                *copied,
                "",
            ]
        ),
        encoding="utf-8",
    )
    ensure_parent(args.complete)
    Path(args.complete).write_text("complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
