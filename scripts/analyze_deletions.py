#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact, ks_2samp, linregress, mannwhitneyu, spearmanr, t

from common import ensure_parent


def bh(p_values: list[float]) -> list[float]:
    indexed = [(i, p) for i, p in enumerate(p_values) if pd.notna(p)]
    adjusted = [float("nan")] * len(p_values)
    if not indexed:
        return adjusted
    indexed.sort(key=lambda item: item[1])
    n = len(indexed)
    prev = 1.0
    for rank, (idx, p) in reversed(list(enumerate(indexed, start=1))):
        q = min(prev, p * n / rank)
        adjusted[idx] = q
        prev = q
    return adjusted


def add_fdr(rows: list[dict], p_col: str = "p_value", q_col: str = "q_value_bh") -> list[dict]:
    q_values = bh([row.get(p_col, float("nan")) for row in rows])
    for row, q in zip(rows, q_values):
        row[q_col] = q
    return rows


def read_fragment_count(path: str) -> dict:
    if not Path(path).exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    if "metric" in df.columns and "value" in df.columns:
        return dict(zip(df["metric"], df["value"]))
    if "fragments" in df.columns and len(df):
        return {"total_fragments": df["fragments"].iloc[0]}
    return {}


def read_json_count(path: str) -> dict:
    if not Path(path).exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def normalize_matrix(raw: pd.DataFrame, samples: pd.DataFrame, denominator: pd.Series, scale: float = 1_000_000.0) -> pd.DataFrame:
    meta = samples.copy()
    feature_cols = [col for col in raw.columns if col not in samples.columns]
    out = raw[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    denom = denominator.replace(0, np.nan).astype(float)
    out = out.div(denom, axis=0).fillna(0.0) * scale
    return pd.concat([meta.reset_index(drop=True), out.reset_index(drop=True)], axis=1)


def count_matrix(samples: pd.DataFrame, rows: pd.DataFrame, column: str) -> pd.DataFrame:
    meta = samples.copy()
    labels = sorted(str(value) for value in rows[column].dropna().unique()) if not rows.empty and column in rows.columns else []
    counts = pd.DataFrame(0, index=meta.index, columns=labels)
    if labels:
        grouped = rows.groupby(["sample", column]).size()
        sample_index = {sample: i for i, sample in enumerate(meta["sample"])}
        for (sample, label), value in grouped.items():
            if sample in sample_index and str(label) in counts.columns:
                counts.loc[sample_index[sample], str(label)] = int(value)
    return pd.concat([meta.reset_index(drop=True), counts.reset_index(drop=True)], axis=1)


def deduplicate_evidence_reads(reads: pd.DataFrame) -> pd.DataFrame:
    if reads.empty:
        return reads
    key = [col for col in ["sample", "read_id", "junction_id"] if col in reads.columns]
    if len(key) < 3:
        return reads
    work = reads.copy()
    if {"left_anchor_length", "right_anchor_length"}.issubset(work.columns):
        left = pd.to_numeric(work["left_anchor_length"], errors="coerce").fillna(0)
        right = pd.to_numeric(work["right_anchor_length"], errors="coerce").fillna(0)
        work["_anchor_support_for_dedup"] = left + right
        work = work.sort_values("_anchor_support_for_dedup", ascending=False)
        work = work.drop_duplicates(key, keep="first").drop(columns=["_anchor_support_for_dedup"])
    else:
        work = work.drop_duplicates(key, keep="first")
    return work.reset_index(drop=True)


def expected_transcript_mask(reads: pd.DataFrame) -> pd.Series:
    if reads.empty or "junction_interpretation" not in reads.columns:
        return pd.Series(False, index=reads.index)
    return reads["junction_interpretation"].astype(str).eq("expected_transcript_junction")


def abundance_change_label(left_mean: float, right_mean: float) -> str:
    if left_mean <= 0 and right_mean <= 0:
        return "absent_in_both_groups"
    if left_mean <= 0:
        return "present_only_in_right_group"
    if right_mean <= 0:
        return "present_only_in_left_group"
    if right_mean > left_mean:
        return "higher_in_right_group"
    if right_mean < left_mean:
        return "higher_in_left_group"
    return "same_group_mean"


def finite_fold_change(left_mean: float, right_mean: float) -> float:
    if left_mean <= 0 or right_mean <= 0:
        return float("nan")
    return float(right_mean / left_mean)


def comparison_rows(
    samples: pd.DataFrame,
    raw: pd.DataFrame,
    normalized: pd.DataFrame,
    group_col: str,
    feature_cols: list[str],
    label_col: str,
    log2fc_pseudocount: float,
    denominator_col: str = "normalization_reads",
) -> list[dict]:
    rows = []
    if not group_col or group_col not in samples.columns or len(feature_cols) == 0:
        return rows
    groups = [value for value in samples[group_col].dropna().unique()]
    for left, right in combinations(groups, 2):
        left_mask = samples[group_col] == left
        right_mask = samples[group_col] == right
        for feature in feature_cols:
            left_values = pd.to_numeric(normalized.loc[left_mask, feature], errors="coerce").fillna(0)
            right_values = pd.to_numeric(normalized.loc[right_mask, feature], errors="coerce").fillna(0)
            left_raw = pd.to_numeric(raw.loc[left_mask, feature], errors="coerce").fillna(0)
            right_raw = pd.to_numeric(raw.loc[right_mask, feature], errors="coerce").fillna(0)
            try:
                p_value = mannwhitneyu(left_values, right_values, alternative="two-sided").pvalue
                test = "mann_whitney_u_normalized_support"
            except ValueError:
                p_value = float("nan")
                test = "not_tested"
            try:
                fisher_p = fisher_exact([[int(left_raw.sum()), int(right_raw.sum())], [int((left_raw == 0).sum()), int((right_raw == 0).sum())]])[1]
            except ValueError:
                fisher_p = float("nan")
            try:
                depth_col = denominator_col if denominator_col in samples.columns else "reads_passed_to_minimap2"
                left_depth = pd.to_numeric(samples.loc[left_mask, depth_col], errors="coerce").fillna(0).sum()
                right_depth = pd.to_numeric(samples.loc[right_mask, depth_col], errors="coerce").fillna(0).sum()
                left_support = int(left_raw.sum())
                right_support = int(right_raw.sum())
                read_depth_fisher_p = fisher_exact(
                    [
                        [left_support, max(int(left_depth) - left_support, 0)],
                        [right_support, max(int(right_depth) - right_support, 0)],
                    ]
                )[1]
            except (KeyError, ValueError):
                read_depth_fisher_p = float("nan")
            left_mean = float(left_values.mean())
            right_mean = float(right_values.mean())
            rows.append(
                {
                    label_col: feature,
                    "group_column": group_col,
                    "left_group": left,
                    "right_group": right,
                    "left_mean_per_million_mt_reads": left_mean,
                    "right_mean_per_million_mt_reads": right_mean,
                    "difference_per_million_mt_reads": float(right_mean - left_mean),
                    "abundance_change": abundance_change_label(left_mean, right_mean),
                    "log2_fold_change_right_over_left": float(np.log2((right_mean + log2fc_pseudocount) / (left_mean + log2fc_pseudocount))),
                    "log2_fold_change_pseudocount_per_million": log2fc_pseudocount,
                    "fold_change_right_over_left": finite_fold_change(left_mean, right_mean),
                    "left_total_supporting_reads": int(left_raw.sum()),
                    "right_total_supporting_reads": int(right_raw.sum()),
                    "samples_with_signal": int(((left_raw.append(right_raw) if hasattr(left_raw, "append") else pd.concat([left_raw, right_raw])) > 0).sum()),
                    "test": test,
                    "p_value": p_value,
                    "fisher_presence_p": fisher_p,
                    "read_depth_fisher_p": read_depth_fisher_p,
                }
            )
    rows = add_fdr(rows)
    return add_fdr(rows, p_col="read_depth_fisher_p", q_col="read_depth_fisher_q_value_bh")


def numeric_metadata(series: pd.Series) -> pd.Series | None:
    direct = pd.to_numeric(series, errors="coerce")
    if direct.notna().sum() >= max(3, len(series) // 2):
        return direct
    extracted = series.astype(str).str.extract(r"([-+]?\d*\.?\d+)")[0]
    numeric = pd.to_numeric(extracted, errors="coerce")
    if numeric.notna().sum() >= max(3, len(series) // 2):
        return numeric
    return None


def burden_associations(samples: pd.DataFrame, burden: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    metrics = ["deletion_support_per_million_mt_reads", "unique_exact_deletions"]
    for col in group_cols:
        if col not in samples.columns or samples[col].nunique(dropna=True) < 2:
            continue
        numeric = numeric_metadata(samples[col])
        if numeric is not None:
            mask = numeric.notna()
            for metric in metrics:
                if mask.sum() >= 3:
                    y = pd.to_numeric(burden.loc[mask, metric], errors="coerce").fillna(0)
                    try:
                        sp = spearmanr(numeric[mask], y)
                        lr = linregress(numeric[mask], y)
                        rows.append(
                            {
                                "metadata_column": col,
                                "metadata_type": "numeric_or_ordered",
                                "outcome": metric,
                                "test": "spearman_and_linear_regression",
                                "spearman_rho": sp.statistic,
                                "linear_slope": lr.slope,
                                "p_value": sp.pvalue,
                                "linear_p": lr.pvalue,
                            }
                        )
                    except ValueError:
                        pass
        elif samples[col].nunique(dropna=True) <= 12:
            groups = sorted(samples[col].dropna().unique())
            if len(groups) == 2:
                left, right = groups
                for metric in metrics:
                    lv = pd.to_numeric(burden.loc[samples[col] == left, metric], errors="coerce").fillna(0)
                    rv = pd.to_numeric(burden.loc[samples[col] == right, metric], errors="coerce").fillna(0)
                    try:
                        p = mannwhitneyu(lv, rv, alternative="two-sided").pvalue
                    except ValueError:
                        p = float("nan")
                    rows.append(
                        {
                            "metadata_column": col,
                            "metadata_type": "categorical",
                            "outcome": metric,
                            "test": "mann_whitney_u",
                            "left_group": left,
                            "right_group": right,
                            "left_mean": float(lv.mean()),
                            "right_mean": float(rv.mean()),
                            "p_value": p,
                        }
                    )
    return pd.DataFrame(add_fdr(rows))


def size_distribution_tests(reads: pd.DataFrame, samples: pd.DataFrame, group_col: str) -> pd.DataFrame:
    rows = []
    if reads.empty or not group_col or group_col not in samples.columns:
        return pd.DataFrame()
    merged = reads.merge(samples[["sample", group_col]], on="sample", how="left")
    groups = [value for value in merged[group_col].dropna().unique()]
    for left, right in combinations(groups, 2):
        left_sizes = pd.to_numeric(merged.loc[merged[group_col] == left, "deleted_size"], errors="coerce").dropna()
        right_sizes = pd.to_numeric(merged.loc[merged[group_col] == right, "deleted_size"], errors="coerce").dropna()
        if len(left_sizes) and len(right_sizes):
            rows.append(
                {
                    "group_column": group_col,
                    "left_group": left,
                    "right_group": right,
                    "test": "kolmogorov_smirnov_unweighted_deletion_sizes",
                    "left_n_deletion_supporting_reads": len(left_sizes),
                    "right_n_deletion_supporting_reads": len(right_sizes),
                    "left_median_deleted_size": float(left_sizes.median()),
                    "right_median_deleted_size": float(right_sizes.median()),
                    "p_value": ks_2samp(left_sizes, right_sizes).pvalue,
                }
            )
            try:
                rows.append(
                    {
                        "group_column": group_col,
                        "left_group": left,
                        "right_group": right,
                        "test": "mann_whitney_u_unweighted_deletion_sizes",
                        "left_n_deletion_supporting_reads": len(left_sizes),
                        "right_n_deletion_supporting_reads": len(right_sizes),
                        "left_median_deleted_size": float(left_sizes.median()),
                        "right_median_deleted_size": float(right_sizes.median()),
                        "p_value": mannwhitneyu(left_sizes, right_sizes, alternative="two-sided").pvalue,
                    }
                )
            except ValueError:
                pass
    return pd.DataFrame(add_fdr(rows))


def sorted_factor_levels(series: pd.Series) -> list[str]:
    values = series.dropna().astype(str).unique().tolist()

    def key(value: str) -> tuple[int, float, str]:
        lower = value.lower()
        number = pd.to_numeric(pd.Series([lower]).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce").iloc[0]
        if "control" in lower or lower in {"ctrl", "vehicle", "untreated"}:
            rank = 0
        else:
            rank = 1
        return (rank, float(number) if pd.notna(number) else 1e9, lower)

    return sorted(values, key=key)


def age_levels(series: pd.Series) -> list[str]:
    values = series.dropna().astype(str).unique().tolist()

    def key(value: str) -> tuple[float, str]:
        number = pd.to_numeric(pd.Series([value]).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce").iloc[0]
        return (float(number) if pd.notna(number) else 1e9, value.lower())

    return sorted(values, key=key)


def sample_size_metrics(annotated_reads: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    denom_col = "normalization_reads" if "normalization_reads" in samples.columns else "reads_passed_to_minimap2"
    denom = samples.set_index("sample")[denom_col].to_dict() if denom_col in samples.columns else {}
    if annotated_reads.empty:
        return pd.DataFrame()
    work = annotated_reads.copy()
    work["deleted_size"] = pd.to_numeric(work["deleted_size"], errors="coerce")
    for sample, sub in work.groupby("sample"):
        sample_denom = float(denom.get(sample, 0) or 0)
        weight = 1_000_000 / sample_denom if sample_denom > 0 else 0.0
        sizes = sub["deleted_size"].dropna()
        rows.append(
            {
                "sample": sample,
                "median_deletion_size": float(sizes.median()) if len(sizes) else float("nan"),
                "large_deletion_support_per_million_mt_reads": float((sizes >= 5000).sum() * weight),
                "small_lt_1kb_support_per_million_mt_reads": float((sizes < 1000).sum() * weight),
                "medium_1kb_to_4999bp_support_per_million_mt_reads": float(((sizes >= 1000) & (sizes < 5000)).sum() * weight),
            }
        )
    return pd.DataFrame(rows)


def size_bin_summary(annotated_reads: pd.DataFrame, samples: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if annotated_reads.empty:
        return pd.DataFrame()
    work = annotated_reads[["sample", "deleted_size"]].copy()
    work["deleted_size"] = pd.to_numeric(work["deleted_size"], errors="coerce")
    work = work.dropna(subset=["deleted_size"])
    denom_col = "normalization_reads" if "normalization_reads" in samples.columns else "reads_passed_to_minimap2"
    meta_cols = ["sample", group_col, denom_col] if group_col in samples.columns else ["sample", denom_col]
    work = work.merge(samples[[col for col in meta_cols if col in samples.columns]], on="sample", how="left")
    denom = pd.to_numeric(work.get(denom_col, pd.Series(0, index=work.index)), errors="coerce")
    work["support_per_million_mt_reads"] = np.where(denom > 0, 1_000_000 / denom, 0)
    bins = [
        ("small_lt_1kb", work["deleted_size"] < 1000),
        ("medium_1kb_to_4999bp", (work["deleted_size"] >= 1000) & (work["deleted_size"] < 5000)),
        ("large_ge_5kb", work["deleted_size"] >= 5000),
    ]
    group_key = group_col if group_col in work.columns else None
    rows = []
    groups = work[group_key].dropna().astype(str).unique().tolist() if group_key else ["all"]
    for group in groups:
        sub = work[work[group_key].astype(str) == group] if group_key else work
        total = float(sub["support_per_million_mt_reads"].sum())
        row = {
            "group": group,
            "n_deletion_supporting_reads": int(len(sub)),
            "median_deletion_size": float(sub["deleted_size"].median()) if len(sub) else float("nan"),
            "total_support_per_million_mt_reads": total,
        }
        for label, mask in bins:
            value = float(sub.loc[mask.loc[sub.index], "support_per_million_mt_reads"].sum())
            row[f"{label}_support_per_million_mt_reads"] = value
            row[f"{label}_fraction_of_support"] = value / total if total > 0 else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def factorial_model_summary(samples: pd.DataFrame, burden: pd.DataFrame, annotated_reads: pd.DataFrame) -> pd.DataFrame:
    if not {"age", "treatment"}.issubset(samples.columns):
        return pd.DataFrame()
    ages = age_levels(samples["age"])
    treatments = sorted_factor_levels(samples["treatment"])
    if len(ages) != 2 or len(treatments) != 2:
        return pd.DataFrame()
    metrics = burden[["sample", "deletion_support_per_million_mt_reads", "unique_exact_deletions"]].copy()
    size_metrics = sample_size_metrics(annotated_reads, samples)
    if not size_metrics.empty:
        metrics = metrics.merge(size_metrics, on="sample", how="left")
    data = samples[["sample", "age", "treatment"]].merge(metrics, on="sample", how="left")
    data["_age_high"] = data["age"].astype(str).eq(ages[1]).astype(float)
    data["_treatment_noncontrol"] = data["treatment"].astype(str).eq(treatments[1]).astype(float)
    data["_interaction"] = data["_age_high"] * data["_treatment_noncontrol"]
    terms = [
        ("age", "_age_high", f"{ages[1]} vs {ages[0]}"),
        ("treatment", "_treatment_noncontrol", f"{treatments[1]} vs {treatments[0]}"),
        ("age:treatment", "_interaction", f"extra {treatments[1]} effect in {ages[1]}"),
    ]
    outcomes = [
        "deletion_support_per_million_mt_reads",
        "unique_exact_deletions",
        "median_deletion_size",
        "large_deletion_support_per_million_mt_reads",
        "small_lt_1kb_support_per_million_mt_reads",
        "medium_1kb_to_4999bp_support_per_million_mt_reads",
    ]
    rows = []
    x_cols = ["_age_high", "_treatment_noncontrol", "_interaction"]
    for outcome in outcomes:
        if outcome not in data.columns:
            continue
        model_data = data[["sample", "age", "treatment", outcome, *x_cols]].dropna(subset=[outcome])
        if len(model_data) < 5:
            continue
        y = pd.to_numeric(model_data[outcome], errors="coerce").to_numpy(dtype=float)
        x = np.column_stack([np.ones(len(model_data)), model_data[x_cols].to_numpy(dtype=float)])
        try:
            beta, *_ = np.linalg.lstsq(x, y, rcond=None)
            resid = y - x @ beta
            dof = len(y) - x.shape[1]
            if dof <= 0:
                continue
            sigma2 = float((resid @ resid) / dof)
            cov = sigma2 * np.linalg.pinv(x.T @ x)
            se = np.sqrt(np.diag(cov))
            for idx, (term, col, contrast) in enumerate(terms, start=1):
                estimate = float(beta[idx])
                stderr = float(se[idx]) if se[idx] > 0 else float("nan")
                statistic = estimate / stderr if pd.notna(stderr) and stderr > 0 else float("nan")
                p_value = float(2 * t.sf(abs(statistic), dof)) if pd.notna(statistic) else float("nan")
                rows.append(
                    {
                        "outcome": outcome,
                        "term": term,
                        "contrast": contrast,
                        "estimate": estimate,
                        "direction": "positive" if estimate > 0 else ("negative" if estimate < 0 else "zero"),
                        "standard_error": stderr,
                        "t_statistic": statistic,
                        "degrees_of_freedom": dof,
                        "p_value": p_value,
                    }
                )
        except np.linalg.LinAlgError:
            continue
    return pd.DataFrame(add_fdr(rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--id-map", required=True)
    parser.add_argument("--all-reads", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--group-column", default="")
    parser.add_argument("--group-columns", default="")
    parser.add_argument("--normalization-denominator", choices=["mt_evidence_reads", "total_usable_reads"], default="total_usable_reads")
    parser.add_argument("--fragment-counts", nargs="*", default=[])
    parser.add_argument("--mt-summaries", nargs="*", default=[])
    parser.add_argument("--out-exact-raw", required=True)
    parser.add_argument("--out-exact-mtpm", required=True)
    parser.add_argument("--out-affected-raw", required=True)
    parser.add_argument("--out-affected-mtpm", required=True)
    parser.add_argument("--out-impact-class-raw", required=True)
    parser.add_argument("--out-impact-class-mtpm", required=True)
    parser.add_argument("--out-per-gene-burden", required=True)
    parser.add_argument("--out-burden", required=True)
    parser.add_argument("--out-exact-comparison", required=True)
    parser.add_argument("--out-affected-comparison", required=True)
    parser.add_argument("--out-impact-class-comparison", required=True)
    parser.add_argument("--out-size-tests", required=True)
    parser.add_argument("--out-size-bin-summary", required=True)
    parser.add_argument("--out-factorial-model-summary", required=True)
    parser.add_argument("--out-metadata-associations", required=True)
    parser.add_argument("--out-qc-summary", required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t")
    clusters = pd.read_csv(args.clusters, sep="\t")
    id_map = pd.read_csv(args.id_map, sep="\t")
    all_reads = pd.read_csv(args.all_reads, sep="\t")
    if "exact_deletion_id" not in clusters.columns and "junction_id" in clusters.columns:
        clusters["exact_deletion_id"] = clusters["junction_id"]
    if "exact_deletion_id" not in all_reads.columns and "junction_id" in all_reads.columns:
        all_reads["exact_deletion_id"] = all_reads["junction_id"]
    with open(args.config, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    log2fc_pseudocount = float(config.get("analysis", {}).get("effect_size_pseudocount_per_million", 0.5))

    fragment_totals = {}
    for path in args.fragment_counts:
        sample = Path(path).parts[-2] if len(Path(path).parts) >= 2 else Path(path).stem
        metrics = read_fragment_count(path)
        fragment_totals[sample] = float(metrics.get("post_filtering_reads", metrics.get("read1_after_filtering", metrics.get("total_fragments", 0))) or 0)
    mt_totals = {}
    for path in args.mt_summaries:
        sample = Path(path).name.replace(".mt_read_summary.json", "")
        metrics = read_json_count(path)
        mt_totals[sample] = float(metrics.get("mt_evidence_fastq_records_written", metrics.get("mt_evidence_reads_selected", 0)) or 0)

    samples = samples.copy()
    samples["total_usable_reads"] = samples["sample"].map(fragment_totals).fillna(0).astype(float)
    samples["reads_passed_to_minimap2"] = samples["sample"].map(mt_totals).fillna(0).astype(float)

    cluster_meta = clusters.drop(columns=[col for col in ["supporting_samples"] if col in clusters.columns])
    annotated_reads = all_reads.merge(cluster_meta, on="junction_id", how="inner", suffixes=("_read", ""))
    for col in ["exact_deletion_id", "left_breakpoint", "right_breakpoint", "deleted_size", "wraps_origin"]:
        cluster_col = col
        read_col = f"{col}_read"
        if cluster_col not in annotated_reads.columns and read_col in annotated_reads.columns:
            annotated_reads[cluster_col] = annotated_reads[read_col]
    annotated_reads["exact_deletion_id"] = annotated_reads["junction_id"]
    annotated_reads = deduplicate_evidence_reads(annotated_reads)
    annotated_reads_before_transcript_filter = annotated_reads.copy()
    transcript_mask = expected_transcript_mask(annotated_reads_before_transcript_filter)
    exclude_expected = bool(config.get("junctions", {}).get("exclude_expected_transcript_junctions", True))
    if exclude_expected and transcript_mask.any():
        annotated_reads = annotated_reads_before_transcript_filter.loc[~transcript_mask].copy()

    exact_raw = count_matrix(samples, annotated_reads, "exact_deletion_id")
    exact_cols = [col for col in exact_raw.columns if col not in samples.columns]
    denominator = samples["total_usable_reads"] if args.normalization_denominator == "total_usable_reads" else samples["reads_passed_to_minimap2"]
    samples["normalization_denominator"] = args.normalization_denominator
    samples["normalization_reads"] = denominator
    exact_mtpm = normalize_matrix(exact_raw, samples, denominator)

    affected_raw = count_matrix(samples, annotated_reads, "affected_feature_label")
    affected_cols = [col for col in affected_raw.columns if col not in samples.columns]
    affected_mtpm = normalize_matrix(affected_raw, samples, denominator)

    impact_raw = count_matrix(samples, annotated_reads, "feature_impact_class")
    impact_cols = [col for col in impact_raw.columns if col not in samples.columns]
    impact_mtpm = normalize_matrix(impact_raw, samples, denominator)

    burden = samples.copy()
    burden["deletion_supporting_reads"] = exact_raw[exact_cols].sum(axis=1) if exact_cols else 0
    burden["deletion_support_per_million_mt_reads"] = exact_mtpm[exact_cols].sum(axis=1) if exact_cols else 0.0
    burden["normalization_denominator"] = args.normalization_denominator
    burden["normalization_reads"] = denominator
    burden["unique_exact_deletions"] = (exact_raw[exact_cols] > 0).sum(axis=1) if exact_cols else 0
    burden["unique_affected_feature_categories"] = (affected_raw[affected_cols] > 0).sum(axis=1) if affected_cols else 0

    per_gene_rows = []
    if not annotated_reads.empty and "per_feature_overlap_details" in annotated_reads.columns:
        mt_denom_by_sample = samples.set_index("sample")["normalization_reads"].to_dict()
        for _, row in annotated_reads.iterrows():
            details = str(row.get("per_feature_overlap_details", ""))
            if not details or details.lower() == "nan":
                continue
            sample = row["sample"]
            denom = float(mt_denom_by_sample.get(sample, 0) or 0)
            for token in details.split(";"):
                if not token:
                    continue
                gene = token.split(":")[0]
                per_gene_rows.append(
                    {
                        "sample": sample,
                        "feature": gene,
                        "supporting_reads": 1,
                        "support_per_million_mt_reads": 1_000_000 / denom if denom > 0 else 0.0,
                    }
                )
    per_gene = pd.DataFrame(per_gene_rows)
    if not per_gene.empty:
        per_gene = per_gene.groupby(["sample", "feature"], as_index=False).sum().merge(samples, on="sample", how="left")

    group_cols = [col for col in args.group_columns.split(",") if col] or ([args.group_column] if args.group_column else [])
    exact_comp = comparison_rows(samples, exact_raw, exact_mtpm, args.group_column, exact_cols, "exact_deletion_id", log2fc_pseudocount)
    if exact_comp and not clusters.empty:
        meta_cols = [
            "exact_deletion_id",
            "left_breakpoint",
            "right_breakpoint",
            "deleted_size",
            "wraps_origin",
            "affected_feature_label",
            "feature_impact_class",
            "known_deletion_label",
        ]
        available = [col for col in meta_cols if col in clusters.columns]
        exact_comp = pd.DataFrame(exact_comp).merge(clusters[available], on="exact_deletion_id", how="left").to_dict(orient="records")
    affected_comp = comparison_rows(samples, affected_raw, affected_mtpm, args.group_column, affected_cols, "affected_feature_label", log2fc_pseudocount)
    impact_comp = comparison_rows(samples, impact_raw, impact_mtpm, args.group_column, impact_cols, "feature_impact_class", log2fc_pseudocount)
    size_tests = size_distribution_tests(annotated_reads, samples, args.group_column)
    size_bins = size_bin_summary(annotated_reads, samples, args.group_column)
    factorial_models = factorial_model_summary(samples, burden, annotated_reads)
    metadata_assoc = burden_associations(samples, burden, group_cols)

    pre_filter_counts = annotated_reads_before_transcript_filter.groupby("sample").size()
    expected_counts = annotated_reads_before_transcript_filter.loc[transcript_mask].groupby("sample").size()
    candidate_counts = annotated_reads_before_transcript_filter.loc[~transcript_mask].groupby("sample").size()

    qc = samples[
        [
            "sample",
            "total_usable_reads",
            "reads_passed_to_minimap2",
        ]
    ].copy()
    qc["clustered_split_read_alignments_before_transcript_filter"] = qc["sample"].map(pre_filter_counts).fillna(0).astype(int)
    qc["expected_transcript_compatible_split_reads"] = qc["sample"].map(expected_counts).fillna(0).astype(int)
    qc["candidate_deletion_split_reads_before_transcript_filter"] = qc["sample"].map(candidate_counts).fillna(0).astype(int)
    qc["expected_transcript_filter_applied"] = exclude_expected
    qc["deletion_supporting_reads"] = burden["deletion_supporting_reads"]
    qc["unique_exact_deletions"] = burden["unique_exact_deletions"]
    qc["warning_flags"] = np.where(qc["reads_passed_to_minimap2"] == 0, "no_mitochondrial_evidence_reads", "")

    outputs = {
        args.out_exact_raw: exact_raw,
        args.out_exact_mtpm: exact_mtpm,
        args.out_affected_raw: affected_raw,
        args.out_affected_mtpm: affected_mtpm,
        args.out_impact_class_raw: impact_raw,
        args.out_impact_class_mtpm: impact_mtpm,
        args.out_per_gene_burden: per_gene,
        args.out_burden: burden,
        args.out_exact_comparison: pd.DataFrame(exact_comp),
        args.out_affected_comparison: pd.DataFrame(affected_comp),
        args.out_impact_class_comparison: pd.DataFrame(impact_comp),
        args.out_size_tests: size_tests,
        args.out_size_bin_summary: size_bins,
        args.out_factorial_model_summary: factorial_models,
        args.out_metadata_associations: metadata_assoc,
        args.out_qc_summary: qc,
    }
    for path, df in outputs.items():
        ensure_parent(path)
        if df is None or df.empty:
            pd.DataFrame().to_csv(path, sep="\t", index=False)
        else:
            df.to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main()
