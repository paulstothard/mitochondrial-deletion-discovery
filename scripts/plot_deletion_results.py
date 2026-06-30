#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors, ticker
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from common import ensure_parent


sns.set_theme(style="whitegrid", context="notebook")

TECHNICAL_COLUMNS = {
    "sample",
    "dataset",
    "species",
    "run_accession",
    "fastq_1",
    "fastq_2",
    "layout",
    "biosample",
    "sra_study",
    "sample_name",
    "biological_replicate",
    "total_usable_reads",
    "reads_passed_to_minimap2",
    "normalization_denominator",
    "normalization_reads",
    "read_count",
    "reads_examined",
    "matching_reads",
    "matching_reads_per_million_examined",
}

MITO_FEATURE_ORDER = {
    name.lower(): i
    for i, name in enumerate(
        [
            "MT-TF",
            "MT-RNR1",
            "MT-TV",
            "MT-RNR2",
            "MT-TL1",
            "MT-ND1",
            "MT-TI",
            "MT-TQ",
            "MT-TM",
            "MT-ND2",
            "MT-TW",
            "MT-TA",
            "MT-TN",
            "MT-TC",
            "MT-TY",
            "MT-CO1",
            "MT-TS1",
            "MT-TD",
            "MT-CO2",
            "MT-TK",
            "MT-ATP8",
            "MT-ATP6",
            "MT-CO3",
            "MT-TG",
            "MT-ND3",
            "MT-TR",
            "MT-ND4L",
            "MT-ND4",
            "MT-TH",
            "MT-TS2",
            "MT-TL2",
            "MT-ND5",
            "MT-ND6",
            "MT-TE",
            "MT-CYB",
            "MT-TT",
            "MT-TP",
        ]
    )
}


def value_columns(matrix: pd.DataFrame, samples: pd.DataFrame | None = None) -> list[str]:
    metadata = set(samples.columns) if samples is not None else set()
    cols = []
    for col in matrix.columns:
        if col in TECHNICAL_COLUMNS or col in metadata:
            continue
        if col.endswith(("_denominator", "_read_count", "_reads_examined")):
            continue
        values = pd.to_numeric(matrix[col], errors="coerce")
        if values.notna().any():
            cols.append(col)
    return cols


def save(fig: plt.Figure, path: str) -> None:
    ensure_parent(path)
    fig.savefig(path, bbox_inches="tight")
    svg = str(Path(path).with_suffix(".svg"))
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def save_multi_page(figures: list[plt.Figure], path: str) -> None:
    ensure_parent(path)
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")


def move_legend_outside(ax: plt.Axes, title: str | None = None) -> None:
    handles, labels = ax.get_legend_handles_labels()
    pairs = [(handle, label) for handle, label in zip(handles, labels) if label and not str(label).startswith("_")]
    if not pairs and ax.legend_:
        handles = getattr(ax.legend_, "legend_handles", getattr(ax.legend_, "legendHandles", []))
        labels = [text.get_text() for text in ax.legend_.get_texts()]
        pairs = [(handle, label) for handle, label in zip(handles, labels) if label and not str(label).startswith("_")]
    if not pairs:
        if ax.legend_:
            ax.legend_.remove()
        return
    ax.legend(
        [handle for handle, _ in pairs],
        [label for _, label in pairs],
        title=title,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        frameon=True,
    )


def empty(path: str, title: str, message: str = "No data available for this plot") -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set_title(title)
    ax.set_axis_off()
    save(fig, path)


def rainfall_support_limits(values: pd.Series | np.ndarray) -> tuple[float, float]:
    support = pd.to_numeric(pd.Series(values), errors="coerce")
    support = support[np.isfinite(support) & (support > 0)]
    if support.empty:
        return 1.0, 1.0
    return float(support.min()), float(support.max())


def rainfall_point_sizes(values: pd.Series | np.ndarray, support_min: float, support_max: float) -> np.ndarray:
    support = pd.to_numeric(pd.Series(values), errors="coerce").fillna(0).to_numpy(dtype=float)
    min_size = 7.0
    max_size = 245.0
    if support_max <= support_min:
        return np.full_like(support, (min_size + max_size) / 2, dtype=float)
    fraction = np.clip((support - support_min) / (support_max - support_min), 0, 1)
    return min_size + (max_size - min_size) * np.power(fraction, 0.92)


def support_legend_values(support_min: float, support_max: float, n: int = 5) -> list[float]:
    if support_max <= 0:
        return []
    if support_min <= 0 or support_min >= support_max:
        return [support_max]
    values = np.geomspace(support_min, support_max, n)
    rounded: list[float] = []
    for value in values:
        if not rounded or not np.isclose(value, rounded[-1], rtol=0.04, atol=0):
            rounded.append(float(value))
    rounded[-1] = support_max
    return rounded


def support_tick_label(value: float) -> str:
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    if value >= 1:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return f"{value:.3g}"


def rainfall_y_axis_min(values: pd.Series | np.ndarray) -> float:
    sizes = pd.to_numeric(pd.Series(values), errors="coerce")
    sizes = sizes[np.isfinite(sizes) & (sizes > 0)]
    if sizes.empty:
        return 0.1
    return max(0.1, float(10 ** np.floor(np.log10(float(sizes.min())))))


def read_tsv_safe(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def normalize_deletion_ids(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    if "exact_deletion_id" not in out.columns and "junction_id" in out.columns:
        out["exact_deletion_id"] = out["junction_id"]
    if "junction_id" not in out.columns and "exact_deletion_id" in out.columns:
        out["junction_id"] = out["exact_deletion_id"]
    return out


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


def compact_label(value: object, max_len: int = 88) -> str:
    text = str(value)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1] + "..."


def compact_feature_label(value: object, max_parts: int = 6, max_len: int = 92) -> str:
    text = str(value)
    parts = [part for part in text.split("+") if part]
    if len(parts) > max_parts:
        text = f"{parts[0]} ... {parts[-1]} ({len(parts)} features)"
    return compact_label(text, max_len=max_len)


def safe_filename(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return text.strip("_") or "group"


def display_label(value: object) -> str:
    text = str(value).replace("_", " ").replace("-", " ").strip()
    if not text:
        return ""
    return " ".join(word.upper() if word.upper() in {"QC", "PCA", "MDS"} else word.capitalize() for word in text.split())


def normalization_mode(table: pd.DataFrame) -> str:
    if "normalization_denominator" not in table.columns:
        return "mt_evidence_reads"
    values = table["normalization_denominator"].dropna().astype(str).unique().tolist()
    return values[0] if values else "mt_evidence_reads"


def per_million_phrase(table: pd.DataFrame) -> str:
    if normalization_mode(table) == "total_usable_reads":
        return "per million usable reads"
    return "per million mitochondrial-evidence reads"


def normalization_denominator_by_sample(table: pd.DataFrame) -> pd.Series | None:
    if "normalization_reads" in table.columns:
        return pd.to_numeric(table.set_index("sample")["normalization_reads"], errors="coerce")
    if "reads_passed_to_minimap2" in table.columns:
        return pd.to_numeric(table.set_index("sample")["reads_passed_to_minimap2"], errors="coerce")
    return None


def mito_order_key(name: object) -> tuple[int, str]:
    text = str(name)
    normalized = text.upper().replace("MT-", "MT-")
    return (MITO_FEATURE_ORDER.get(normalized.lower(), 999), text)


def age_sort_key(value: object) -> tuple[float, str]:
    text = str(value)
    number = pd.to_numeric(pd.Series([text]).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce").iloc[0]
    return (float(number) if pd.notna(number) else 1e9, text.lower())


def treatment_sort_key(value: object) -> tuple[int, str]:
    text = str(value)
    lower = text.lower()
    control_rank = 0 if lower in {"control", "ctrl", "vehicle", "untreated"} or "control" in lower else 1
    return (control_rank, lower)


def ordered_groups(samples: pd.DataFrame, group_col: str) -> list[str]:
    if not group_col or group_col not in samples.columns:
        return []
    work = samples[[group_col] + [col for col in ["age", "treatment"] if col in samples.columns]].copy()
    work[group_col] = work[group_col].fillna("missing").astype(str)
    if {"age", "treatment"}.issubset(work.columns):
        reps = work.drop_duplicates(group_col).copy()
        reps["_age_key"] = reps["age"].map(age_sort_key)
        reps["_treatment_key"] = reps["treatment"].map(treatment_sort_key)
        return reps.sort_values(["_age_key", "_treatment_key", group_col])[group_col].astype(str).tolist()
    return sorted(work[group_col].unique(), key=lambda value: (treatment_sort_key(value), age_sort_key(value), str(value)))


def ordered_ages(samples: pd.DataFrame) -> list[str]:
    return sorted(samples["age"].dropna().astype(str).unique().tolist(), key=age_sort_key) if "age" in samples.columns else []


def ordered_treatments(samples: pd.DataFrame) -> list[str]:
    return sorted(samples["treatment"].dropna().astype(str).unique().tolist(), key=treatment_sort_key) if "treatment" in samples.columns else []


def is_noncontrol(value: object) -> bool:
    return treatment_sort_key(value)[0] > 0


def feature_kind(row: pd.Series) -> str:
    text = " ".join(str(row.get(col, "")) for col in row.index).lower()
    name = str(row.get("feature", row.get("gene_name", row.get("name", "")))).lower()
    if "trna" in text or name.startswith(("mt-t", "trn")):
        return "tRNA"
    if "rrna" in text or name.startswith(("mt-r", "rrn")):
        return "rRNA"
    if "protein_coding" in text or "cds" in text or name.startswith(("mt-co", "mt-cy", "mt-nd", "mt-atp")):
        return "protein-coding"
    return "other"


def feature_name(row: pd.Series) -> str:
    for col in ["gene_name", "feature", "name", "gene_id"]:
        value = str(row.get(col, "")).strip()
        if value and value.lower() != "nan":
            return value
    return ""


def mitochondrial_axis_bounds(features: pd.DataFrame) -> tuple[float, float]:
    """Return a stable full-mtDNA x-axis span for coordinate plots."""
    if features.empty or "end" not in features.columns:
        return (1.0, 1.0)
    ends = pd.to_numeric(features["end"], errors="coerce").dropna()
    if ends.empty:
        return (1.0, 1.0)
    return (1.0, float(max(ends.max(), 1.0)))


def draw_feature_tracks(ax: plt.Axes, features: pd.DataFrame, y_base: float, height: float) -> None:
    if features.empty or "start" not in features.columns or "end" not in features.columns:
        return
    colors = {
        "protein-coding": "#2563eb",
        "rRNA": "#7c3aed",
        "tRNA": "#059669",
        "other": "#64748b",
    }
    lanes = {"protein-coding": 0, "rRNA": 1, "tRNA": 2, "other": 3}
    label_min_width = 360
    for _, row in features.iterrows():
        start = pd.to_numeric(row.get("start"), errors="coerce")
        end = pd.to_numeric(row.get("end"), errors="coerce")
        if pd.isna(start) or pd.isna(end) or end <= start:
            continue
        kind = feature_kind(row)
        lane = lanes.get(kind, 3)
        y = y_base - lane * height * 1.18
        ax.add_patch(
            Rectangle(
                (float(start), y),
                float(end - start),
                height,
                facecolor=colors.get(kind, "#64748b"),
                edgecolor="none",
                alpha=0.22,
                clip_on=True,
            )
        )
        name = feature_name(row)
        if name and end - start >= label_min_width and kind in {"protein-coding", "rRNA"}:
            ax.text(
                (float(start) + float(end)) / 2,
                y + height / 2,
                name,
                ha="center",
                va="center",
                fontsize=6,
                color="#1f2933",
                clip_on=True,
            )
    for kind, lane in lanes.items():
        y = y_base - lane * height * 1.18 + height / 2
        ax.text(
            -0.01,
            y,
            kind,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=6,
            color="#52606d",
            clip_on=False,
        )


def draw_feature_track_axis(ax: plt.Axes, features: pd.DataFrame) -> tuple[float, float]:
    if features.empty or "start" not in features.columns or "end" not in features.columns:
        ax.set_axis_off()
        return (1.0, 1.0)
    x_min, x_max = mitochondrial_axis_bounds(features)
    colors = {
        "protein-coding": "#2563eb",
        "rRNA": "#7c3aed",
        "tRNA": "#059669",
        "other": "#64748b",
    }
    lanes = {"protein-coding": 3, "rRNA": 2, "tRNA": 1, "other": 0}
    label_min_width = 520
    for _, row in features.iterrows():
        start = pd.to_numeric(row.get("start"), errors="coerce")
        end = pd.to_numeric(row.get("end"), errors="coerce")
        if pd.isna(start) or pd.isna(end) or end <= start:
            continue
        kind = feature_kind(row)
        lane = lanes.get(kind, 0)
        ax.add_patch(
            Rectangle(
                (float(start), lane - 0.32),
                float(end - start),
                0.64,
                facecolor=colors.get(kind, "#64748b"),
                edgecolor="none",
                alpha=0.28,
            )
        )
        name = feature_name(row)
        if name and end - start >= label_min_width and kind in {"protein-coding", "rRNA"}:
            ax.text(
                (float(start) + float(end)) / 2,
                lane,
                name,
                ha="center",
                va="center",
                fontsize=6.5,
                color="#1f2933",
                clip_on=True,
            )
    ax.set_yticks(list(lanes.values()))
    ax.set_yticklabels(list(lanes.keys()), fontsize=7)
    ax.set_ylim(-0.7, 3.7)
    ax.set_xlim(x_min, x_max)
    ax.grid(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return (x_min, x_max)


def format_deletion_size_log_axis(ax: plt.Axes) -> None:
    ax.set_yscale("log")
    major_ticks = [0.1, 1, 10, 100, 1000, 10000]
    ax.yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda value, _: f"{int(value):,}" if value >= 1 else f"{value:g}")
    )
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(axis="y", which="major", length=6, width=0.8)
    ax.tick_params(axis="y", which="minor", length=3, width=0.55)
    ax.grid(True, which="major", axis="y", alpha=0.32, linewidth=0.8)
    ax.grid(True, which="minor", axis="y", alpha=0.18, linewidth=0.45)
    ax.grid(True, which="major", axis="x", alpha=0.18, linewidth=0.6)


def palette(samples: pd.DataFrame, group_col: str) -> dict:
    if not group_col or group_col not in samples.columns:
        return {}
    groups = ordered_groups(samples, group_col)
    if {"age", "treatment"}.issubset(samples.columns):
        reps = samples.drop_duplicates(group_col).set_index(group_col)
        control_colors = ["#4c78a8", "#59a14f", "#76b7b2", "#9c9ede"]
        treatment_colors = ["#f28e8b", "#d62728", "#b22222", "#8c1d1d", "#e15759", "#af7aa1"]
        control_i = 0
        treatment_i = 0
        out = {}
        for group in groups:
            treatment = reps.loc[group, "treatment"] if group in reps.index else group
            if is_noncontrol(treatment):
                out[group] = treatment_colors[min(treatment_i, len(treatment_colors) - 1)]
                treatment_i += 1
            else:
                out[group] = control_colors[min(control_i, len(control_colors) - 1)]
                control_i += 1
        return out
    base = ["#4c78a8", "#e15759", "#59a14f", "#af7aa1", "#f28e2b", "#76b7b2", "#edc949", "#9c755f"]
    return {group: base[i % len(base)] for i, group in enumerate(groups)}


def treatment_palette(samples: pd.DataFrame) -> dict:
    treatments = ordered_treatments(samples)
    control_colors = ["#4c78a8", "#59a14f", "#76b7b2"]
    treatment_colors = ["#d62728", "#b22222", "#af7aa1", "#e15759"]
    out = {}
    control_i = 0
    treatment_i = 0
    for treatment in treatments:
        if is_noncontrol(treatment):
            out[treatment] = treatment_colors[min(treatment_i, len(treatment_colors) - 1)]
            treatment_i += 1
        else:
            out[treatment] = control_colors[min(control_i, len(control_colors) - 1)]
            control_i += 1
    return out


def sample_group_order(samples: pd.DataFrame, group_col: str) -> list[str]:
    if group_col and group_col in samples.columns:
        return samples.sort_values([group_col, "sample"])["sample"].tolist()
    return samples.sort_values("sample")["sample"].tolist()


def burden_plot(burden: pd.DataFrame, group_col: str, path: str, y: str, title: str, ylabel: str) -> None:
    if burden.empty or y not in burden.columns:
        empty(path, title, "No sample-level burden values are available")
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    pal = palette(burden, group_col)
    if group_col and group_col in burden.columns:
        order = ordered_groups(burden, group_col)
        sns.stripplot(data=burden, x=group_col, y=y, hue=group_col, order=order, hue_order=order, palette=pal, size=7, jitter=0.18, ax=ax, legend=False)
        group_sizes = burden.groupby(group_col)["sample"].count() if "sample" in burden.columns else burden.groupby(group_col).size()
        show_ci = not group_sizes.empty and int(group_sizes.min()) >= 3
        sns.pointplot(
            data=burden,
            x=group_col,
            y=y,
            order=order,
            color="black",
            errorbar=("ci", 95) if show_ci else None,
            markers="D",
            linestyles="none",
            ax=ax,
        )
        ax.scatter([], [], color="#4c78a8", s=45, label="sample")
        ax.scatter([], [], color="black", marker="D", s=55, label="group mean" + (", 95% CI" if show_ci else ""))
        move_legend_outside(ax)
        ax.set_xlabel(display_label(group_col))
    else:
        order = sample_group_order(burden, group_col)
        sns.barplot(data=burden, x="sample", y=y, order=order, color="#4c78a8", ax=ax)
        ax.tick_params(axis="x", rotation=60)
        ax.set_xlabel("sample")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    save(fig, path)


def factorial_interaction_plot(burden: pd.DataFrame, path: str, y: str, title: str, ylabel: str) -> None:
    if burden.empty or y not in burden.columns or not {"age", "treatment"}.issubset(burden.columns):
        empty(path, title, "This plot requires age and treatment metadata")
        return
    ages = ordered_ages(burden)
    treatments = ordered_treatments(burden)
    if len(ages) < 2 or len(treatments) < 2:
        empty(path, title, "This plot requires at least two ages and two treatment groups")
        return
    fig, ax = plt.subplots(figsize=(8.5, 5))
    pal = treatment_palette(burden)
    sns.stripplot(
        data=burden,
        x="age",
        y=y,
        hue="treatment",
        order=ages,
        hue_order=treatments,
        dodge=True,
        jitter=0.12,
        size=7,
        palette=pal,
        ax=ax,
    )
    summary = burden.groupby(["age", "treatment"], as_index=False)[y].mean()
    sns.pointplot(
        data=summary,
        x="age",
        y=y,
        hue="treatment",
        order=ages,
        hue_order=treatments,
        dodge=0.42,
        markers="D",
        linestyles="-",
        errorbar=None,
        palette=pal,
        ax=ax,
        legend=False,
    )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[: len(treatments)], labels[: len(treatments)], title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, frameon=True)
    ax.set_xlabel("Age")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    save(fig, path)


def size_distribution(reads: pd.DataFrame, samples: pd.DataFrame, group_col: str, path: str, title: str, weighted: bool, log_y: bool = False, size_min: int | None = None, size_max: int | None = None) -> None:
    if reads.empty or "deleted_size" not in reads.columns:
        empty(path, title, "No deletion-supporting reads are available for this size distribution")
        return
    df = reads.copy()
    df["deleted_size"] = pd.to_numeric(df["deleted_size"], errors="coerce")
    df = df.dropna(subset=["deleted_size"])
    if size_min is not None:
        df = df[df["deleted_size"] >= size_min]
    if size_max is not None:
        df = df[df["deleted_size"] <= size_max]
    if df.empty:
        empty(path, title, "No deletion-supporting reads fall in this size range")
        return
    df = df.merge(samples[["sample", group_col]] if group_col in samples.columns else samples[["sample"]], on="sample", how="left")
    fig, ax = plt.subplots(figsize=(10, 5))
    stat = "count"
    weights = None
    ylabel = "Deletion-supporting read count"
    if weighted:
        denom = normalization_denominator_by_sample(samples)
        if denom is not None:
            df["weight"] = df["sample"].map(lambda s: 1_000_000 / denom.get(s, np.nan) if denom.get(s, 0) > 0 else 0)
            weights = df["weight"]
            ylabel = f"Support {per_million_phrase(samples)}"
    hue = group_col if group_col in df.columns else None
    sns.histplot(
        data=df,
        x="deleted_size",
        hue=hue,
        hue_order=ordered_groups(samples, group_col) if hue else None,
        palette=palette(samples, group_col) if hue else None,
        weights=weights,
        bins=60,
        element="step",
        stat=stat,
        common_norm=False,
        ax=ax,
    )
    if log_y:
        ax.set_yscale("log")
        ax.set_ylim(bottom=max(0.1, ax.get_ylim()[0]))
    if hue:
        move_legend_outside(ax, title=display_label(group_col))
    ax.set_xlabel("Deleted size (bp)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    save(fig, path)


def rainfall(reads: pd.DataFrame, samples: pd.DataFrame, features: pd.DataFrame, group_col: str, path: str) -> None:
    if reads.empty:
        empty(path, "Deletion Position/Size Abundance", "No deletion-supporting reads are available for position/size plotting")
        return
    df = reads.copy()
    df["left_breakpoint"] = pd.to_numeric(df["left_breakpoint"], errors="coerce")
    df["right_breakpoint"] = pd.to_numeric(df["right_breakpoint"], errors="coerce")
    df["deleted_size"] = pd.to_numeric(df["deleted_size"], errors="coerce")
    df = df.dropna(subset=["left_breakpoint", "right_breakpoint", "deleted_size"])
    df["midpoint"] = (df["left_breakpoint"] + df["right_breakpoint"]) / 2
    sample_cols = ["sample", group_col] if group_col in samples.columns else ["sample"]
    for col in ["normalization_denominator", "normalization_reads", "reads_passed_to_minimap2"]:
        if col in samples.columns:
            sample_cols.append(col)
    df = df.merge(samples[sample_cols], on="sample", how="left")
    if group_col in df.columns:
        df["_plot_group"] = df[group_col].fillna("missing").astype(str)
    else:
        df["_plot_group"] = "all"
    if "normalization_reads" in df.columns or "reads_passed_to_minimap2" in df.columns:
        denom_col = "normalization_reads" if "normalization_reads" in df.columns else "reads_passed_to_minimap2"
        denom = pd.to_numeric(df[denom_col], errors="coerce")
        df["_support_weight"] = np.where(denom > 0, 1_000_000 / denom, 0)
        support_col = "support_per_million_mt_reads"
        support_label = f"support {per_million_phrase(samples)}"
    else:
        df["_support_weight"] = 1.0
        support_col = "supporting_reads"
        support_label = "supporting reads"
    grouped = (
        df.groupby(["_plot_group", "left_breakpoint", "right_breakpoint", "deleted_size"], as_index=False)
        .agg(supporting_reads=("sample", "size"), support_per_million_mt_reads=("_support_weight", "sum"))
    )
    grouped["midpoint"] = (grouped["left_breakpoint"] + grouped["right_breakpoint"]) / 2
    grouped["_plot_support"] = pd.to_numeric(grouped[support_col], errors="coerce").fillna(0)
    groups = [group for group in ordered_groups(samples, group_col) if group in set(grouped["_plot_group"])] or sorted(grouped["_plot_group"].unique())
    for old in Path(path).parent.glob(f"{Path(path).stem}__*.pdf"):
        old.unlink()
    for old in Path(path).parent.glob(f"{Path(path).stem}__*.svg"):
        old.unlink()
    support_min, support_max = rainfall_support_limits(grouped["_plot_support"])
    if support_min < support_max:
        support_norm: colors.Normalize = colors.LogNorm(vmin=support_min, vmax=support_max)
    else:
        support_norm = colors.Normalize(vmin=0, vmax=max(support_max, 1.0))
    legend_values = support_legend_values(support_min, support_max)
    y_axis_min = rainfall_y_axis_min(grouped["deleted_size"])
    figures: list[plt.Figure] = []
    sidecars = []
    for group_index, group in enumerate(groups):
        sub = grouped[grouped["_plot_group"] == group].sort_values("_plot_support", ascending=False).head(300)
        fig = plt.figure(figsize=(14.4, 6.9), constrained_layout=True)
        grid = fig.add_gridspec(
            2,
            2,
            width_ratios=[1.0, 0.38],
            height_ratios=[5.2, 1.0],
            hspace=0.05,
            wspace=0.14,
        )
        ax = fig.add_subplot(grid[0, 0])
        feature_ax = fig.add_subplot(grid[1, 0], sharex=ax)
        legend_grid = grid[:, 1].subgridspec(
            6,
            1,
            height_ratios=[0.34, 0.62, 0.34, 0.18, 1.18, 2.6],
            hspace=0.18,
        )
        scatter = None
        if sub.empty:
            ax.text(0.5, 0.5, "No deletion-supporting reads in this group", ha="center", va="center", transform=ax.transAxes)
        else:
            sizes = rainfall_point_sizes(sub["_plot_support"], support_min, support_max)
            scatter = ax.scatter(
                sub["midpoint"],
                sub["deleted_size"],
                c=sub["_plot_support"],
                s=sizes,
                cmap="magma",
                norm=support_norm,
                alpha=0.78,
                edgecolors="#1f2933",
                linewidths=0.35,
            )
        ax.set_title(f"Deletion Position/Size Abundance: {group}")
        ax.set_ylabel("Deleted size (bp)")
        y_min = y_axis_min
        y_max = max(10_000, float(sub["deleted_size"].max()) * 1.18) if not sub.empty else 10_000
        ax.set_ylim(y_min, y_max)
        format_deletion_size_log_axis(ax)
        x_min, x_max = draw_feature_track_axis(feature_ax, features)
        ax.set_xlim(x_min, x_max)
        feature_ax.set_xlim(x_min, x_max)
        feature_ax.set_xlabel("Deletion midpoint on mitochondrial genome (bp)")
        if scatter is not None:
            color_title_ax = fig.add_subplot(legend_grid[0, 0])
            color_title_ax.set_axis_off()
            color_title_ax.text(0.5, 0.5, "Point color", ha="center", va="center", fontsize=9, fontweight="bold")

            color_label_ax = fig.add_subplot(legend_grid[1, 0])
            color_label_ax.set_axis_off()
            color_label = textwrap.fill(support_label, width=42)
            color_label_ax.text(0.5, 0.5, color_label, ha="center", va="center", fontsize=8)

            cbar_ax = fig.add_subplot(legend_grid[2, 0])
            cbar = fig.colorbar(scatter, cax=cbar_ax, orientation="horizontal")
            if legend_values:
                cbar.set_ticks(legend_values)
                cbar.set_ticklabels([support_tick_label(value) for value in legend_values])
            cbar.ax.tick_params(labelsize=8, length=3)
            cbar.outline.set_linewidth(0.5)
            handles = []
            for value in legend_values:
                if value <= 0:
                    continue
                legend_color = scatter.cmap(scatter.norm(value))
                handles.append(
                    plt.Line2D(
                        [],
                        [],
                        linestyle="",
                        marker="o",
                        markersize=np.sqrt(rainfall_point_sizes([value], support_min, support_max)[0]) / 1.45,
                        markerfacecolor=legend_color,
                        markeredgecolor="#1f2933",
                        alpha=0.7,
                        label=support_tick_label(value),
                    )
                )
            if handles:
                size_ax = fig.add_subplot(legend_grid[4, 0])
                size_ax.set_axis_off()
                size_ax.legend(
                    handles=handles,
                    title=f"Point size:\n{support_label}",
                    loc="center",
                    borderaxespad=0,
                    fontsize=8,
                    title_fontsize=8,
                    frameon=True,
                )
        else:
            blank_legend_ax = fig.add_subplot(grid[:, 1])
            blank_legend_ax.set_axis_off()
        sidecar_pdf = Path(path).with_name(f"{Path(path).stem}__{safe_filename(group)}.pdf")
        sidecar_svg = sidecar_pdf.with_suffix(".svg")
        ensure_parent(sidecar_pdf)
        fig.savefig(sidecar_pdf, bbox_inches="tight")
        fig.savefig(sidecar_svg, bbox_inches="tight")
        if group_index == 0:
            fig.savefig(path, bbox_inches="tight")
            fig.savefig(Path(path).with_suffix(".svg"), bbox_inches="tight")
        figures.append(fig)
        sidecars.append(sidecar_pdf)
    save_multi_page(figures, path)
    for fig in figures:
        plt.close(fig)


def category_bar(matrix: pd.DataFrame, samples: pd.DataFrame, group_col: str, path: str, title: str, ylabel: str, top_n: int = 14, proportional: bool = False) -> None:
    feature_cols = value_columns(matrix, samples)
    if matrix.empty or not feature_cols:
        empty(path, title, "No affected-feature categories are available for this plot")
        return
    df = matrix.merge(samples[["sample", group_col]] if group_col in samples.columns else samples[["sample"]], on="sample", how="left", suffixes=("", "_sample"))
    group_key = group_col if group_col in df.columns else None
    long = df.melt(id_vars=[group_key] if group_key else [], value_vars=feature_cols, var_name="category", value_name="value")
    long["value"] = pd.to_numeric(long["value"], errors="coerce").fillna(0)
    if group_key:
        grouped = long.groupby([group_key, "category"], as_index=False)["value"].sum()
        if proportional:
            grouped["value"] = grouped["value"] / grouped.groupby(group_key)["value"].transform("sum").replace(0, np.nan) * 100
            ylabel = "Percent of deletion support"
    else:
        grouped = long.groupby(["category"], as_index=False)["value"].sum()
    top = grouped.groupby("category")["value"].sum().sort_values(ascending=False).head(top_n).index
    if len(top) < grouped["category"].nunique():
        other = grouped[~grouped["category"].isin(top)].copy()
        grouped = grouped[grouped["category"].isin(top)].copy()
        if not other.empty:
            if group_key:
                other = other.groupby(group_key, as_index=False)["value"].sum()
                other["category"] = "Other categories"
            else:
                other = pd.DataFrame([{"category": "Other categories", "value": other["value"].sum()}])
            grouped = pd.concat([grouped, other], ignore_index=True)
    else:
        grouped = grouped[grouped["category"].isin(top)].copy()
    grouped["display_category"] = grouped["category"].map(compact_feature_label)
    fig, ax = plt.subplots(figsize=(12, max(5, 0.42 * len(top))))
    sns.barplot(
        data=grouped,
        y="display_category",
        x="value",
        hue=group_key,
        hue_order=ordered_groups(samples, group_col) if group_key else None,
        palette=palette(samples, group_col),
        ax=ax,
    )
    ax.set_xlabel(ylabel)
    ax.set_ylabel("")
    ax.set_title(title)
    if ax.legend_:
        move_legend_outside(ax, title=display_label(group_key))
    save(fig, path)


def per_gene_plot(per_gene: pd.DataFrame, features: pd.DataFrame, group_col: str, path: str, ylabel: str) -> None:
    if per_gene.empty:
        empty(path, "Per-Gene Affected Burden", "No per-feature deletion burden values are available")
        return
    group_key = group_col if group_col in per_gene.columns else None
    grouped = per_gene.groupby(([group_key] if group_key else []) + ["feature"], as_index=False)["support_per_million_mt_reads"].sum()
    top = grouped.groupby("feature")["support_per_million_mt_reads"].sum().sort_values(ascending=False).head(25).index
    grouped = grouped[grouped["feature"].isin(top)]
    order = None
    if not features.empty and {"gene_name", "start"}.issubset(features.columns):
        f = features.copy()
        f["start"] = pd.to_numeric(f["start"], errors="coerce")
        f = f.dropna(subset=["start"])
        order = [name for name in f.sort_values("start")["gene_name"].astype(str).drop_duplicates().tolist() if name in set(grouped["feature"])]
    if not order:
        order = sorted(grouped["feature"].drop_duplicates().tolist(), key=mito_order_key)
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.3 * len(top))))
    sns.barplot(
        data=grouped,
        y="feature",
        x="support_per_million_mt_reads",
        hue=group_key,
        hue_order=ordered_groups(per_gene, group_col) if group_key else None,
        order=order,
        palette=palette(per_gene, group_col) if group_key else None,
        ax=ax,
    )
    ax.set_xlabel(ylabel)
    ax.set_ylabel("")
    ax.set_title("Per-Gene Affected Burden")
    if ax.legend_:
        move_legend_outside(ax, title=display_label(group_key))
    save(fig, path)


def exact_recurrence(clusters: pd.DataFrame, exact_mtpm: pd.DataFrame, samples: pd.DataFrame, group_col: str, path: str, ylabel: str) -> None:
    if clusters.empty:
        empty(path, "Exact Deletion Recurrence", "No exact deletion calls are available")
        return
    df = clusters.copy()
    df["left_breakpoint"] = pd.to_numeric(df["left_breakpoint"], errors="coerce").astype("Int64")
    df["right_breakpoint"] = pd.to_numeric(df["right_breakpoint"], errors="coerce").astype("Int64")
    df["deleted_size"] = pd.to_numeric(df["deleted_size"], errors="coerce").astype("Int64")
    feature_label = df.get("affected_feature_label", pd.Series([""] * len(df), index=df.index)).map(compact_feature_label)
    df["label"] = (
        df["left_breakpoint"].astype(str)
        + "->"
        + df["right_breakpoint"].astype(str)
        + " ("
        + df["deleted_size"].astype(str)
        + " bp) | "
        + feature_label.astype(str)
    ).map(lambda value: compact_label(value, max_len=110))
    df = df.sort_values("total_supporting_reads", ascending=False).head(22)
    if not exact_mtpm.empty and group_col in samples.columns:
        value_cols = [col for col in df["exact_deletion_id"].astype(str).tolist() if col in exact_mtpm.columns]
        if value_cols:
            merged = exact_mtpm[["sample", *value_cols]].merge(samples[["sample", group_col]], on="sample", how="left")
            long = merged.melt(id_vars=["sample", group_col], value_vars=value_cols, var_name="exact_deletion_id", value_name="support_per_million_mt_reads")
            grouped = long.groupby([group_col, "exact_deletion_id"], as_index=False)["support_per_million_mt_reads"].sum()
            label_map = df.set_index("exact_deletion_id")["label"].to_dict()
            grouped["label"] = grouped["exact_deletion_id"].map(label_map)
            piv = grouped.pivot_table(index="label", columns=group_col, values="support_per_million_mt_reads", aggfunc="sum", fill_value=0)
            order_labels = df["label"].tolist()
            piv = piv.reindex(order_labels)
            group_order = [group for group in ordered_groups(samples, group_col) if group in piv.columns]
            fig, ax = plt.subplots(figsize=(12, max(5, 0.38 * len(df))))
            pal = palette(samples, group_col)
            y = np.arange(len(piv))
            bar_height = min(0.82 / max(1, len(group_order)), 0.22)
            offsets = (np.arange(len(group_order)) - (len(group_order) - 1) / 2) * bar_height
            for offset, group in zip(offsets, group_order):
                values = piv[group].to_numpy(dtype=float)
                ax.barh(y + offset, values, height=bar_height * 0.9, color=pal.get(group), label=group)
            ax.set_yticks(y)
            ax.set_yticklabels(piv.index)
            ax.invert_yaxis()
            ax.set_xlabel(ylabel)
            ax.legend(title=display_label(group_col), loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, frameon=True)
            ax.set_ylabel("")
            ax.set_title("Exact Deletion Recurrence")
            save(fig, path)
            return
    fig, ax = plt.subplots(figsize=(12, max(5, 0.38 * len(df))))
    sns.barplot(data=df, y="label", x="total_supporting_reads", color="#4c78a8", ax=ax)
    ax.set_xlabel("Total supporting reads")
    ax.set_ylabel("")
    ax.set_title("Exact Deletion Recurrence")
    save(fig, path)


def ordination(matrix: pd.DataFrame, samples: pd.DataFrame, group_col: str, path: str, title: str, method: str) -> None:
    feature_cols = value_columns(matrix, samples)
    if len(matrix) < 2:
        empty(path, title, f"{title} requires at least two samples; this dataset has {len(matrix)}")
        return
    if len(feature_cols) < 2:
        empty(path, title, f"{title} requires at least two deletion features with numeric support")
        return
    x = matrix[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    if np.all(x == 0):
        empty(path, title, "All deletion feature values are zero, so ordination is not informative")
        return
    if method == "pca":
        coords = PCA(n_components=2, random_state=1).fit_transform(x)
        xlabel, ylabel = "PC1", "PC2"
    else:
        dist = squareform(pdist(x, metric="braycurtis"))
        dist = np.nan_to_num(dist, nan=0.0, posinf=0.0, neginf=0.0)
        coords = MDS(
            n_components=2,
            metric="precomputed",
            random_state=1,
            n_init=4,
            init="random",
            normalized_stress="auto",
        ).fit_transform(dist)
        xlabel, ylabel = "MDS1", "MDS2"
    df = samples.copy()
    df[xlabel] = coords[:, 0]
    df[ylabel] = coords[:, 1]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sns.scatterplot(
        data=df,
        x=xlabel,
        y=ylabel,
        hue=group_col if group_col in df.columns else None,
        hue_order=ordered_groups(samples, group_col) if group_col in df.columns else None,
        palette=palette(samples, group_col),
        s=85,
        ax=ax,
    )
    ax.set_title(title)
    if ax.legend_:
        move_legend_outside(ax, title=display_label(group_col))
    save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--all-reads", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--burden", required=True)
    parser.add_argument("--exact-mtpm", required=True)
    parser.add_argument("--affected-raw", required=True)
    parser.add_argument("--affected-mtpm", required=True)
    parser.add_argument("--impact-class-mtpm", required=True)
    parser.add_argument("--per-gene-burden", required=True)
    parser.add_argument("--exact-comparison", required=True)
    parser.add_argument("--group-column", default="")
    parser.add_argument("--out-burden", required=True)
    parser.add_argument("--out-unique-count", required=True)
    parser.add_argument("--out-burden-factorial", required=True)
    parser.add_argument("--out-unique-factorial", required=True)
    parser.add_argument("--out-size-unweighted", required=True)
    parser.add_argument("--out-size-weighted", required=True)
    parser.add_argument("--out-size-weighted-log", required=True)
    parser.add_argument("--out-size-small", required=True)
    parser.add_argument("--out-size-medium", required=True)
    parser.add_argument("--out-size-large", required=True)
    parser.add_argument("--out-rainfall", required=True)
    parser.add_argument("--out-affected-support", required=True)
    parser.add_argument("--out-affected-counts", required=True)
    parser.add_argument("--out-affected-proportions", required=True)
    parser.add_argument("--out-impact-class", required=True)
    parser.add_argument("--out-per-gene", required=True)
    parser.add_argument("--out-exact-recurrence", required=True)
    parser.add_argument("--out-exact-pca", required=True)
    parser.add_argument("--out-exact-mds", required=True)
    parser.add_argument("--out-affected-pca", required=True)
    parser.add_argument("--out-affected-mds", required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t")
    features = read_tsv_safe(args.features)
    reads = read_tsv_safe(args.all_reads)
    reads = deduplicate_evidence_reads(normalize_deletion_ids(reads))
    clusters = normalize_deletion_ids(read_tsv_safe(args.clusters))
    burden = read_tsv_safe(args.burden)
    if burden.empty:
        burden = samples.copy()
    exact_mtpm = read_tsv_safe(args.exact_mtpm)
    affected_raw = read_tsv_safe(args.affected_raw)
    affected_mtpm = read_tsv_safe(args.affected_mtpm)
    impact_mtpm = read_tsv_safe(args.impact_class_mtpm)
    per_gene = read_tsv_safe(args.per_gene_burden)
    comparison = normalize_deletion_ids(read_tsv_safe(args.exact_comparison))
    support_label = f"Support {per_million_phrase(burden)}"
    burden_label = f"Deletion-supporting reads {per_million_phrase(burden)}"

    burden_plot(burden, args.group_column, args.out_burden, "deletion_support_per_million_mt_reads", "Total Deletion Burden", burden_label)
    burden_plot(burden, args.group_column, args.out_unique_count, "unique_exact_deletions", "Distinct Exact Deletions", "Distinct exact deletion calls")
    factorial_interaction_plot(burden, args.out_burden_factorial, "deletion_support_per_million_mt_reads", "Deletion Burden: Age By Treatment", burden_label)
    factorial_interaction_plot(burden, args.out_unique_factorial, "unique_exact_deletions", "Distinct Exact Deletions: Age By Treatment", "Distinct exact deletion calls")
    size_distribution(reads, burden, args.group_column, args.out_size_unweighted, "Deletion Size Distribution, Unweighted", weighted=False)
    size_distribution(reads, burden, args.group_column, args.out_size_weighted, "Deletion Size Distribution, Support-Weighted", weighted=True)
    size_distribution(reads, burden, args.group_column, args.out_size_weighted_log, "Deletion Size Distribution, Support-Weighted Log Scale", weighted=True, log_y=True)
    size_distribution(reads, burden, args.group_column, args.out_size_small, "Small Deletion Size Distribution (<1 kb)", weighted=True, size_max=999)
    size_distribution(reads, burden, args.group_column, args.out_size_medium, "Medium Deletion Size Distribution (1-5 kb)", weighted=True, size_min=1000, size_max=4999)
    size_distribution(reads, burden, args.group_column, args.out_size_large, "Large Deletion Size Distribution (>=5 kb)", weighted=True, size_min=5000)
    rainfall(reads, burden, features, args.group_column, args.out_rainfall)
    category_bar(affected_mtpm, samples, args.group_column, args.out_affected_support, "Affected Features: Normalized Abundance", support_label)
    category_bar(affected_raw, samples, args.group_column, args.out_affected_counts, "Affected Features: Raw Supporting Reads", "Supporting reads")
    category_bar(affected_mtpm, samples, args.group_column, args.out_affected_proportions, "Affected Features: Within-Group Percent", "Percent", proportional=True)
    category_bar(impact_mtpm, samples, args.group_column, args.out_impact_class, "Collapsed Feature-Impact Classes", support_label)
    per_gene_plot(per_gene, features, args.group_column, args.out_per_gene, support_label)
    exact_recurrence(clusters, exact_mtpm, samples, args.group_column, args.out_exact_recurrence, support_label)
    ordination(exact_mtpm, samples, args.group_column, args.out_exact_pca, "Exact Deletion PCA", "pca")
    ordination(exact_mtpm, samples, args.group_column, args.out_exact_mds, "Exact Deletion Bray-Curtis MDS", "mds")
    ordination(affected_mtpm, samples, args.group_column, args.out_affected_pca, "Affected-Feature PCA", "pca")
    ordination(affected_mtpm, samples, args.group_column, args.out_affected_mds, "Affected-Feature Bray-Curtis MDS", "mds")


if __name__ == "__main__":
    main()
