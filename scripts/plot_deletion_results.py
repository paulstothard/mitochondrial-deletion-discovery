#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib import colors, ticker
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.textpath import TextPath
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from annotate_junctions import apply_feature_aliases
from circular_deletions import circular_distance
from common import ensure_parent


sns.set_theme(style="whitegrid", context="notebook")

# Shared mitochondrial annotation palette used in every coordinate-based feature depiction.
# The class colors match the report's circular annotation ring.
MITOCHONDRIAL_FEATURE_COLORS = {
    "protein_coding": "#7CAE00",
    "rRNA": "#00BFC4",
    "tRNA": "#C77CFF",
    "region": "#F8766D",
    "other": "#9AA3AF",
}


def mitochondrial_feature_color(feature_class: str) -> str:
    canonical = "protein_coding" if feature_class == "protein-coding" else feature_class
    return MITOCHONDRIAL_FEATURE_COLORS.get(canonical, MITOCHONDRIAL_FEATURE_COLORS["other"])

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
        if col.endswith(("_denominator", "_read_count", "_reads_examined", "_reads", "_read_total")):
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


def save_sample_points_interactive(
    fig: plt.Figure,
    ax: plt.Axes,
    path: str,
    points: list[dict[str, object]],
    plot_type: str = "sample-points",
) -> None:
    """Save a sample-point plot with transparent SVG hit targets and metadata."""
    ensure_parent(path)
    fig.canvas.draw()
    fig.savefig(path, bbox_inches="tight")
    svg_path = Path(path).with_suffix(".svg")
    temporary = svg_path.with_suffix(".tmp.svg")
    fig.savefig(temporary, format="svg")
    svg = temporary.read_text(encoding="utf-8")
    temporary.unlink(missing_ok=True)
    viewbox_match = re.search(
        r'<svg\b[^>]*\bviewBox="([0-9.eE+-]+) ([0-9.eE+-]+) '
        r'([0-9.eE+-]+) ([0-9.eE+-]+)"',
        svg,
    )
    root_match = re.search(r"<svg\b[^>]*>", svg)
    if not viewbox_match or not root_match:
        svg_path.write_text(svg, encoding="utf-8")
        plt.close(fig)
        return
    viewbox_width = float(viewbox_match.group(3))
    viewbox_height = float(viewbox_match.group(4))
    canvas_width, canvas_height = fig.canvas.get_width_height()
    scale_x = viewbox_width / float(canvas_width)
    scale_y = viewbox_height / float(canvas_height)
    circles = []
    for point in points:
        try:
            x_value = float(point.get("x", point.get("x_value")))
            y_value = float(point.get("y", point.get("y_value")))
        except (TypeError, ValueError, KeyError):
            continue
        display_x, display_y = ax.transData.transform((x_value, y_value))
        attrs = {
            "class": "sample-point",
            "cx": f"{display_x * scale_x:.3f}",
            "cy": f"{viewbox_height - display_y * scale_y:.3f}",
            "r": f"{max(5.5, np.sqrt(95.0 / np.pi)) * scale_x:.3f}",
            "fill": "#285f8f",
            "fill-opacity": "0",
            "stroke": "transparent",
            "data-sample": point.get("sample", ""),
            "data-group": point.get("group", ""),
            "data-x-label": point.get("x_label", "X"),
            "data-x-value": point.get("x_value", x_value),
            "data-y-label": point.get("y_label", "Y"),
            "data-y-value": point.get("y_value", y_value),
            "data-biological-replicate": point.get("biological_replicate", ""),
            "data-layout": point.get("layout", ""),
            "data-tissue": point.get("tissue", ""),
            "data-age": point.get("age", ""),
            "data-treatment": point.get("treatment", ""),
        }
        rendered = " ".join(f'{key}="{svg_attribute_value(value)}"' for key, value in attrs.items())
        circles.append(f"<circle {rendered}/>")
    root_end = root_match.end() - 1
    metadata = ' data-plot-type="{}" data-point-count="{}"'.format(
        svg_attribute_value(plot_type), len(circles)
    )
    svg = svg[:root_end] + metadata + svg[root_end:]
    style = '<style>.sample-point{cursor:help;}.sample-point:hover{stroke:#172b4d;stroke-width:2;fill-opacity:.16;}</style>'
    svg = svg[: root_end + len(metadata) + 1] + style + svg[root_end + len(metadata) + 1 :]
    svg = svg.replace("</svg>", '<g id="sample-interactive-points">' + "".join(circles) + "</g></svg>")
    svg_path.write_text(svg, encoding="utf-8")
    plt.close(fig)


def save_ordination_interactive(
    fig: plt.Figure,
    ax: plt.Axes,
    path: str,
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_col: str,
) -> None:
    """Save an ordination with transparent sample hit targets in the SVG report view."""
    ensure_parent(path)
    fig.canvas.draw()
    fig.savefig(path, bbox_inches="tight")
    svg_path = Path(path).with_suffix(".svg")
    temporary = svg_path.with_suffix(".tmp.svg")
    # Keep the SVG in the figure's native coordinate system. The report embeds
    # this sidecar separately from the tight-cropped PDF, so this keeps display
    # coordinates aligned even when the legend is outside the axes.
    fig.savefig(temporary, format="svg")
    svg = temporary.read_text(encoding="utf-8")
    temporary.unlink(missing_ok=True)
    viewbox_match = re.search(
        r'<svg\b[^>]*\bviewBox="([0-9.eE+-]+) ([0-9.eE+-]+) '
        r'([0-9.eE+-]+) ([0-9.eE+-]+)"',
        svg,
    )
    root_match = re.search(r"<svg\b[^>]*>", svg)
    if not viewbox_match or not root_match:
        svg_path.write_text(svg, encoding="utf-8")
        plt.close(fig)
        return
    viewbox_width = float(viewbox_match.group(3))
    viewbox_height = float(viewbox_match.group(4))
    canvas_width, canvas_height = fig.canvas.get_width_height()
    scale_x = viewbox_width / float(canvas_width)
    scale_y = viewbox_height / float(canvas_height)
    circles = []
    for _, row in data.iterrows():
        x_value = float(row[x_col])
        y_value = float(row[y_col])
        display_x, display_y = ax.transData.transform((x_value, y_value))
        group = row.get(group_col, "") if group_col else ""
        group = "" if pd.isna(group) else str(group)
        attrs = {
            "class": "ordination-point",
            "cx": f"{display_x * scale_x:.3f}",
            "cy": f"{viewbox_height - display_y * scale_y:.3f}",
            "r": f"{max(4.5, np.sqrt(85.0 / np.pi)) * scale_x:.3f}",
            "fill": "#285f8f",
            "fill-opacity": "0",
            "stroke": "transparent",
            "data-sample": row.get("sample", ""),
            "data-group": group,
            "data-x-label": x_col,
            "data-y-label": y_col,
            "data-x-value": x_value,
            "data-y-value": y_value,
            "data-biological-replicate": row.get("biological_replicate", ""),
            "data-layout": row.get("layout", ""),
            "data-tissue": row.get("tissue", ""),
        }
        rendered = " ".join(f'{key}="{svg_attribute_value(value)}"' for key, value in attrs.items())
        circles.append(f"<circle {rendered}/>")
    metadata = ' data-plot-type="ordination" data-point-count="{}"'.format(len(circles))
    root_end = root_match.end() - 1
    svg = svg[:root_end] + metadata + svg[root_end:]
    style = '<style>.ordination-point{cursor:help;}.ordination-point:hover{stroke:#172b4d;stroke-width:2;fill-opacity:.16;}</style>'
    svg = svg[: root_end + len(metadata) + 1] + style + svg[root_end + len(metadata) + 1 :]
    svg = svg.replace("</svg>", '<g id="ordination-interactive-points">' + "".join(circles) + "</g></svg>")
    svg_path.write_text(svg, encoding="utf-8")
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


def gene_pair_pca_enabled(config: dict) -> bool:
    """Return whether the configured evidence model supports STAR gene-pair PCA."""
    return bool(
        (config.get("quality", {}) or {})
        .get("short_read_rna_dual_caller", {})
        .get("enabled", False)
    )


def rainfall_support_limits(values: pd.Series | np.ndarray) -> tuple[float, float]:
    support = pd.to_numeric(pd.Series(values), errors="coerce")
    support = support[np.isfinite(support) & (support > 0)]
    if support.empty:
        return 1.0, 1.0
    return float(support.min()), float(support.max())


def rainfall_point_sizes(values: pd.Series | np.ndarray, support_min: float, support_max: float) -> np.ndarray:
    support = pd.to_numeric(pd.Series(values), errors="coerce").fillna(0).to_numpy(dtype=float)
    min_size = 1.2
    max_size = 600.0
    if support_max <= support_min:
        return np.full_like(support, (min_size + max_size) / 2, dtype=float)
    fraction = np.clip((support - support_min) / (support_max - support_min), 0, 1)
    return min_size + (max_size - min_size) * fraction


def nice_support_value(value: float) -> float:
    if value <= 0 or not np.isfinite(value):
        return 0.0
    exponent = np.floor(np.log10(value))
    scaled = value / (10**exponent)
    if scaled <= 1.5:
        mantissa = 1.0
    elif scaled <= 3.5:
        mantissa = 2.0
    elif scaled <= 7.5:
        mantissa = 5.0
    else:
        mantissa = 10.0
    return float(mantissa * (10**exponent))


def nice_support_floor(value: float) -> float:
    if value <= 0 or not np.isfinite(value):
        return 0.0
    exponent = np.floor(np.log10(value))
    scaled = value / (10**exponent)
    if scaled < 2:
        mantissa = 1.0
    elif scaled < 5:
        mantissa = 2.0
    else:
        mantissa = 5.0
    return float(mantissa * (10**exponent))


def nice_support_ceiling(value: float) -> float:
    if value <= 0 or not np.isfinite(value):
        return 0.0
    exponent = np.floor(np.log10(value))
    scaled = value / (10**exponent)
    if scaled <= 1:
        mantissa = 1.0
    elif scaled <= 2:
        mantissa = 2.0
    elif scaled <= 5:
        mantissa = 5.0
    else:
        mantissa = 10.0
    return float(mantissa * (10**exponent))


def support_scale_limits(support_min: float, support_max: float) -> tuple[float, float]:
    if support_max <= 0:
        return 1.0, 1.0
    if support_min <= 0 or support_min >= support_max:
        value = nice_support_value(support_max)
        return value, value
    return nice_support_floor(support_min), nice_support_ceiling(support_max)


def nice_support_candidates(support_min: float, support_max: float) -> list[float]:
    if support_max <= 0 or support_min <= 0:
        return []
    start_exp = int(np.floor(np.log10(support_min))) - 1
    end_exp = int(np.ceil(np.log10(support_max))) + 1
    candidates = []
    for exponent in range(start_exp, end_exp + 1):
        for mantissa in (1, 2, 5):
            value = float(mantissa * (10**exponent))
            if support_min <= value <= support_max:
                candidates.append(value)
    return sorted(set(candidates))


def support_legend_values(support_min: float, support_max: float) -> list[float]:
    if support_max <= 0:
        return []
    if support_min <= 0 or support_min >= support_max:
        return [nice_support_value(support_max)]
    candidates = nice_support_candidates(support_min, support_max)
    if not candidates:
        return sorted(set([nice_support_value(support_min), nice_support_value(support_max)]))
    order_span = np.log10(support_max) - np.log10(support_min)
    if len(candidates) <= 8 and order_span <= 2.7:
        return candidates
    start_exp = int(np.ceil(np.log10(support_min)))
    end_exp = int(np.floor(np.log10(support_max)))
    powers = [float(10**exponent) for exponent in range(start_exp, end_exp + 1)]
    values = []
    if not powers or not np.isclose(powers[0], support_min):
        values.append(support_min)
    values.extend(powers)
    if not np.isclose(values[-1], support_max):
        values.append(support_max)
    if len(values) <= 8:
        return values
    keep_indexes = np.linspace(0, len(values) - 1, 8).round().astype(int)
    return [values[index] for index in sorted(set(keep_indexes))]


def support_size_legend_values(support_min: float, support_max: float) -> list[float]:
    if support_max <= 0:
        return []
    if support_min <= 0 or support_min >= support_max:
        return [nice_support_value(support_max)]
    targets = np.geomspace(support_min, support_max, num=6)
    values = [nice_support_value(float(value)) for value in targets]
    values[0] = support_min
    values[-1] = support_max
    deduped = []
    for value in values:
        if support_min <= value <= support_max and value not in deduped:
            deduped.append(float(value))
    return deduped


def support_tick_label(value: float) -> str:
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    if value >= 1:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return f"{value:.3g}"


def draw_support_size_scale(
    ax: plt.Axes,
    values: list[float],
    support_min: float,
    support_max: float,
    support_norm: colors.Normalize,
    cmap,
    support_label: str,
) -> None:
    ax.set_xscale("log")
    if support_min > 0 and support_max > support_min:
        log_pad = max(0.22, (np.log10(support_max) - np.log10(support_min)) * 0.06)
        ax.set_xlim(support_min / (10**log_pad), support_max * (10**log_pad))
    else:
        ax.set_xlim(support_min, support_max)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.scatter(
        values,
        [0.56] * len(values),
        s=rainfall_point_sizes(values, support_min, support_max),
        c=values,
        cmap=cmap,
        norm=support_norm,
        alpha=0.78,
        edgecolors="#1f2933",
        linewidths=0.45,
        clip_on=False,
        zorder=3,
    )
    ax.set_title("Point size", fontsize=9, fontweight="bold", pad=3)
    ax.set_xticks(values)
    ax.set_xticklabels([support_tick_label(value) for value in values], fontsize=8)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.tick_params(axis="x", length=3, width=0.7, pad=2)
    ax.set_xlabel(f"{support_label}\narea proportional to support", fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(True, axis="x", which="major", alpha=0.25, linewidth=0.6)


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
                facecolor=mitochondrial_feature_color(kind),
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
                facecolor=mitochondrial_feature_color(kind),
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


def sample_point_metadata(
    row: pd.Series,
    *,
    group: object,
    x_label: object,
    x_value: object,
    y_label: object,
    y_value: object,
) -> dict[str, object]:
    """Build the common metadata payload used by sample-point hover targets."""
    return {
        "sample": row.get("sample", ""),
        "group": group,
        "x_label": x_label,
        "x_value": x_value,
        "y_label": y_label,
        "y_value": y_value,
        "biological_replicate": row.get("biological_replicate", ""),
        "layout": row.get("layout", ""),
        "tissue": row.get("tissue", ""),
        "age": row.get("age", ""),
        "treatment": row.get("treatment", ""),
    }


def burden_plot(burden: pd.DataFrame, group_col: str, path: str, y: str, title: str, ylabel: str) -> None:
    if burden.empty or y not in burden.columns:
        empty(path, title, "No sample-level burden values are available")
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    pal = palette(burden, group_col)
    sample_points: list[dict[str, object]] = []
    if group_col and group_col in burden.columns:
        order = ordered_groups(burden, group_col)
        sns.stripplot(data=burden, x=group_col, y=y, hue=group_col, order=order, hue_order=order, palette=pal, size=7, jitter=0.18, ax=ax, legend=False)
        # One collection is emitted per group; retain the actual jittered
        # coordinates so the hit target follows the visible sample point.
        for group_index, group in enumerate(order):
            collection = ax.collections[group_index] if group_index < len(ax.collections) else None
            offsets = collection.get_offsets() if collection is not None else []
            rows = burden[burden[group_col].fillna("missing").astype(str) == str(group)].reset_index(drop=True)
            for row, offset in zip(rows.to_dict("records"), offsets):
                point = sample_point_metadata(
                    pd.Series(row),
                    group=group,
                    x_label=group_col,
                    x_value=group,
                    y_label=ylabel,
                    y_value=offset[1],
                )
                point.update({"x": offset[0], "y": offset[1]})
                sample_points.append(point)
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
    if sample_points:
        save_sample_points_interactive(fig, ax, path, sample_points)
    else:
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
    sample_points: list[dict[str, object]] = []
    collections = list(ax.collections)
    collection_index = 0
    for age in ages:
        for treatment in treatments:
            collection = collections[collection_index] if collection_index < len(collections) else None
            offsets = collection.get_offsets() if collection is not None else []
            rows = burden[
                (burden["age"].astype(str) == str(age))
                & (burden["treatment"].astype(str) == str(treatment))
            ].reset_index(drop=True)
            for row, offset in zip(rows.to_dict("records"), offsets):
                point = sample_point_metadata(
                    pd.Series(row),
                    group=treatment,
                    x_label="age",
                    x_value=age,
                    y_label=ylabel,
                    y_value=offset[1],
                )
                point.update({"x": offset[0], "y": offset[1]})
                sample_points.append(point)
            collection_index += 1
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
    save_sample_points_interactive(fig, ax, path, sample_points)


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


ORIGIN_OUTLINE_COLOR = "#00c8ff"
RANK_LABEL_FONT_SIZES = (6.0, 5.5, 5.0, 4.5, 4.0)
RANK_LABEL_MIN_MARKER_DIAMETER = 5.5
RANK_LABEL_FIT_FRACTION = 0.76


def apply_cluster_coordinates(reads: pd.DataFrame, clusters: pd.DataFrame | None, mt_length: int = 0) -> pd.DataFrame:
    """Replace read-level deletion coordinates with merged cluster representatives."""
    if reads.empty:
        return reads.copy()
    df = reads.copy()
    if clusters is not None and not clusters.empty:
        id_col = "junction_id" if "junction_id" in df.columns and "junction_id" in clusters.columns else ""
        if not id_col and "exact_deletion_id" in df.columns and "exact_deletion_id" in clusters.columns:
            id_col = "exact_deletion_id"
        required = {"left_breakpoint", "right_breakpoint", "deleted_size"}
        if id_col and required.issubset(clusters.columns):
            cluster_coords = clusters[[id_col, "left_breakpoint", "right_breakpoint", "deleted_size"]].drop_duplicates(id_col).copy()
            cluster_coords = cluster_coords.rename(
                columns={
                    "left_breakpoint": "_cluster_left_breakpoint",
                    "right_breakpoint": "_cluster_right_breakpoint",
                    "deleted_size": "_cluster_deleted_size",
                }
            )
            df = df.merge(cluster_coords, on=id_col, how="left")
            for col in ["left_breakpoint", "right_breakpoint", "deleted_size"]:
                cluster_col = f"_cluster_{col}"
                if cluster_col in df.columns:
                    df[col] = df[cluster_col].where(df[cluster_col].notna(), df[col])
            df = df.drop(columns=[col for col in df.columns if col.startswith("_cluster_")])
            # Cluster-level annotations are the canonical labels for a deletion.
            # Preserve any read-level value, while filling missing values from the
            # representative cluster row for interactive plot metadata.
            metadata_columns = [
                "affected_feature_label",
                "affected_features",
                "replication_arc_context",
                "major_arc_deleted_bp",
                "minor_arc_deleted_bp",
                "known_deletion_label",
                "deleted_interval",
                "wraps_origin",
            ]
            available_metadata = [column for column in metadata_columns if column in clusters.columns]
            if available_metadata:
                cluster_metadata = clusters[[id_col, *available_metadata]].drop_duplicates(id_col).copy()
                cluster_metadata = cluster_metadata.rename(
                    columns={column: f"_cluster_meta_{column}" for column in available_metadata}
                )
                df = df.merge(cluster_metadata, on=id_col, how="left")
                for column in available_metadata:
                    cluster_column = f"_cluster_meta_{column}"
                    if column not in df.columns:
                        df[column] = df[cluster_column]
                    else:
                        missing = df[column].isna() | df[column].astype(str).eq("")
                        df.loc[missing, column] = df.loc[missing, cluster_column]
                df = df.drop(columns=[f"_cluster_meta_{column}" for column in available_metadata])
    for col in ["left_breakpoint", "right_breakpoint", "deleted_size"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if mt_length and {"left_breakpoint", "right_breakpoint"}.issubset(df.columns):
        valid = df["left_breakpoint"].notna() & df["right_breakpoint"].notna()
        df.loc[valid, "deleted_size"] = [
            circular_distance(int(left), int(right), int(mt_length))
            for left, right in zip(df.loc[valid, "left_breakpoint"], df.loc[valid, "right_breakpoint"])
        ]
    return df


def location_support_colormap() -> colors.LinearSegmentedColormap:
    return colors.LinearSegmentedColormap.from_list(
        "support_dark_to_orange",
        ["#171126", "#3d1a6e", "#7c238e", "#bd3b78", "#ee654f", "#f6a15a"],
        N=256,
    )


def location_feature_display_name(name: str) -> str:
    if name == "mitochondrial_control_region":
        return "D-loop/control"
    return name.replace("_", " ")


def location_feature_class(row: pd.Series) -> str:
    kind = feature_kind(row)
    if kind == "protein-coding":
        return "protein_coding"
    if kind in {"rRNA", "tRNA"}:
        return kind
    feature_type = str(row.get("feature_type", "")).lower()
    name = feature_name(row)
    if feature_type == "region" or name in {"mitochondrial_control_region", "D-loop/control"}:
        return "region"
    return "other"


def read_yaml_safe(path: str) -> dict:
    if not path or not Path(path).exists():
        return {}
    with Path(path).open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def location_features(features: pd.DataFrame, config: dict | None = None) -> pd.DataFrame:
    if features.empty or not {"start", "end"}.issubset(features.columns):
        work = pd.DataFrame(columns=["start", "end", "name", "class"])
    else:
        work = features.copy()
        if "feature_type" in work.columns:
            feature_type = work["feature_type"].fillna("").astype(str).str.lower()
            keep = feature_type.isin({"gene", "region"})
            if keep.any():
                work = work[keep].copy()
        work = apply_feature_aliases(work, config or {})
        work["start"] = pd.to_numeric(work["start"], errors="coerce")
        work["end"] = pd.to_numeric(work["end"], errors="coerce")
        work["name"] = work.apply(lambda row: location_feature_display_name(feature_name(row)), axis=1)
        work["class"] = work.apply(location_feature_class, axis=1)
        work = work.dropna(subset=["start", "end"])
        work = work[work["end"] > work["start"]]
        work = work[["start", "end", "name", "class"]]
    region_rows = []
    for region in ((config or {}).get("analysis", {}) or {}).get("mt_regions", []) or []:
        try:
            region_rows.append(
                {
                    "start": int(region["start"]),
                    "end": int(region["end"]),
                    "name": location_feature_display_name(str(region["name"])),
                    "class": "region",
                }
            )
        except (KeyError, TypeError, ValueError):
            continue
    if region_rows:
        work = pd.concat([work, pd.DataFrame(region_rows)], ignore_index=True)
    return work.drop_duplicates().sort_values(["start", "end", "name"]).reset_index(drop=True)


def location_genome_length(mt_length: int, features: pd.DataFrame, df: pd.DataFrame | None = None) -> int:
    if mt_length > 0:
        return int(mt_length)
    values = []
    if not features.empty and "end" in features.columns:
        values.extend(pd.to_numeric(features["end"], errors="coerce").dropna().tolist())
    if df is not None and not df.empty:
        for col in ["left_breakpoint", "right_breakpoint"]:
            if col in df.columns:
                values.extend(pd.to_numeric(df[col], errors="coerce").dropna().tolist())
    return int(max(values)) if values else 1


def fitted_location_label_fontsize(ax: plt.Axes, label: str, visible_width_bp: float, x_span_bp: float, max_size: float, min_size: float) -> float | None:
    if visible_width_bp <= 0 or x_span_bp <= 0:
        return None
    fig_width_pt = ax.figure.get_size_inches()[0] * 72
    axis_width_pt = max(240.0, ax.get_position().width * fig_width_pt)
    available_pt = visible_width_bp / x_span_bp * axis_width_pt
    estimated_size = available_pt / max(1.0, len(label) * 0.60)
    fontsize = min(max_size, estimated_size)
    if fontsize < min_size:
        return None
    return max(min_size, fontsize)


def draw_location_feature_track(ax: plt.Axes, features: pd.DataFrame, genome_length: int, x_min: float = 1.0, x_max: float | None = None) -> None:
    if x_max is None:
        x_max = float(genome_length)
    rows = {"protein_coding": 3.25, "rRNA": 2.1, "tRNA": 1.0, "region": -0.15, "other": -1.0}
    labels = {"protein_coding": "protein-coding", "rRNA": "rRNA", "tRNA": "tRNA", "region": "regions"}
    x_span = max(1.0, float(x_max) - float(x_min))
    dloop_segments: list[tuple[float, float]] = []
    labeled_names: set[str] = set()
    for _, feat in features.iterrows():
        cls = feat["class"] if feat["class"] in rows else "other"
        y = rows[cls]
        start = float(feat["start"])
        end = float(feat["end"])
        if end <= x_min or start >= x_max:
            continue
        visible_start = max(start, float(x_min))
        visible_end = min(end, float(x_max))
        visible_width = max(0.0, visible_end - visible_start)
        ax.add_patch(
            Rectangle(
                (start, y - 0.28),
                max(1.0, end - start + 1),
                0.56,
                color=mitochondrial_feature_color(cls),
                alpha=0.9,
            )
        )
        name = str(feat["name"])
        if cls == "region" and name == "D-loop/control" and visible_width > 20:
            dloop_segments.append((visible_start, visible_end))
            continue
        if cls not in {"protein_coding", "rRNA", "region"} or visible_width <= 45 or name in labeled_names:
            continue
        text_x = (visible_start + visible_end) / 2
        ha = "center"
        label_width = visible_width
        max_font = 6.2 if cls == "protein_coding" else 7.0
        min_font = 3.4 if cls == "protein_coding" else 4.0
        if cls == "region":
            label_width = max(visible_width, 380.0)
            if visible_start <= x_min + 75:
                text_x = x_min + 45
                ha = "left"
            elif visible_end >= x_max - 75:
                text_x = x_max - 45
                ha = "right"
        fontsize = fitted_location_label_fontsize(ax, name, label_width, x_span, max_font, min_font)
        if fontsize is None:
            continue
        ax.text(text_x, y, name, ha=ha, va="center", fontsize=fontsize, color="#111827", clip_on=True)
        labeled_names.add(name)
    for start, end in dloop_segments:
        if end <= start:
            continue
        if start <= x_min + 100:
            text_x = x_min + 45
            ha = "left"
        elif end >= x_max - 100:
            text_x = x_max - 45
            ha = "right"
        else:
            text_x = (start + end) / 2
            ha = "center"
        fontsize = fitted_location_label_fontsize(ax, "D-loop/control", max(end - start, 1300.0), x_span, 7.0, 3.6) or 3.6
        ax.text(text_x, rows["region"], "D-loop/control", ha=ha, va="center", fontsize=fontsize, color="#111827", clip_on=True)
    ax.set_ylim(-1.0, 4.18)
    ax.set_yticks([rows[key] for key in ["protein_coding", "rRNA", "tRNA", "region"]])
    ax.set_yticklabels([labels[key] for key in ["protein_coding", "rRNA", "tRNA", "region"]], fontsize=8)
    ax.set_xlim(x_min, x_max)
    ax.grid(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)


def circular_deleted_interval_midpoint(df: pd.DataFrame, genome_length: int) -> pd.Series:
    left = pd.to_numeric(df["left_breakpoint"], errors="coerce")
    right = pd.to_numeric(df["right_breakpoint"], errors="coerce")
    crosses = right < left
    adjusted_right = right.where(~crosses, right + genome_length)
    midpoint = (left + adjusted_right) / 2.0
    return ((midpoint - 1) % genome_length) + 1


def origin_outline_linewidths(marker_sizes: np.ndarray) -> tuple[float, float]:
    if len(marker_sizes) == 0:
        return 0.45, 0.18
    size = float(np.nanmedian(marker_sizes))
    if size < 6:
        return 0.42, 0.18
    if size < 18:
        return 0.55, 0.24
    if size < 45:
        return 0.75, 0.32
    if size < 120:
        return 0.95, 0.42
    if size < 280:
        return 1.20, 0.55
    return 1.45, 0.70


def support_ordered_groups(df: pd.DataFrame, support_col: str = "_plot_support") -> list[pd.DataFrame]:
    if df.empty:
        return []
    ordered = df.sort_values(support_col, ascending=True, kind="mergesort")
    return [group.copy() for _, group in ordered.groupby(support_col, sort=True)]


def assign_group_support_ranks(df: pd.DataFrame) -> pd.DataFrame:
    """Assign stable, unique support ranks independently within each plot group."""
    if df.empty:
        out = df.copy()
        out["_support_rank"] = pd.Series(dtype="Int64")
        return out
    out = df.copy()
    tie_columns = [col for col in ["exact_deletion_id", "left_breakpoint", "right_breakpoint", "deleted_size"] if col in out.columns]
    sort_columns = ["_plot_group", "_plot_support", *tie_columns]
    ascending = [True, False, *([True] * len(tie_columns))]
    out = out.sort_values(sort_columns, ascending=ascending, kind="mergesort")
    out["_support_rank"] = out.groupby("_plot_group", sort=False).cumcount() + 1
    return out


def rank_label_font_size(marker_area: float, rank: int) -> float | None:
    """Return the largest allowed font size whose worst-case digits fit a marker."""
    if not np.isfinite(marker_area) or marker_area <= 0 or int(rank) < 1:
        return None
    marker_diameter = float(np.sqrt(marker_area))
    if marker_diameter < RANK_LABEL_MIN_MARKER_DIAMETER:
        return None
    digit_template = "8" * len(str(int(rank)))
    available = marker_diameter * RANK_LABEL_FIT_FRACTION
    font = FontProperties(family="DejaVu Sans", weight="bold")
    for font_size in RANK_LABEL_FONT_SIZES:
        bounds = TextPath((0, 0), digit_template, size=font_size, prop=font).get_extents()
        if bounds.width <= available and bounds.height <= available:
            return font_size
    return None


def rank_label_color(support: float, support_norm: colors.Normalize, cmap) -> str:
    red, green, blue, _ = cmap(support_norm(float(support)))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#111827" if luminance >= 0.53 else "#ffffff"


def rank_label_boxes_overlap(first: tuple[float, float, float], second: tuple[float, float, float]) -> bool:
    """Return true when two display-coordinate label clearance circles overlap."""
    x1, y1, radius1 = first
    x2, y2, radius2 = second
    return float(np.hypot(x1 - x2, y1 - y2)) < 0.72 * (radius1 + radius2)


def adjusted_breakpoint_ticks(genome_length: int, y_max: float) -> tuple[list[float], list[str]]:
    base_ticks = [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, genome_length]
    ticks = [float(tick) for tick in base_ticks if tick <= y_max]
    labels = [f"{int(tick):,}" if tick != genome_length else f"{genome_length:,} / 0" for tick in ticks]
    for tick in [2000, 4000, 6000, 8000]:
        adjusted = genome_length + tick
        if adjusted <= y_max:
            ticks.append(float(adjusted))
            labels.append(f"{tick:,}")
    return ticks, labels


def circular_rolling_mean(values: np.ndarray, window: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if len(arr) == 0 or window <= 1:
        return arr.copy()
    if window % 2 == 0:
        window += 1
    radius = window // 2
    total = np.zeros_like(arr, dtype=float)
    for shift in range(-radius, radius + 1):
        total += np.roll(arr, shift)
    return total / float(window)


def pooled_endpoint_density(display_grouped: pd.DataFrame, genome_length: int, bin_size: int, smooth_bins: int) -> pd.DataFrame:
    if display_grouped.empty or genome_length <= 0:
        return pd.DataFrame(columns=["bin_index", "endpoint_count", "left_endpoint_count", "right_endpoint_count", "summed_support", "left_support", "right_support", "raw_supporting_reads", "left_raw_supporting_reads", "right_raw_supporting_reads", "bin_start", "bin_end", "bin_midpoint", "smoothed_summed_support", "smoothed_endpoint_count"])
    bin_size = max(1, int(bin_size))
    def endpoint_rows(breakpoint_column: str, side: str) -> pd.DataFrame:
        frame = display_grouped[[breakpoint_column, "_plot_support"]].rename(columns={breakpoint_column: "coordinate"}).copy()
        if "supporting_reads" in display_grouped.columns:
            frame["raw_supporting_reads"] = pd.to_numeric(display_grouped["supporting_reads"], errors="coerce").to_numpy()
        else:
            frame["raw_supporting_reads"] = np.nan
        frame["endpoint_side"] = side
        return frame

    endpoints = pd.concat(
        [endpoint_rows("left_breakpoint", "left"), endpoint_rows("right_breakpoint", "right")],
        ignore_index=True,
    )
    endpoints["coordinate"] = pd.to_numeric(endpoints["coordinate"], errors="coerce")
    endpoints["support"] = pd.to_numeric(endpoints["_plot_support"], errors="coerce")
    endpoints["raw_supporting_reads"] = pd.to_numeric(endpoints["raw_supporting_reads"], errors="coerce")
    endpoints = endpoints.dropna(subset=["coordinate", "support"])
    endpoints = endpoints[(endpoints["coordinate"] >= 1) & (endpoints["coordinate"] <= genome_length)]
    endpoints["bin_index"] = ((endpoints["coordinate"] - 1) // bin_size).astype(int)
    n_bins = int(np.ceil(genome_length / bin_size))
    base = pd.DataFrame({"bin_index": np.arange(n_bins, dtype=int)})
    side = (
        endpoints.groupby(["bin_index", "endpoint_side"], as_index=False)
        .agg(endpoint_count=("coordinate", "size"), support=("support", "sum"), raw_supporting_reads=("raw_supporting_reads", "sum"))
    )
    counts = side.pivot_table(index="bin_index", columns="endpoint_side", values="endpoint_count", aggfunc="sum", fill_value=0)
    support = side.pivot_table(index="bin_index", columns="endpoint_side", values="support", aggfunc="sum", fill_value=0)
    agg = base.copy()
    agg["left_endpoint_count"] = counts.get("left", pd.Series(0, index=counts.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    agg["right_endpoint_count"] = counts.get("right", pd.Series(0, index=counts.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    agg["left_support"] = support.get("left", pd.Series(0, index=support.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    agg["right_support"] = support.get("right", pd.Series(0, index=support.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    raw_support = side.pivot_table(index="bin_index", columns="endpoint_side", values="raw_supporting_reads", aggfunc="sum", fill_value=0)
    agg["left_raw_supporting_reads"] = raw_support.get("left", pd.Series(0, index=raw_support.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    agg["right_raw_supporting_reads"] = raw_support.get("right", pd.Series(0, index=raw_support.index)).reindex(agg["bin_index"], fill_value=0).to_numpy()
    agg["endpoint_count"] = agg["left_endpoint_count"] + agg["right_endpoint_count"]
    agg["summed_support"] = agg["left_support"] + agg["right_support"]
    agg["raw_supporting_reads"] = agg["left_raw_supporting_reads"] + agg["right_raw_supporting_reads"]
    out = agg
    out["bin_start"] = out["bin_index"] * bin_size + 1
    out["bin_end"] = np.minimum(out["bin_start"] + bin_size - 1, genome_length)
    out["bin_midpoint"] = (out["bin_start"] + out["bin_end"]) / 2.0
    out["smoothed_summed_support"] = circular_rolling_mean(out["summed_support"].to_numpy(dtype=float), max(1, int(smooth_bins)))
    out["smoothed_endpoint_count"] = circular_rolling_mean(out["endpoint_count"].to_numpy(dtype=float), max(1, int(smooth_bins)))
    return out


def endpoint_density_hotspots(
    density: pd.DataFrame,
    genome_length: int,
    max_labels: int = 8,
    min_spacing_bp: int | None = None,
) -> list[dict[str, float]]:
    if density.empty or "smoothed_summed_support" not in density.columns:
        return []
    work = density.copy()
    y = pd.to_numeric(work["smoothed_summed_support"], errors="coerce").fillna(0).to_numpy(dtype=float)
    if len(y) < 3 or np.nanmax(y) <= 0:
        return []
    left = np.roll(y, 1)
    right = np.roll(y, -1)
    candidates = work[(y >= left) & (y >= right) & (work["smoothed_summed_support"] > 0)].copy()
    candidates = candidates.sort_values("smoothed_summed_support", ascending=False)
    selected: list[dict[str, float]] = []
    min_spacing = max(1, int(min_spacing_bp if min_spacing_bp is not None else genome_length * 0.075))
    for _, row in candidates.iterrows():
        coord = float(row["bin_midpoint"])
        if any(min(abs(coord - item["coord"]), genome_length - abs(coord - item["coord"])) < min_spacing for item in selected):
            continue
        selected.append({"coord": coord, "height": float(row["smoothed_summed_support"])})
        if len(selected) >= max_labels:
            break
    return selected


def endpoint_label_levels(hotspots: list[dict[str, float]], genome_length: int, max_levels: int = 4) -> list[tuple[dict[str, float], int]]:
    if not hotspots:
        return []
    min_spacing = max(850, int(genome_length * 0.05))
    last_coord_by_level = [-float("inf")] * max_levels
    assigned: list[tuple[dict[str, float], int]] = []
    for hotspot in sorted(hotspots, key=lambda item: item["coord"]):
        coord = float(hotspot["coord"])
        chosen = None
        for level in range(max_levels):
            if coord - last_coord_by_level[level] >= min_spacing:
                chosen = level
                break
        if chosen is None:
            chosen = min(range(max_levels), key=lambda level: last_coord_by_level[level])
        last_coord_by_level[chosen] = coord
        assigned.append((hotspot, chosen))
    return assigned


def endpoint_density_figure(
    display_grouped: pd.DataFrame,
    features: pd.DataFrame,
    config: dict,
    mt_length: int,
    title: str,
    support_label: str,
    bin_size: int = 50,
    smooth_bins: int = 7,
    capped: bool = False,
) -> plt.Figure:
    plot_features = location_features(features, config)
    genome_length = location_genome_length(mt_length, plot_features, display_grouped)
    fig = plt.figure(figsize=(15.2, 6.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[4.65, 1.05], hspace=0.04)
    ax = fig.add_subplot(grid[0, 0])
    feature_ax = fig.add_subplot(grid[1, 0], sharex=ax)
    if display_grouped.empty:
        ax.text(0.5, 0.5, "No exact deletions meet the endpoint-density display threshold", ha="center", va="center", wrap=True, transform=ax.transAxes)
        ax.set_title(title, pad=15)
        draw_location_feature_track(feature_ax, plot_features, genome_length, x_min=1, x_max=genome_length)
        feature_ax.set_xlabel("Mitochondrial genome coordinate (bp)")
        return fig
    density = pooled_endpoint_density(display_grouped, genome_length, bin_size, smooth_bins)
    if density.empty:
        ax.text(0.5, 0.5, "No breakpoint endpoints are available for this plot", ha="center", va="center", wrap=True, transform=ax.transAxes)
        ax.set_title(title, pad=15)
        draw_location_feature_track(feature_ax, plot_features, genome_length, x_min=1, x_max=genome_length)
        feature_ax.set_xlabel("Mitochondrial genome coordinate (bp)")
        return fig
    # Retain the exact plotted aggregation for the report's interactive SVG
    # sidecar; the static PDF and its smoothing remain unchanged.
    fig._endpoint_density = density.copy()
    fig._endpoint_density_bin_size = int(bin_size)
    fig._endpoint_density_support_label = support_label
    fig._endpoint_density_capped = bool(capped)
    x = density["bin_midpoint"].to_numpy(dtype=float)
    raw_left = density["left_support"].to_numpy(dtype=float)
    raw_right = density["right_support"].to_numpy(dtype=float)
    raw = density["summed_support"].to_numpy(dtype=float)
    y = density["smoothed_summed_support"].to_numpy(dtype=float)
    width = max(1, int(bin_size)) * 0.88
    ax.bar(x, raw_left, width=width, color="#4E79A7", alpha=0.42, edgecolor="none", label=f"left endpoint support, {int(bin_size)} bp bins", zorder=1)
    ax.bar(x, raw_right, width=width, bottom=raw_left, color="#F28E2B", alpha=0.42, edgecolor="none", label=f"right endpoint support, {int(bin_size)} bp bins", zorder=1)
    ax.fill_between(x, y, color="#111111", alpha=0.08, linewidth=0, zorder=2)
    ax.plot(x, y, color="#222222", linewidth=1.6, label="circular-smoothed pooled support", zorder=3)
    finite_heights = np.r_[raw, y]
    finite_heights = finite_heights[np.isfinite(finite_heights)]
    full_ymax = float(np.nanmax(finite_heights)) if len(finite_heights) else 0.0
    if full_ymax <= 0:
        full_ymax = 1.0
    cap_note = ""
    if capped:
        nonzero = finite_heights[finite_heights > 0]
        cap = float(np.nanpercentile(nonzero, 96)) if len(nonzero) else full_ymax
        ymax = min(full_ymax, cap) if cap > 0 else full_ymax
        ax.set_ylim(0, ymax * 1.12)
        cap_note = f"; y-axis capped at {support_tick_label(ymax)} to show secondary hotspots"
    else:
        ax.set_ylim(0, full_ymax * 1.12)
    ax.set_xlim(1, genome_length)
    ax.set_title(title, pad=15)
    normalized = "per million" in support_label.lower()
    ax.text(
        0.5,
        1.01,
        f"{'support per million usable reads; ' if normalized else ''}left and right endpoints pooled; {int(bin_size)} bp bins; {int(smooth_bins) * int(bin_size)} bp circular smoothing window{cap_note}",
        ha="center",
        va="bottom",
        transform=ax.transAxes,
        fontsize=9,
        color="#4b5563",
    )
    ax.set_ylabel(f"Summed {support_label.lower()} at endpoints")
    ax.grid(axis="both", color="#d9dee7", linewidth=0.65, alpha=0.55)
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    label_ceiling = ax.get_ylim()[1]
    smoothing_window_bp = max(1, int(smooth_bins)) * max(1, int(bin_size))
    label_spacing_bp = max(2 * smoothing_window_bp, int(genome_length * 0.075))
    hotspots = endpoint_density_hotspots(density, genome_length, min_spacing_bp=label_spacing_bp)
    for hotspot, label_level in endpoint_label_levels(hotspots, genome_length):
        coord = hotspot["coord"]
        height = min(hotspot["height"], label_ceiling * 0.80)
        half_window = max(1, int(smooth_bins) * int(bin_size) / 2)
        start = max(1, int(round(coord - half_window)))
        end = min(genome_length, int(round(coord + half_window)))
        ax.annotate(
            f"{start:,}-{end:,} bp",
            xy=(coord, height),
            xytext=(0, 10 + label_level * 12),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#1f2937",
            arrowprops={"arrowstyle": "-", "color": "#6b7280", "linewidth": 0.6, "alpha": 0.8},
            clip_on=True,
        )
    draw_location_feature_track(feature_ax, plot_features, genome_length, x_min=1, x_max=genome_length)
    feature_ax.set_xlabel("Mitochondrial genome coordinate (bp)")
    return fig


def endpoint_density_pages(
    display_grouped: pd.DataFrame,
    groups: list[str],
    features: pd.DataFrame,
    config: dict,
    mt_length: int,
    path: str,
    title_prefix: str,
    support_label: str,
    bin_size: int = 50,
    smooth_bins: int = 7,
    capped: bool = False,
) -> None:
    clear_location_sidecars(path)
    if display_grouped.empty:
        empty(path, title_prefix, "No exact deletions meet the endpoint-density display threshold")
        return
    figures: list[plt.Figure] = []
    for group in groups:
        sub = display_grouped[display_grouped["_plot_group"] == group].copy()
        fig = endpoint_density_figure(
            sub,
            features,
            config,
            mt_length,
            f"{title_prefix}: {group}",
            support_label,
            bin_size=bin_size,
            smooth_bins=smooth_bins,
            capped=capped,
        )
        save_location_sidecar(fig, path, group)
        density = getattr(fig, "_endpoint_density", pd.DataFrame())
        if not density.empty:
            save_endpoint_density_interactive_sidecar(
                fig,
                ax=fig.axes[0],
                path=path,
                group=group,
                density=density,
                support_label=support_label,
            )
        figures.append(fig)
    write_location_plot_pages(figures, path)
    for fig in figures:
        plt.close(fig)


def write_location_plot_pages(figures: list[plt.Figure], path: str) -> None:
    ensure_parent(path)
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    if figures:
        figures[0].savefig(Path(path).with_suffix(".svg"), bbox_inches="tight")


def save_location_sidecar(fig: plt.Figure, path: str, group: str) -> None:
    sidecar_pdf = Path(path).with_name(f"{Path(path).stem}__{safe_filename(group)}.pdf")
    ensure_parent(sidecar_pdf)
    fig.savefig(sidecar_pdf, bbox_inches="tight")
    fig.savefig(sidecar_pdf.with_suffix(".svg"), bbox_inches="tight")


def svg_attribute_value(value: object) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return html.escape(str(value), quote=True)


def save_rainfall_interactive_sidecar(
    fig: plt.Figure,
    ax: plt.Axes,
    path: str,
    group: str,
    data: pd.DataFrame,
    x_col: str,
    support_min: float,
    support_max: float,
    support_norm: colors.Normalize,
    cmap,
    support_label: str = "Deletion support",
) -> None:
    """Write a full-call SVG with point metadata for the report controls.

    The static Matplotlib point collections are hidden in this sidecar and replaced
    by individually addressable circles. This keeps the PDF unchanged while making
    support filtering and hover inspection independent of the static display.
    """
    if data.empty:
        return
    sidecar = Path(path).with_name(f"{Path(path).stem}__{safe_filename(group)}__interactive.svg")
    ensure_parent(sidecar)
    fig.canvas.draw()
    canvas_width, canvas_height = fig.canvas.get_width_height()
    temporary = sidecar.with_suffix(".tmp.svg")
    fig.savefig(temporary, format="svg")
    svg = temporary.read_text(encoding="utf-8")
    temporary.unlink(missing_ok=True)
    viewbox_match = re.search(r'<svg\b[^>]*\bviewBox="0 0 ([0-9.eE+-]+) ([0-9.eE+-]+)"', svg)
    if not viewbox_match:
        return
    viewbox_width = float(viewbox_match.group(1))
    viewbox_height = float(viewbox_match.group(2))
    scale_x = viewbox_width / float(canvas_width)
    scale_y = viewbox_height / float(canvas_height)
    circles = []
    sort_columns = ["_plot_support", *[col for col in ["exact_deletion_id", "left_breakpoint", "right_breakpoint"] if col in data.columns]]
    work = data.sort_values(sort_columns, ascending=[True] * len(sort_columns), kind="mergesort")
    marker_areas = rainfall_point_sizes(work["_plot_support"], support_min, support_max)
    for (_, row), marker_area in zip(work.iterrows(), marker_areas):
        x_value = float(row[x_col])
        y_value = float(row["deleted_size"])
        display_x, display_y = ax.transData.transform((x_value, y_value))
        color = colors.to_hex(cmap(support_norm(float(row["_plot_support"]))), keep_alpha=False)
        radius = max(1.7, float(np.sqrt(float(marker_area) / np.pi))) * scale_x
        attrs = {
            "class": "rainfall-point",
            "cx": f"{display_x * scale_x:.3f}",
            "cy": f"{viewbox_height - display_y * scale_y:.3f}",
            "r": f"{radius:.3f}",
            "fill": color,
            "fill-opacity": "0.86",
            "stroke": ORIGIN_OUTLINE_COLOR if bool(row.get("crosses_origin", False)) else "#17202a",
            "stroke-width": "1.15" if bool(row.get("crosses_origin", False)) else "0.45",
            "data-group": group,
            "data-exact-deletion-id": row.get("exact_deletion_id", ""),
            "data-left-breakpoint": row.get("left_breakpoint", ""),
            "data-right-breakpoint": row.get("right_breakpoint", ""),
            "data-deleted-size": row.get("deleted_size", ""),
            "data-support": row.get("_plot_support", ""),
            "data-support-label": support_label,
            "data-supporting-reads": row.get("supporting_reads", ""),
            "data-crosses-origin": "yes" if bool(row.get("crosses_origin", False)) else "no",
            "data-affected-features": row.get("affected_feature_label", row.get("affected_features", "")),
            "data-arc-context": row.get("replication_arc_context", ""),
            "data-major-arc-bp": row.get("major_arc_deleted_bp", ""),
            "data-minor-arc-bp": row.get("minor_arc_deleted_bp", ""),
            "data-known-deletion": row.get("known_deletion_label", ""),
        }
        rendered = " ".join(f'{key}="{svg_attribute_value(value)}"' for key, value in attrs.items())
        circles.append(f"<circle {rendered}/>")
    root_match = re.search(r"<svg\b[^>]*>", svg)
    if not root_match:
        return
    metadata = (
        f' data-plot-type="rainfall" data-group="{svg_attribute_value(group)}"'
        f' data-call-count="{len(work)}" data-support-min="{svg_attribute_value(support_min)}"'
        f' data-support-max="{svg_attribute_value(support_max)}" data-support-label="{svg_attribute_value(support_label)}"'
    )
    root_end = root_match.end() - 1
    svg = svg[:root_end] + metadata + svg[root_end:]
    style = (
        "<style>g[id^=\"rainfall-static-points-\"]{display:none;}.rainfall-point{cursor:help;}"
        ".rainfall-point:hover{stroke:#172b4d;stroke-width:1.8;fill-opacity:1;}</style>"
    )
    svg = svg[: root_end + len(metadata) + 1] + style + svg[root_end + len(metadata) + 1 :]
    overlay = '<g id="rainfall-interactive-points">' + "".join(circles) + "</g>"
    svg = svg.replace("</svg>", overlay + "</svg>")
    sidecar.write_text(svg, encoding="utf-8")


def save_endpoint_density_interactive_sidecar(
    fig: plt.Figure,
    ax: plt.Axes,
    path: str,
    group: str,
    density: pd.DataFrame,
    support_label: str,
) -> None:
    """Write transparent bin hit-targets over the static density plot SVG."""
    if density.empty:
        return
    sidecar = Path(path).with_name(f"{Path(path).stem}__{safe_filename(group)}__interactive.svg")
    ensure_parent(sidecar)
    fig.canvas.draw()
    canvas_width, canvas_height = fig.canvas.get_width_height()
    temporary = sidecar.with_suffix(".tmp.svg")
    fig.savefig(temporary, format="svg")
    svg = temporary.read_text(encoding="utf-8")
    temporary.unlink(missing_ok=True)
    viewbox_match = re.search(r'<svg\b[^>]*\bviewBox="0 0 ([0-9.eE+-]+) ([0-9.eE+-]+)"', svg)
    root_match = re.search(r"<svg\b[^>]*>", svg)
    if not viewbox_match or not root_match:
        return
    viewbox_width = float(viewbox_match.group(1))
    viewbox_height = float(viewbox_match.group(2))
    scale_x = viewbox_width / float(canvas_width)
    scale_y = viewbox_height / float(canvas_height)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_min_px = ax.transData.transform((x_min, y_min))[0]
    x_max_px = ax.transData.transform((x_max, y_min))[0]
    y_bottom_px = ax.transData.transform((x_min, y_min))[1]
    y_top_px = ax.transData.transform((x_min, y_max))[1]
    y_svg = viewbox_height - max(y_bottom_px, y_top_px) * scale_y
    height_svg = abs(y_top_px - y_bottom_px) * scale_y
    rectangles = []
    for _, row in density.iterrows():
        start = float(row["bin_start"])
        end = float(row["bin_end"])
        left_px = ax.transData.transform((max(x_min, start), y_min))[0]
        right_px = ax.transData.transform((min(x_max, end + 1), y_min))[0]
        left_px = max(x_min_px, min(x_max_px, left_px))
        right_px = max(x_min_px, min(x_max_px, right_px))
        if right_px <= left_px:
            continue
        attrs = {
            "class": "endpoint-density-bin",
            "x": f"{left_px * scale_x:.3f}",
            "y": f"{y_svg:.3f}",
            "width": f"{(right_px - left_px) * scale_x:.3f}",
            "height": f"{height_svg:.3f}",
            "fill": "#285f8f",
            "fill-opacity": "0",
            "data-group": group,
            "data-bin-start": row.get("bin_start", ""),
            "data-bin-end": row.get("bin_end", ""),
            "data-bin-midpoint": row.get("bin_midpoint", ""),
            "data-left-endpoint-count": row.get("left_endpoint_count", ""),
            "data-right-endpoint-count": row.get("right_endpoint_count", ""),
            "data-endpoint-count": row.get("endpoint_count", ""),
            "data-left-support": row.get("left_support", ""),
            "data-right-support": row.get("right_support", ""),
            "data-summed-support": row.get("summed_support", ""),
            "data-left-raw-supporting-reads": row.get("left_raw_supporting_reads", ""),
            "data-right-raw-supporting-reads": row.get("right_raw_supporting_reads", ""),
            "data-raw-supporting-reads": row.get("raw_supporting_reads", ""),
            "data-smoothed-support": row.get("smoothed_summed_support", ""),
            "data-smoothed-endpoint-count": row.get("smoothed_endpoint_count", ""),
            "data-support-label": support_label,
        }
        rendered = " ".join(f'{key}="{svg_attribute_value(value)}"' for key, value in attrs.items())
        rectangles.append(f"<rect {rendered}/>")
    metadata = (
        f' data-plot-type="endpoint-density" data-group="{svg_attribute_value(group)}"'
        f' data-bin-count="{len(rectangles)}" data-support-label="{svg_attribute_value(support_label)}"'
    )
    root_end = root_match.end() - 1
    svg = svg[:root_end] + metadata + svg[root_end:]
    style = '<style>.endpoint-density-bin{cursor:help;}.endpoint-density-bin:hover{fill-opacity:.08;}</style>'
    svg = svg[: root_end + len(metadata) + 1] + style + svg[root_end + len(metadata) + 1 :]
    overlay = '<g id="endpoint-density-interactive-bins">' + "".join(rectangles) + "</g>"
    svg = svg.replace("</svg>", overlay + "</svg>")
    sidecar.write_text(svg, encoding="utf-8")


def save_breakpoint_pair_interactive_sidecar(
    fig: plt.Figure,
    ax: plt.Axes,
    path: str,
    group: str,
    pairs: pd.DataFrame,
    support_min: float,
    support_max: float,
    support_label: str = "Deletion support",
) -> None:
    """Write transparent point hit targets over the breakpoint-pair map SVG."""
    if pairs.empty:
        return
    sidecar = Path(path).with_name(f"{Path(path).stem}__{safe_filename(group)}__interactive.svg")
    ensure_parent(sidecar)
    fig.canvas.draw()
    canvas_width, canvas_height = fig.canvas.get_width_height()
    temporary = sidecar.with_suffix(".tmp.svg")
    # Keep the SVG in the figure's native coordinate system so transparent hit
    # targets use the same transform as the plotted points.
    fig.savefig(temporary, format="svg")
    svg = temporary.read_text(encoding="utf-8")
    temporary.unlink(missing_ok=True)
    viewbox_match = re.search(
        r'<svg\b[^>]*\bviewBox="([0-9.eE+-]+) ([0-9.eE+-]+) '
        r'([0-9.eE+-]+) ([0-9.eE+-]+)"',
        svg,
    )
    root_match = re.search(r"<svg\b[^>]*>", svg)
    if not viewbox_match or not root_match:
        return
    viewbox_width = float(viewbox_match.group(3))
    viewbox_height = float(viewbox_match.group(4))
    scale_x = viewbox_width / float(canvas_width)
    scale_y = viewbox_height / float(canvas_height)
    circles = []
    work = pairs.sort_values("_plot_support", ascending=True, kind="mergesort")
    marker_areas = rainfall_point_sizes(work["_plot_support"], support_min, support_max)
    for (_, row), marker_area in zip(work.iterrows(), marker_areas):
        x_value = float(row["left_breakpoint"])
        y_value = float(row["adjusted_right_breakpoint"])
        display_x, display_y = ax.transData.transform((x_value, y_value))
        affected = row.get("affected_feature_label", row.get("affected_features", ""))
        attrs = {
            "class": "breakpoint-pair-point",
            "cx": f"{display_x * scale_x:.3f}",
            "cy": f"{viewbox_height - display_y * scale_y:.3f}",
            "r": f"{max(4.0, np.sqrt(float(marker_area) / np.pi)) * scale_x:.3f}",
            "fill": "#285f8f",
            "fill-opacity": "0",
            "stroke": "transparent",
            "data-group": group,
            "data-exact-deletion-id": row.get("exact_deletion_id", ""),
            "data-left-breakpoint": row.get("left_breakpoint", ""),
            "data-right-breakpoint": row.get("right_breakpoint", ""),
            "data-deleted-size": row.get("deleted_size", ""),
            "data-adjusted-right-breakpoint": row.get("adjusted_right_breakpoint", ""),
            "data-support": row.get("_plot_support", ""),
            "data-supporting-observations": row.get("supporting_reads", ""),
            "data-support-label": support_label,
            "data-pair-count": row.get("pair_count", ""),
            "data-rank": row.get("_support_rank", ""),
            "data-crosses-origin": "yes" if bool(row.get("crosses_origin", False)) else "no",
            "data-affected-features": "" if pd.isna(affected) else affected,
            "data-arc-context": row.get("replication_arc_context", ""),
            "data-major-arc-bp": row.get("major_arc_deleted_bp", ""),
            "data-minor-arc-bp": row.get("minor_arc_deleted_bp", ""),
        }
        rendered = " ".join(f'{key}="{svg_attribute_value(value)}"' for key, value in attrs.items())
        circles.append(f"<circle {rendered}/>")
    metadata = (
        f' data-plot-type="breakpoint-pair-map" data-group="{svg_attribute_value(group)}"'
        f' data-point-count="{len(circles)}" data-support-label="{svg_attribute_value(support_label)}"'
    )
    root_end = root_match.end() - 1
    svg = svg[:root_end] + metadata + svg[root_end:]
    style = '<style>.breakpoint-pair-point{cursor:help;}.breakpoint-pair-point:hover{stroke:#172b4d;stroke-width:2;fill-opacity:.16;}</style>'
    svg = svg[: root_end + len(metadata) + 1] + style + svg[root_end + len(metadata) + 1 :]
    svg = svg.replace("</svg>", '<g id="breakpoint-pair-interactive-points">' + "".join(circles) + "</g></svg>")
    sidecar.write_text(svg, encoding="utf-8")


def clear_location_sidecars(path: str) -> None:
    for old in Path(path).parent.glob(f"{Path(path).stem}__*.pdf"):
        old.unlink()
    for old in Path(path).parent.glob(f"{Path(path).stem}__*.svg"):
        old.unlink()


def add_location_legend(
    fig: plt.Figure,
    legend_ax: plt.Axes,
    scatter,
    legend_values: list[float],
    support_min: float,
    support_max: float,
    support_label: str,
    show_origin_outline: bool,
    origin_label: str,
    note: str | None = None,
    show_rank_note: bool = True,
) -> None:
    legend_ax.set_axis_off()
    legend_ax.text(0.5, 0.82, "Point color", ha="center", va="center", fontsize=9, fontweight="bold", transform=legend_ax.transAxes)
    legend_ax.text(0.5, 0.76, f"{support_label}\nlog color scale", ha="center", va="center", fontsize=8, transform=legend_ax.transAxes)
    cax = legend_ax.inset_axes([0.16, 0.68, 0.72, 0.045])
    cb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
    cb.set_ticks(legend_values)
    cb.ax.set_xticklabels([support_tick_label(float(value)) for value in legend_values], fontsize=7)
    cb.ax.xaxis.set_minor_locator(ticker.NullLocator())
    cb.ax.tick_params(axis="x", length=2, pad=1)
    size_ax = legend_ax.inset_axes([0.13, 0.38, 0.80, 0.18])
    draw_support_size_scale(size_ax, support_size_legend_values(support_min, support_max), support_min, support_max, scatter.norm, scatter.cmap, support_label)
    if show_origin_outline:
        legend_ax.scatter([0.30], [0.21], s=120, facecolors="none", edgecolors=ORIGIN_OUTLINE_COLOR, linewidths=1.5, transform=legend_ax.transAxes, clip_on=False)
        legend_ax.text(0.44, 0.21, origin_label, va="center", fontsize=8, transform=legend_ax.transAxes)
    rank_note = "Numbers identify the same support-ranked calls across the breakpoint-pair view." if show_rank_note else "Interactive HTML report: hover a point for exact deletion details."
    if note:
        rank_note += f"\n{note}"
    legend_ax.text(0.5, 0.075, rank_note, ha="center", va="center", fontsize=7.1, linespacing=1.18, transform=legend_ax.transAxes)


def prepare_location_plot_data(
    reads: pd.DataFrame,
    samples: pd.DataFrame,
    group_col: str,
    clusters: pd.DataFrame | None = None,
    mt_length: int = 0,
    min_support_per_million: float = 0.0,
    max_points_per_group: int = 0,
) -> tuple[pd.DataFrame, list[str], str]:
    if reads.empty:
        return pd.DataFrame(), [], "Deletion-supporting reads"
    df = apply_cluster_coordinates(reads, clusters, mt_length=mt_length)
    df["left_breakpoint"] = pd.to_numeric(df["left_breakpoint"], errors="coerce")
    df["right_breakpoint"] = pd.to_numeric(df["right_breakpoint"], errors="coerce")
    df["deleted_size"] = pd.to_numeric(df["deleted_size"], errors="coerce")
    df = df.dropna(subset=["left_breakpoint", "right_breakpoint", "deleted_size"])
    sample_cols = ["sample", group_col] if group_col in samples.columns else ["sample"]
    for col in ["normalization_denominator", "normalization_reads", "reads_passed_to_minimap2"]:
        if col in samples.columns:
            sample_cols.append(col)
    df = df.merge(samples[sample_cols], on="sample", how="left")
    df["_plot_group"] = df[group_col].fillna("missing").astype(str) if group_col in df.columns else "all"
    if "normalization_reads" in df.columns or "reads_passed_to_minimap2" in df.columns:
        denom_col = "normalization_reads" if "normalization_reads" in df.columns else "reads_passed_to_minimap2"
        denom = pd.to_numeric(df[denom_col], errors="coerce")
        df["_support_weight"] = np.where(denom > 0, 1_000_000 / denom, 0)
        support_col = "support_per_million_mt_reads"
        support_label = f"Deletion support {per_million_phrase(samples)}"
    else:
        df["_support_weight"] = 1.0
        support_col = "supporting_reads"
        support_label = "Deletion-supporting reads"
    group_columns = ["_plot_group", "left_breakpoint", "right_breakpoint", "deleted_size"]
    if "exact_deletion_id" in df.columns:
        group_columns.append("exact_deletion_id")
    metadata_columns = [
        "affected_feature_label",
        "affected_features",
        "replication_arc_context",
        "major_arc_deleted_bp",
        "minor_arc_deleted_bp",
        "known_deletion_label",
        "deleted_interval",
        "wraps_origin",
    ]
    aggregation = {
        "supporting_reads": ("sample", "size"),
        "support_per_million_mt_reads": ("_support_weight", "sum"),
    }
    aggregation.update({column: (column, "first") for column in metadata_columns if column in df.columns})
    grouped = (
        df.groupby(group_columns, as_index=False, dropna=False)
        .agg(**aggregation)
    )
    grouped["_plot_support"] = pd.to_numeric(grouped[support_col], errors="coerce").fillna(0)
    if support_col == "support_per_million_mt_reads" and min_support_per_million > 0:
        grouped = grouped[grouped["_plot_support"] >= min_support_per_million].copy()
    groups = [group for group in ordered_groups(samples, group_col) if group in set(grouped["_plot_group"])] or sorted(grouped["_plot_group"].unique())
    capped = []
    for group in groups:
        tie_columns = [col for col in ["exact_deletion_id", "left_breakpoint", "right_breakpoint", "deleted_size"] if col in grouped.columns]
        sub = grouped[grouped["_plot_group"] == group].sort_values(
            ["_plot_support", *tie_columns],
            ascending=[False, *([True] * len(tie_columns))],
            kind="mergesort",
        )
        if max_points_per_group > 0:
            sub = sub.head(max_points_per_group)
        capped.append(sub)
    grouped = pd.concat(capped, ignore_index=True) if capped else grouped.iloc[0:0].copy()
    grouped = assign_group_support_ranks(grouped)
    return grouped, groups, support_label


def location_support_scale(display_grouped: pd.DataFrame) -> tuple[float, float, colors.Normalize, list[float]]:
    observed_min, observed_max = rainfall_support_limits(display_grouped["_plot_support"] if "_plot_support" in display_grouped.columns else [])
    support_min, support_max = support_scale_limits(observed_min, observed_max)
    if support_min < support_max:
        support_norm: colors.Normalize = colors.LogNorm(vmin=support_min, vmax=support_max)
    else:
        support_norm = colors.Normalize(vmin=0, vmax=max(support_max, 1.0))
    return support_min, support_max, support_norm, support_legend_values(support_min, support_max)


def draw_location_points(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    support_min: float,
    support_max: float,
    support_norm: colors.Normalize,
    cmap,
    outline_crossing: bool = True,
    artist_prefix: str | None = None,
):
    artist_index = 0

    def mark_artist(artist) -> None:
        nonlocal artist_index
        if artist_prefix:
            artist.set_gid(f"{artist_prefix}-{artist_index}")
            artist_index += 1

    non_origin = data[~data["crosses_origin"]]
    origin = data[data["crosses_origin"]]
    scatter = None
    if not non_origin.empty:
        ordered = non_origin.sort_values("_plot_support", ascending=True, kind="mergesort")
        scatter = ax.scatter(
            ordered[x_col],
            ordered[y_col],
            c=ordered["_plot_support"],
            s=rainfall_point_sizes(ordered["_plot_support"], support_min, support_max),
            cmap=cmap,
            norm=support_norm,
            edgecolors="#17202a",
            linewidths=0.35,
            alpha=1.0,
            zorder=3,
        )
        mark_artist(scatter)
    if not origin.empty:
        for chunk_index, chunk in enumerate(support_ordered_groups(origin)):
            chunk_sizes = rainfall_point_sizes(chunk["_plot_support"], support_min, support_max)
            zbase = 4.0 + chunk_index * 0.02
            origin_scatter = ax.scatter(
                chunk[x_col],
                chunk[y_col],
                c=chunk["_plot_support"],
                s=chunk_sizes,
                cmap=cmap,
                norm=support_norm,
                edgecolors="#17202a",
                linewidths=0.35,
                alpha=1.0,
                zorder=zbase,
            )
            mark_artist(origin_scatter)
            if scatter is None:
                scatter = origin_scatter
            if outline_crossing:
                outer_lw, _ = origin_outline_linewidths(chunk_sizes)
                outline = ax.scatter(chunk[x_col], chunk[y_col], s=chunk_sizes, facecolors="none", edgecolors=ORIGIN_OUTLINE_COLOR, linewidths=outer_lw, alpha=0.98, zorder=zbase + 0.01)
                mark_artist(outline)
    return scatter


def draw_location_rank_labels(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    support_min: float,
    support_max: float,
    support_norm: colors.Normalize,
    cmap,
) -> int:
    if data.empty or "_support_rank" not in data.columns:
        return 0
    work = data.copy()
    work["_marker_area"] = rainfall_point_sizes(work["_plot_support"], support_min, support_max)
    labeled = 0
    occupied = []
    ax.figure.canvas.draw()
    # Highest-support ranks claim space first. Lower ranks are omitted when labels would merge visually.
    for _, row in work.sort_values("_support_rank", ascending=True, kind="mergesort").iterrows():
        rank = int(row["_support_rank"])
        font_size = rank_label_font_size(float(row["_marker_area"]), rank)
        if font_size is None:
            continue
        display_x, display_y = ax.transData.transform((float(row[x_col]), float(row[y_col])))
        marker_radius = np.sqrt(float(row["_marker_area"])) * ax.figure.dpi / 72.0 / 2.0
        text_radius = font_size * ax.figure.dpi / 72.0 * max(1.0, 0.42 * len(str(rank)))
        label_box = (float(display_x), float(display_y), max(text_radius, marker_radius * 0.38))
        if any(rank_label_boxes_overlap(label_box, prior) for prior in occupied):
            continue
        ax.text(
            float(row[x_col]),
            float(row[y_col]),
            str(rank),
            ha="center",
            va="center",
            fontsize=font_size,
            fontweight="bold",
            color=rank_label_color(float(row["_plot_support"]), support_norm, cmap),
            clip_on=True,
            zorder=20,
        )
        ax.texts[-1].set_gid(f"breakpoint-pair-rank-{rank}")
        occupied.append(label_box)
        labeled += 1
    return labeled


def location_rainfall(
    display_grouped: pd.DataFrame,
    groups: list[str],
    features: pd.DataFrame,
    config: dict,
    mt_length: int,
    path: str,
    title_prefix: str,
    x_col: str,
    x_label: str,
    support_label: str,
) -> None:
    clear_location_sidecars(path)
    if display_grouped.empty:
        empty(path, title_prefix, "No exact deletions meet the location plot display threshold")
        return
    plot_features = location_features(features, config)
    genome_length = location_genome_length(mt_length, plot_features, display_grouped)
    work = display_grouped.copy()
    work["crosses_origin"] = work["right_breakpoint"] < work["left_breakpoint"]
    work["circular_midpoint"] = circular_deleted_interval_midpoint(work, genome_length)
    work["_plot_x"] = work[x_col] if x_col in {"left_breakpoint", "right_breakpoint"} else work["circular_midpoint"]
    support_min, support_max, support_norm, legend_values = location_support_scale(work)
    y_axis_min = rainfall_y_axis_min(work["deleted_size"])
    figures: list[plt.Figure] = []
    for group in groups:
        sub = work[work["_plot_group"] == group].copy()
        fig = plt.figure(figsize=(15.2, 7.7), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=[5.1, 1.05], width_ratios=[1, 0.28], hspace=0.04, wspace=0.10)
        ax = fig.add_subplot(grid[0, 0])
        feature_ax = fig.add_subplot(grid[1, 0], sharex=ax)
        legend_ax = fig.add_subplot(grid[:, 1])
        scatter = None
        if sub.empty:
            ax.text(0.5, 0.5, "No exact deletions meet the location plot display threshold", ha="center", va="center", wrap=True, transform=ax.transAxes)
            legend_ax.set_axis_off()
        else:
            cmap = location_support_colormap()
            scatter = draw_location_points(
                ax,
                sub,
                "_plot_x",
                "deleted_size",
                support_min,
                support_max,
                support_norm,
                cmap,
                outline_crossing=True,
                artist_prefix="rainfall-static-points",
            )
        ax.set_title(f"{title_prefix}: {group}")
        ax.set_ylabel("Deleted size (bp)")
        y_max = max(10_000, float(sub["deleted_size"].max()) * 1.18) if not sub.empty else 10_000
        ax.set_ylim(y_axis_min, y_max)
        format_deletion_size_log_axis(ax)
        ax.set_xlim(1, genome_length)
        draw_location_feature_track(feature_ax, plot_features, genome_length, x_min=1, x_max=genome_length)
        feature_ax.set_xlabel(x_label)
        if scatter is not None:
            add_location_legend(
                fig,
                legend_ax,
                scatter,
                legend_values,
                support_min,
                support_max,
                support_label,
                show_origin_outline=bool(sub["crosses_origin"].any()),
                origin_label="cyan outline =\norigin-spanning deletion",
                show_rank_note=False,
            )
        save_location_sidecar(fig, path, group)
        save_rainfall_interactive_sidecar(
            fig,
            ax,
            path,
            group,
            sub,
            "_plot_x",
            support_min,
            support_max,
            support_norm,
            cmap,
            support_label,
        )
        figures.append(fig)
    write_location_plot_pages(figures, path)
    for fig in figures:
        plt.close(fig)


def breakpoint_pair_support_map(display_grouped: pd.DataFrame, groups: list[str], features: pd.DataFrame, config: dict, mt_length: int, path: str, support_label: str) -> None:
    clear_location_sidecars(path)
    if display_grouped.empty:
        empty(path, "Breakpoint-Pair Support Map", "No exact deletions meet the breakpoint-pair display threshold")
        return
    plot_features = location_features(features, config)
    genome_length = location_genome_length(mt_length, plot_features, display_grouped)
    work = display_grouped.copy()
    work["crosses_origin"] = work["right_breakpoint"] < work["left_breakpoint"]
    work["adjusted_right_breakpoint"] = np.where(work["crosses_origin"], work["right_breakpoint"] + genome_length, work["right_breakpoint"])
    support_min, support_max, support_norm, legend_values = location_support_scale(work)
    figures: list[plt.Figure] = []
    for group in groups:
        sub = work[work["_plot_group"] == group].copy()
        fig = plt.figure(figsize=(15.2, 8.0), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=[5.4, 1.05], width_ratios=[1, 0.28], hspace=0.04, wspace=0.10)
        ax = fig.add_subplot(grid[0, 0])
        feature_ax = fig.add_subplot(grid[1, 0], sharex=ax)
        legend_ax = fig.add_subplot(grid[:, 1])
        scatter = None
        if sub.empty:
            ax.text(0.5, 0.5, "No exact deletions meet the breakpoint-pair display threshold", ha="center", va="center", wrap=True, transform=ax.transAxes)
            legend_ax.set_axis_off()
        else:
            grouping = ["left_breakpoint", "right_breakpoint", "adjusted_right_breakpoint", "crosses_origin"]
            aggregations = {
                "_plot_support": ("_plot_support", "sum"),
                "supporting_reads": ("supporting_reads", "sum"),
                "_support_rank": ("_support_rank", "min"),
                "deleted_size": ("deleted_size", "first"),
            }
            if "exact_deletion_id" in sub.columns:
                aggregations["pair_count"] = ("exact_deletion_id", "nunique")
                aggregations["exact_deletion_id"] = (
                    "exact_deletion_id",
                    lambda values: "+".join(sorted(set(values.dropna().astype(str)))),
                )
            else:
                aggregations["pair_count"] = ("left_breakpoint", "size")
            for column in [
                "affected_feature_label",
                "affected_features",
                "replication_arc_context",
                "major_arc_deleted_bp",
                "minor_arc_deleted_bp",
            ]:
                if column in sub.columns:
                    aggregations[column] = (column, "first")
            pairs = sub.groupby(grouping, as_index=False).agg(**aggregations).sort_values(
                "_plot_support", ascending=True, kind="mergesort"
            )
            cmap = location_support_colormap()
            scatter = draw_location_points(ax, pairs, "left_breakpoint", "adjusted_right_breakpoint", support_min, support_max, support_norm, cmap, outline_crossing=True)
        y_max = max(float(genome_length), float(sub["adjusted_right_breakpoint"].max()) if not sub.empty else float(genome_length))
        y_axis_max = y_max + max(250.0, (y_max - 1.0) * 0.035)
        ax.set_xlim(0, genome_length)
        ax.set_ylim(0, y_axis_max)
        if not sub.empty:
            draw_location_rank_labels(ax, pairs, "left_breakpoint", "adjusted_right_breakpoint", support_min, support_max, support_norm, cmap)
        y_ticks, y_labels = adjusted_breakpoint_ticks(genome_length, y_axis_max)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels)
        ax.plot([0, genome_length], [0, genome_length], color="#64748b", linewidth=0.8, linestyle="--", alpha=0.45, zorder=1)
        ax.axhline(genome_length, color="#334155", linewidth=0.9, linestyle=":", alpha=0.75, zorder=2)
        ax.set_title(f"Breakpoint-Pair Support Map: {group}")
        ax.set_ylabel("Deletion end / right breakpoint (bp; labels restart after origin)")
        ax.grid(axis="both", color="#d9dee7", linewidth=0.65, alpha=0.55)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", labelbottom=True)
        draw_location_feature_track(feature_ax, plot_features, genome_length, x_min=1, x_max=genome_length)
        feature_ax.set_xlim(0, genome_length)
        feature_ax.set_xlabel("Deletion start / left breakpoint on mitochondrial genome (bp)")
        if scatter is not None:
            add_location_legend(
                fig,
                legend_ax,
                scatter,
                legend_values,
                support_min,
                support_max,
                support_label,
                show_origin_outline=bool(sub["crosses_origin"].any()),
                origin_label="cyan outline =\norigin-crossing deletion",
                note="Each dot is one unique left/right breakpoint pair.\nPoints above the horizontal line cross the origin.",
            )
        save_location_sidecar(fig, path, group)
        save_breakpoint_pair_interactive_sidecar(
            fig,
            ax,
            path,
            group,
            pairs if not sub.empty else pd.DataFrame(),
            support_min,
            support_max,
            support_label,
        )
        figures.append(fig)
    write_location_plot_pages(figures, path)
    for fig in figures:
        plt.close(fig)


def location_plots(
    reads: pd.DataFrame,
    samples: pd.DataFrame,
    clusters: pd.DataFrame,
    features: pd.DataFrame,
    config: dict,
    mt_length: int,
    group_col: str,
    out_left: str,
    out_right: str,
    out_midpoint: str,
    out_pair_map: str,
    out_endpoint_density: str,
    out_endpoint_density_capped: str,
    min_support_per_million: float = 0.0,
    max_points_per_group: int = 0,
    endpoint_density_bin_size: int = 50,
    endpoint_density_smooth_bins: int = 7,
) -> None:
    display_grouped, groups, support_label = prepare_location_plot_data(
        reads,
        samples,
        group_col,
        clusters=clusters,
        mt_length=mt_length,
        min_support_per_million=min_support_per_million,
        max_points_per_group=max_points_per_group,
    )
    location_rainfall(display_grouped, groups, features, config, mt_length, out_left, "Deletion Rainfall By Directed Left Breakpoint", "left_breakpoint", "Directed left breakpoint on mitochondrial genome (bp)", support_label)
    location_rainfall(display_grouped, groups, features, config, mt_length, out_right, "Deletion Rainfall By Directed Right Breakpoint", "right_breakpoint", "Directed right breakpoint on mitochondrial genome (bp)", support_label)
    location_rainfall(display_grouped, groups, features, config, mt_length, out_midpoint, "Deletion Rainfall By Circular Deleted-Interval Midpoint", "circular_midpoint", "Circular midpoint of deleted interval on mitochondrial genome (bp)", support_label)
    breakpoint_pair_support_map(display_grouped, groups, features, config, mt_length, out_pair_map, support_label)
    endpoint_density_pages(display_grouped, groups, features, config, mt_length, out_endpoint_density, "Pooled Breakpoint Support Density", support_label, bin_size=endpoint_density_bin_size, smooth_bins=endpoint_density_smooth_bins, capped=False)
    endpoint_density_pages(display_grouped, groups, features, config, mt_length, out_endpoint_density_capped, "Pooled Breakpoint Support Density, Capped Scale", support_label, bin_size=endpoint_density_bin_size, smooth_bins=endpoint_density_smooth_bins, capped=True)


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
    save_ordination_interactive(fig, ax, path, df, xlabel, ylabel, group_col)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--config", default="")
    parser.add_argument("--mt-length", type=int, default=0)
    parser.add_argument("--all-reads", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--burden", required=True)
    parser.add_argument("--exact-mtpm", required=True)
    parser.add_argument("--affected-raw", required=True)
    parser.add_argument("--affected-mtpm", required=True)
    parser.add_argument("--impact-class-mtpm", required=True)
    parser.add_argument("--gene-pair-mtpm", default="")
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
    parser.add_argument("--out-rainfall-left", required=True)
    parser.add_argument("--out-rainfall-right", required=True)
    parser.add_argument("--out-rainfall-midpoint", required=True)
    parser.add_argument("--out-breakpoint-pair-map", required=True)
    parser.add_argument("--out-endpoint-density", required=True)
    parser.add_argument("--out-endpoint-density-capped", required=True)
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
    parser.add_argument("--out-gene-pair-pca", default="")
    parser.add_argument("--rainfall-min-support-per-million", type=float, default=0.0)
    parser.add_argument("--rainfall-max-points-per-group", type=int, default=0)
    parser.add_argument("--endpoint-density-bin-size", type=int, default=50)
    parser.add_argument("--endpoint-density-smooth-bins", type=int, default=7)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t")
    config = read_yaml_safe(args.config)
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
    gene_pair_mtpm = read_tsv_safe(args.gene_pair_mtpm) if args.gene_pair_mtpm else pd.DataFrame()
    per_gene = read_tsv_safe(args.per_gene_burden)
    comparison = normalize_deletion_ids(read_tsv_safe(args.exact_comparison))
    support_label = f"Support {per_million_phrase(burden)}"
    burden_label = f"Deletion-supporting reads {per_million_phrase(burden)}"

    burden_plot(burden, args.group_column, args.out_burden, "deletion_support_per_million_mt_reads", "Total Deletion Burden", burden_label)
    burden_plot(burden, args.group_column, args.out_unique_count, "unique_exact_deletions", "Distinct Exact Deletions", "Distinct exact deletion calls")
    factorial_interaction_plot(burden, args.out_burden_factorial, "deletion_support_per_million_mt_reads", "Deletion Burden: Age By Treatment", burden_label)
    factorial_interaction_plot(burden, args.out_unique_factorial, "unique_exact_deletions", "Distinct Exact Deletions: Age By Treatment", "Distinct exact deletion calls")
    cluster_coordinate_reads = apply_cluster_coordinates(reads, clusters, mt_length=args.mt_length)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_unweighted, "Deletion Size Distribution, Unweighted", weighted=False)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_weighted, "Deletion Size Distribution, Support-Weighted", weighted=True)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_weighted_log, "Deletion Size Distribution, Support-Weighted Log Scale", weighted=True, log_y=True)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_small, "Small Deletion Size Distribution (<1 kb)", weighted=True, size_max=999)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_medium, "Medium Deletion Size Distribution (1-5 kb)", weighted=True, size_min=1000, size_max=4999)
    size_distribution(cluster_coordinate_reads, burden, args.group_column, args.out_size_large, "Large Deletion Size Distribution (>=5 kb)", weighted=True, size_min=5000)
    location_plots(
        reads,
        burden,
        clusters,
        features,
        config,
        args.mt_length,
        args.group_column,
        args.out_rainfall_left,
        args.out_rainfall_right,
        args.out_rainfall_midpoint,
        args.out_breakpoint_pair_map,
        args.out_endpoint_density,
        args.out_endpoint_density_capped,
        min_support_per_million=args.rainfall_min_support_per_million,
        max_points_per_group=args.rainfall_max_points_per_group,
        endpoint_density_bin_size=args.endpoint_density_bin_size,
        endpoint_density_smooth_bins=args.endpoint_density_smooth_bins,
    )
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
    if args.out_gene_pair_pca:
        if gene_pair_pca_enabled(config):
            ordination(gene_pair_mtpm, samples, args.group_column, args.out_gene_pair_pca, "Mitochondrial Gene-Pair PCA", "pca")
        else:
            empty(
                args.out_gene_pair_pca,
                "Mitochondrial Gene-Pair PCA",
                "This view applies only when the short-read RNA STAR evidence stream is enabled",
            )


if __name__ == "__main__":
    main()
