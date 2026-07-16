#!/usr/bin/env python3
"""Generate circular mitochondrial deletion chord plots for workflow reports."""

from __future__ import annotations

import argparse
import copy
import math
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors, ticker
from matplotlib.cm import ScalarMappable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, PathPatch, Patch, Wedge
from matplotlib.path import Path as MplPath

from circular_deletions import replication_arc_annotation
from plot_deletion_results import (
    MITOCHONDRIAL_FEATURE_COLORS,
    deduplicate_evidence_reads,
    location_features,
    location_support_colormap,
    location_support_scale,
    normalize_deletion_ids,
    prepare_location_plot_data,
    read_tsv_safe,
    read_yaml_safe,
    support_tick_label,
)


FEATURE_LABELS = {
    "protein_coding": "protein-coding gene",
    "rRNA": "rRNA gene",
    "tRNA": "tRNA gene",
    "region": "D-loop/control region",
}
COMPARISON_NEGATIVE_COLOR = "#2b6cb0"
COMPARISON_NEUTRAL_COLOR = "#e5e7eb"
COMPARISON_POSITIVE_COLOR = "#e4572e"


def safe_token(value: object) -> str:
    return "".join(character if character.isalnum() or character in "_.-" else "_" for character in str(value))


def integer_or_default(value: object, default: int = 0) -> int:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return default if pd.isna(numeric) else int(numeric)


def coordinate_angle(position: float, genome_length: int) -> float:
    """Map 1-based mtDNA coordinates clockwise with coordinate 1 at 12 o'clock."""
    return math.pi / 2 - 2 * math.pi * (float(position) - 1) / float(genome_length)


def circle_point(position: float, radius: float, genome_length: int) -> tuple[float, float]:
    angle = coordinate_angle(position, genome_length)
    return radius * math.cos(angle), radius * math.sin(angle)


def clockwise_wedge_angles(start: float, end: float, genome_length: int) -> tuple[float, float]:
    """Return Matplotlib's counter-clockwise wedge angles for a forward feature."""
    start_angle = math.degrees(coordinate_angle(start, genome_length))
    end_angle = math.degrees(coordinate_angle(end, genome_length))
    return end_angle, start_angle


def chord_path(left: float, right: float, radius: float, genome_length: int) -> MplPath:
    start = np.asarray(circle_point(left, radius, genome_length), dtype=float)
    end = np.asarray(circle_point(right, radius, genome_length), dtype=float)
    control_fraction = 0.10
    vertices = [start, start * control_fraction, end * control_fraction, end]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    return MplPath(vertices, codes)


def circular_feature_rows(features: pd.DataFrame) -> pd.DataFrame:
    """Keep biological features plus the configured D-loop/control segments."""
    if features.empty:
        return features.copy()
    keep = features[features["class"].isin({"protein_coding", "rRNA", "tRNA"})].copy()
    dloop = features[(features["class"] == "region") & (features["name"] == "D-loop/control")].copy()
    return pd.concat([keep, dloop], ignore_index=True).drop_duplicates().sort_values(["start", "end", "name"])


def add_feature_ring(ax: plt.Axes, features: pd.DataFrame, genome_length: int) -> None:
    outer_radius = 1.0
    ring_width = 0.085
    for feature_index, (_, row) in enumerate(features.iterrows()):
        start = float(row["start"])
        end = float(row["end"])
        theta1, theta2 = clockwise_wedge_angles(start, end, genome_length)
        feature_class = str(row["class"])
        patch = Wedge(
            (0, 0),
            outer_radius,
            theta1,
            theta2,
            width=ring_width,
            facecolor=MITOCHONDRIAL_FEATURE_COLORS.get(
                feature_class, MITOCHONDRIAL_FEATURE_COLORS["other"]
            ),
            edgecolor="#ffffff",
            linewidth=0.45,
            zorder=6,
        )
        patch.set_gid(feature_dom_id(feature_index, row))
        ax.add_patch(patch)
    ax.add_patch(Circle((0, 0), outer_radius, fill=False, edgecolor="#334155", linewidth=0.8, zorder=7))
    ax.add_patch(Circle((0, 0), outer_radius - ring_width, fill=False, edgecolor="#334155", linewidth=0.55, zorder=7))


def add_coordinate_ticks(ax: plt.Axes, genome_length: int) -> None:
    ticks = [1, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000]
    for coordinate in ticks:
        if coordinate > genome_length:
            continue
        inner = circle_point(coordinate, 1.005, genome_length)
        outer = circle_point(coordinate, 1.026, genome_length)
        label = circle_point(coordinate, 1.058, genome_length)
        ax.plot([inner[0], outer[0]], [inner[1], outer[1]], color="#475569", linewidth=0.65, zorder=8)
        text = "1" if coordinate == 1 else f"{coordinate // 1000}k"
        ax.text(label[0], label[1], text, ha="center", va="center", fontsize=7.0, color="#475569", zorder=8)


def add_chords(
    ax: plt.Axes,
    calls: pd.DataFrame,
    genome_length: int,
    support_norm: colors.Normalize,
    cmap: colors.Colormap,
    subset_label: str,
    group: str,
) -> None:
    for _, row in calls.sort_values("_plot_support", ascending=True, kind="mergesort").iterrows():
        color = cmap(support_norm(float(row["_plot_support"])))
        patch = PathPatch(
            chord_path(row["left_breakpoint"], row["right_breakpoint"], 0.908, genome_length),
            fill=False,
            edgecolor=color,
            linewidth=0.82,
            alpha=0.70,
            capstyle="round",
            zorder=2,
        )
        patch.set_gid(chord_dom_id(subset_label, group, row))
        ax.add_patch(patch)


def comparison_colormap() -> colors.Colormap:
    return colors.LinearSegmentedColormap.from_list(
        "comparison_difference",
        [COMPARISON_NEGATIVE_COLOR, COMPARISON_NEUTRAL_COLOR, COMPARISON_POSITIVE_COLOR],
    )


def comparison_chord_dom_id(subset_label: str, left_group: str, right_group: str, row: pd.Series) -> str:
    return "comparison_chord__{}__{}_vs_{}__rank_{}__{}".format(
        safe_token(subset_label),
        safe_token(left_group),
        safe_token(right_group),
        int(row["_comparison_rank"]),
        safe_token(row["exact_deletion_id"]),
    )


def add_comparison_chords(
    ax: plt.Axes,
    calls: pd.DataFrame,
    genome_length: int,
    comparison_norm: colors.Normalize,
    subset_label: str,
    left_group: str,
    right_group: str,
) -> None:
    cmap = comparison_colormap()
    for _, row in calls.sort_values("_absolute_difference", ascending=True, kind="mergesort").iterrows():
        patch = PathPatch(
            chord_path(row["left_breakpoint"], row["right_breakpoint"], 0.908, genome_length),
            fill=False,
            edgecolor=cmap(comparison_norm(float(row["difference_per_million_mt_reads"]))),
            linewidth=0.88,
            alpha=0.76,
            capstyle="round",
            zorder=2,
        )
        patch.set_gid(comparison_chord_dom_id(subset_label, left_group, right_group, row))
        ax.add_patch(patch)


def chord_dom_id(subset_label: str, group: str, row: pd.Series) -> str:
    return "chord__{}__{}__rank_{}__{}".format(
        safe_token(subset_label),
        safe_token(group),
        int(row["_support_rank"]),
        safe_token(row["exact_deletion_id"]),
    )


def feature_dom_id(feature_index: int, row: pd.Series) -> str:
    return f"feature__{feature_index}__{safe_token(row['name'])}__{int(row['start'])}_{int(row['end'])}"


def add_svg_chord_metadata(
    svg_path: Path,
    calls: pd.DataFrame,
    subset_label: str,
    group: str,
    baseline_ids: set[str] | None = None,
    support_label: str = "Deletion support",
) -> None:
    """Attach filter values and hover text to the Matplotlib chord groups."""
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    root.set("data-group", group)
    root.set("data-support-label", support_label)
    nodes_by_id = {node.get("id"): node for node in root.iter() if node.get("id")}
    for _, row in calls.iterrows():
        node_id = chord_dom_id(subset_label, group, row)
        node = nodes_by_id.get(node_id)
        if node is None:
            raise ValueError(f"missing SVG chord group {node_id}")
        support = float(row["_plot_support"])
        observations = int(row["supporting_reads"])
        node.set("class", "deletion-chord")
        node.set("data-support", f"{support:.12g}")
        node.set("data-support-label", support_label)
        node.set("data-observations", str(observations))
        node.set("data-rank", str(int(row["_support_rank"])))
        node.set("data-deletion-id", str(row["exact_deletion_id"]))
        node.set("data-left-breakpoint", str(int(row["left_breakpoint"])))
        node.set("data-right-breakpoint", str(int(row["right_breakpoint"])))
        node.set("data-deleted-size", str(int(row["deleted_size"])))
        affected_features = row.get("affected_feature_label", row.get("affected_features", ""))
        node.set("data-affected-features", "" if pd.isna(affected_features) else str(affected_features))
        node.set("data-arc-context", str(row.get("replication_arc_context", "not annotated")))
        node.set("data-major-arc-bp", str(integer_or_default(row.get("major_arc_deleted_bp", 0))))
        node.set("data-minor-arc-bp", str(integer_or_default(row.get("minor_arc_deleted_bp", 0))))
        node.set(
            "data-baseline",
            "1" if baseline_ids is None or str(row["exact_deletion_id"]) in baseline_ids else "0",
        )
        accessible_label = (
            f"rank {int(row['_support_rank'])}: {row['exact_deletion_id']}; "
            f"{support_label.lower()} {support:.4g}; {observations} supporting observations"
        )
        node.set("aria-label", accessible_label)
        description = ET.Element("{http://www.w3.org/2000/svg}desc")
        description.text = accessible_label
        visible_path = next(
            (child for child in node if child.tag == "{http://www.w3.org/2000/svg}path"),
            None,
        )
        node.insert(0, description)
        if visible_path is not None:
            hit_path = copy.deepcopy(visible_path)
            hit_path.set("class", "chord-hit-target")
            hit_path.set("style", "fill: none; stroke: transparent; stroke-width: 9; pointer-events: stroke")
            node.insert(1, hit_path)
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)


def add_svg_comparison_metadata(
    svg_path: Path,
    calls: pd.DataFrame,
    subset_label: str,
    left_group: str,
    right_group: str,
) -> None:
    """Attach group-comparison values and accessible hover targets to chords."""
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    root.set("data-left-group", left_group)
    root.set("data-right-group", right_group)
    nodes_by_id = {node.get("id"): node for node in root.iter() if node.get("id")}
    numeric_fields = {
        "left-mean": "left_mean_per_million_mt_reads",
        "right-mean": "right_mean_per_million_mt_reads",
        "difference": "difference_per_million_mt_reads",
        "absolute-difference": "_absolute_difference",
        "replicate-p": "p_value",
        "replicate-q": "q_value_bh",
        "depth-p": "read_depth_fisher_p",
        "depth-q": "read_depth_fisher_q_value_bh",
    }
    integer_fields = {
        "left-observations": "left_total_supporting_reads",
        "right-observations": "right_total_supporting_reads",
        "total-observations": "_total_supporting_observations",
        "samples-with-signal": "samples_with_signal",
    }
    for _, row in calls.iterrows():
        node_id = comparison_chord_dom_id(subset_label, left_group, right_group, row)
        node = nodes_by_id.get(node_id)
        if node is None:
            raise ValueError(f"missing SVG comparison chord group {node_id}")
        node.set("class", "comparison-chord")
        node.set("data-rank", str(int(row["_comparison_rank"])))
        node.set("data-deletion-id", str(row["exact_deletion_id"]))
        node.set("data-left-group", left_group)
        node.set("data-right-group", right_group)
        node.set("data-left-breakpoint", str(int(row["left_breakpoint"])))
        node.set("data-right-breakpoint", str(int(row["right_breakpoint"])))
        node.set("data-deleted-size", str(int(row["deleted_size"])))
        affected_features = row.get("affected_feature_label", row.get("affected_features", ""))
        node.set("data-affected-features", "" if pd.isna(affected_features) else str(affected_features))
        node.set("data-arc-context", str(row.get("replication_arc_context", "not annotated")))
        node.set("data-major-arc-bp", str(integer_or_default(row.get("major_arc_deleted_bp", 0))))
        node.set("data-minor-arc-bp", str(integer_or_default(row.get("minor_arc_deleted_bp", 0))))
        node.set("data-abundance-change", str(row.get("abundance_change", "not annotated")))
        known_deletion = row.get("known_deletion_label", "")
        node.set("data-known-deletion", "" if pd.isna(known_deletion) else str(known_deletion))
        for data_name, column in numeric_fields.items():
            value = pd.to_numeric(pd.Series([row.get(column)]), errors="coerce").iloc[0]
            node.set(f"data-{data_name}", "" if pd.isna(value) else f"{float(value):.12g}")
        for data_name, column in integer_fields.items():
            value = pd.to_numeric(pd.Series([row.get(column)]), errors="coerce").iloc[0]
            node.set(f"data-{data_name}", "" if pd.isna(value) else str(int(value)))
        accessible_label = (
            f"comparison rank {int(row['_comparison_rank'])}: {row['exact_deletion_id']}; "
            f"{right_group} minus {left_group} normalized mean difference "
            f"{float(row['difference_per_million_mt_reads']):.4g}"
        )
        node.set("aria-label", accessible_label)
        description = ET.Element("{http://www.w3.org/2000/svg}desc")
        description.text = accessible_label
        visible_path = next(
            (child for child in node if child.tag == "{http://www.w3.org/2000/svg}path"),
            None,
        )
        node.insert(0, description)
        if visible_path is not None:
            hit_path = copy.deepcopy(visible_path)
            hit_path.set("class", "chord-hit-target")
            hit_path.set("style", "fill: none; stroke: transparent; stroke-width: 9; pointer-events: stroke")
            node.insert(1, hit_path)
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)


def add_svg_feature_metadata(svg_path: Path, features: pd.DataFrame) -> None:
    """Make mitochondrial annotation segments addressable by the HTML tooltip."""
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    nodes_by_id = {node.get("id"): node for node in root.iter() if node.get("id")}
    for feature_index, (_, row) in enumerate(features.iterrows()):
        node_id = feature_dom_id(feature_index, row)
        node = nodes_by_id.get(node_id)
        if node is None:
            raise ValueError(f"missing SVG feature group {node_id}")
        feature_class = str(row["class"])
        feature_label = FEATURE_LABELS.get(feature_class, feature_class.replace("_", " "))
        accessible_label = f"{row['name']}; {feature_label}; coordinates {int(row['start'])} to {int(row['end'])}"
        node.set("class", "mt-feature")
        node.set("data-feature-name", str(row["name"]))
        node.set("data-feature-type", feature_label)
        node.set("data-feature-start", str(int(row["start"])))
        node.set("data-feature-end", str(int(row["end"])))
        node.set("aria-label", accessible_label)
        description = ET.Element("{http://www.w3.org/2000/svg}desc")
        description.text = accessible_label
        node.insert(0, description)
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)


def feature_legend_handles(features: pd.DataFrame) -> list[Patch]:
    present = set(features["class"].astype(str))
    return [
        Patch(
            facecolor=MITOCHONDRIAL_FEATURE_COLORS[feature_class],
            edgecolor="#ffffff",
            label=FEATURE_LABELS[feature_class],
        )
        for feature_class in ["protein_coding", "rRNA", "tRNA", "region"]
        if feature_class in present
    ]


def compact_colorbar_ticks(support_norm: colors.Normalize) -> list[float]:
    """Use scale endpoints and intervening powers of ten for a compact legend."""
    support_min = float(support_norm.vmin)
    support_max = float(support_norm.vmax)
    if support_min <= 0 or support_max <= 0 or support_min >= support_max:
        return [support_max]
    start_exponent = int(math.ceil(math.log10(support_min)))
    end_exponent = int(math.floor(math.log10(support_max)))
    values = [support_min]
    values.extend(float(10**exponent) for exponent in range(start_exponent, end_exponent + 1))
    values.append(support_max)
    compact = []
    for value in values:
        if support_min <= value <= support_max and not any(np.isclose(value, prior) for prior in compact):
            compact.append(value)
    return compact


def comparison_color_scale(calls: pd.DataFrame) -> colors.SymLogNorm:
    values = pd.to_numeric(calls.get("difference_per_million_mt_reads", pd.Series(dtype=float)), errors="coerce")
    observed = float(values.abs().max()) if values.notna().any() else 1.0
    observed = max(observed, 1e-6)
    exponent = math.floor(math.log10(observed))
    fraction = observed / (10**exponent)
    nice_fraction = 1 if fraction <= 1 else 2 if fraction <= 2 else 5 if fraction <= 5 else 10
    limit = nice_fraction * (10**exponent)
    return colors.SymLogNorm(
        linthresh=max(limit / 100, 1e-8),
        linscale=0.75,
        vmin=-limit,
        vmax=limit,
        base=10,
    )


def signed_tick_label(value: float) -> str:
    if np.isclose(value, 0):
        return "0"
    sign = "-" if value < 0 else "+"
    return f"{sign}{support_tick_label(abs(value))}"


def plot_comparison_chord_page(
    calls: pd.DataFrame,
    features: pd.DataFrame,
    genome_length: int,
    comparison_norm: colors.SymLogNorm,
    left_group: str,
    right_group: str,
    subset_label: str,
    output_stem: Path,
    multipage_pdf: PdfPages | None = None,
) -> None:
    fig = plt.figure(figsize=(10.6, 8.7), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.36], wspace=0.09)
    ax = fig.add_subplot(grid[0, 0])
    legend_ax = fig.add_subplot(grid[0, 1])
    cmap = comparison_colormap()

    add_comparison_chords(
        ax,
        calls,
        genome_length,
        comparison_norm,
        subset_label,
        left_group,
        right_group,
    )
    add_feature_ring(ax, features, genome_length)
    add_coordinate_ticks(ax, genome_length)
    ax.set_aspect("equal")
    ax.set_xlim(-1.14, 1.14)
    ax.set_ylim(-1.14, 1.14)
    ax.set_axis_off()
    ax.set_title(
        f"Exact Deletion Group Comparison Chords: {subset_label}\n{right_group} compared with {left_group}",
        fontsize=13,
        pad=10,
    )
    ax.text(
        0,
        -1.125,
        f"{len(calls):,} exact deletion comparisons loaded; HTML controls filter chords",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#475569",
    )

    legend_ax.set_axis_off()
    legend_ax.legend(
        handles=feature_legend_handles(features),
        title="Mitochondrial annotation",
        loc="upper left",
        bbox_to_anchor=(0.0, 0.96),
        borderaxespad=0,
        frameon=False,
        fontsize=8.5,
        title_fontsize=9.0,
    )
    legend_center = 0.46
    legend_ax.text(
        legend_center,
        0.72,
        "Chord color",
        ha="center",
        fontsize=9.0,
        fontweight="bold",
        transform=legend_ax.transAxes,
    )
    legend_ax.text(
        legend_center,
        0.655,
        f"Mean normalized support difference\n{right_group} minus {left_group}\n(symmetric log scale)",
        ha="center",
        va="center",
        fontsize=8.2,
        linespacing=1.2,
        transform=legend_ax.transAxes,
    )
    color_ax = legend_ax.inset_axes([0.03, 0.565, 0.86, 0.042])
    colorbar = fig.colorbar(ScalarMappable(norm=comparison_norm, cmap=cmap), cax=color_ax, orientation="horizontal")
    limit = float(comparison_norm.vmax)
    colorbar_ticks = [-limit, -limit / 10, 0.0, limit / 10, limit]
    colorbar.set_ticks(colorbar_ticks)
    colorbar.ax.set_xticklabels([signed_tick_label(value) for value in colorbar_ticks], fontsize=7.2)
    colorbar.ax.xaxis.set_minor_locator(ticker.NullLocator())
    colorbar.ax.tick_params(axis="x", which="major", length=3, width=0.7, pad=2)
    legend_ax.text(
        0.03,
        0.50,
        f"Blue: higher in {left_group}\nOrange: higher in {right_group}",
        fontsize=8.2,
        linespacing=1.25,
        transform=legend_ax.transAxes,
    )
    display_note = "All exact deletion comparison rows."
    how_to_read = (
        "Chord\n"
        "Joins the directed left and right breakpoints of\n"
        "one exact deletion comparison.\n\n"
        "Display set\n"
        f"{display_note}\n\n"
        "Coordinates\n"
        "Coordinate 1 is at 12 o'clock; coordinates\n"
        "increase clockwise."
    )
    legend_ax.text(0.03, 0.425, "How to read", fontsize=9.0, fontweight="bold", transform=legend_ax.transAxes)
    legend_ax.text(
        0.03,
        0.39,
        how_to_read,
        fontsize=8.2,
        linespacing=1.3,
        va="top",
        transform=legend_ax.transAxes,
    )

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    if multipage_pdf is not None:
        multipage_pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    add_svg_comparison_metadata(
        output_stem.with_suffix(".svg"),
        calls,
        subset_label,
        left_group,
        right_group,
    )
    add_svg_feature_metadata(output_stem.with_suffix(".svg"), features)


def plot_chord_page(
    calls: pd.DataFrame,
    features: pd.DataFrame,
    genome_length: int,
    support_norm: colors.Normalize,
    support_label: str,
    group: str,
    subset_label: str,
    output_stem: Path,
    multipage_pdf: PdfPages | None,
    *,
    interactive: bool = False,
    baseline_ids: set[str] | None = None,
) -> None:
    fig = plt.figure(figsize=(10.6, 8.7), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.36], wspace=0.09)
    ax = fig.add_subplot(grid[0, 0])
    legend_ax = fig.add_subplot(grid[0, 1])
    cmap = location_support_colormap()

    add_chords(ax, calls, genome_length, support_norm, cmap, subset_label, group)
    add_feature_ring(ax, features, genome_length)
    add_coordinate_ticks(ax, genome_length)
    ax.set_aspect("equal")
    ax.set_xlim(-1.14, 1.14)
    ax.set_ylim(-1.14, 1.14)
    ax.set_axis_off()
    ax.set_title(f"Circular Breakpoint Chords: {subset_label}\n{group}", fontsize=13, pad=10)
    ax.text(
        0,
        -1.125,
        (
            f"{len(calls):,} exact deletions loaded; HTML controls filter chords"
            if interactive
            else f"{len(calls):,} exact deletions displayed"
        ),
        ha="center",
        va="center",
        fontsize=8.5,
        color="#475569",
    )

    legend_ax.set_axis_off()
    legend_ax.legend(
        handles=feature_legend_handles(features),
        title="Mitochondrial annotation",
        loc="upper left",
        bbox_to_anchor=(0.0, 0.96),
        borderaxespad=0,
        frameon=False,
        fontsize=8.5,
        title_fontsize=9.0,
    )
    legend_center = 0.46
    legend_ax.text(legend_center, 0.72, "Chord color", ha="center", fontsize=9.0, fontweight="bold", transform=legend_ax.transAxes)
    legend_ax.text(
        legend_center,
        0.665,
        f"{support_label}\n(log scale)",
        ha="center",
        va="center",
        fontsize=8.2,
        linespacing=1.2,
        transform=legend_ax.transAxes,
    )
    color_ax = legend_ax.inset_axes([0.03, 0.585, 0.86, 0.042])
    colorbar = fig.colorbar(ScalarMappable(norm=support_norm, cmap=cmap), cax=color_ax, orientation="horizontal")
    colorbar_ticks = compact_colorbar_ticks(support_norm)
    colorbar.set_ticks(colorbar_ticks)
    colorbar.ax.set_xticklabels([support_tick_label(value) for value in colorbar_ticks], fontsize=7.2)
    colorbar.ax.xaxis.set_minor_locator(ticker.NullLocator())
    colorbar.ax.tick_params(axis="x", which="major", length=3, width=0.7, pad=2)
    if interactive:
        display_note = "All calls meeting the rainfall support threshold;\nthe HTML view does not apply the rainfall count cap."
    else:
        display_note = "Same group-specific support threshold and\nmaximum count as the rainfall plots."
    how_to_read = (
        "Chord\n"
        "Joins the directed left and right breakpoints of\n"
        "one exact deletion.\n\n"
        "Display set\n"
        f"{display_note}\n\n"
        "Coordinates\n"
        "Coordinate 1 is at 12 o'clock; coordinates\n"
        "increase clockwise."
    )
    legend_ax.text(0.03, 0.49, "How to read", fontsize=9.0, fontweight="bold", transform=legend_ax.transAxes)
    legend_ax.text(
        0.03,
        0.455,
        how_to_read,
        fontsize=8.2,
        linespacing=1.3,
        va="top",
        transform=legend_ax.transAxes,
    )

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    if multipage_pdf is not None:
        multipage_pdf.savefig(fig, bbox_inches="tight")
    if not interactive:
        fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    add_svg_chord_metadata(
        output_stem.with_suffix(".svg"),
        calls,
        subset_label,
        group,
        baseline_ids=baseline_ids,
        support_label=support_label,
    )
    add_svg_feature_metadata(output_stem.with_suffix(".svg"), features)


def annotate_replication_arcs(calls: pd.DataFrame, config: dict, genome_length: int) -> pd.DataFrame:
    out = calls.copy()
    arc_rows = [
        replication_arc_annotation(config, int(row.left_breakpoint), int(row.right_breakpoint), genome_length)
        for row in out.itertuples()
    ]
    arc_columns = ["replication_arc_context", "minor_arc_deleted_bp", "major_arc_deleted_bp"]
    out = out.drop(columns=[column for column in arc_columns if column in out.columns])
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(arc_rows)], axis=1)


def prepare_comparison_calls(comparison: pd.DataFrame, config: dict, genome_length: int) -> pd.DataFrame:
    required_columns = [
        "exact_deletion_id",
        "left_group",
        "right_group",
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "difference_per_million_mt_reads",
        "left_total_supporting_reads",
        "right_total_supporting_reads",
    ]
    if comparison.empty and not len(comparison.columns):
        comparison = pd.DataFrame(columns=required_columns)
    missing = sorted(set(required_columns) - set(comparison.columns))
    if missing:
        raise ValueError(f"comparison table is missing required columns: {', '.join(missing)}")
    out = comparison.copy()
    numeric_columns = [
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "left_mean_per_million_mt_reads",
        "right_mean_per_million_mt_reads",
        "difference_per_million_mt_reads",
        "left_total_supporting_reads",
        "right_total_supporting_reads",
        "samples_with_signal",
        "p_value",
        "q_value_bh",
        "read_depth_fisher_p",
        "read_depth_fisher_q_value_bh",
    ]
    for column in numeric_columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    out = out.dropna(
        subset=[
            "left_breakpoint",
            "right_breakpoint",
            "deleted_size",
            "difference_per_million_mt_reads",
        ]
    ).copy()
    out["_absolute_difference"] = out["difference_per_million_mt_reads"].abs()
    out["_total_supporting_observations"] = (
        out["left_total_supporting_reads"].fillna(0) + out["right_total_supporting_reads"].fillna(0)
    ).astype(int)
    out = out[out["_total_supporting_observations"] >= 1].copy()
    ranked = []
    for _, pair in out.groupby(["left_group", "right_group"], sort=False, dropna=False):
        pair = pair.sort_values(
            ["_absolute_difference", "_total_supporting_observations", "exact_deletion_id"],
            ascending=[False, False, True],
            kind="mergesort",
        ).copy()
        pair["_comparison_rank"] = np.arange(1, len(pair) + 1)
        ranked.append(pair)
    out = pd.concat(ranked, ignore_index=True) if ranked else out.iloc[0:0].copy()
    return annotate_replication_arcs(out, config, genome_length)


def validate_coordinates(genome_length: int) -> None:
    top = circle_point(1, 1.0, genome_length)
    quarter = circle_point(1 + genome_length / 4, 1.0, genome_length)
    assert np.allclose(top, (0.0, 1.0), atol=1e-9)
    assert np.allclose(quarter, (1.0, 0.0), atol=1e-9)


def validate_calls(calls: pd.DataFrame, genome_length: int) -> None:
    assert calls["exact_deletion_id"].notna().all()
    assert calls["_plot_support"].gt(0).all()
    for row in calls.itertuples():
        path = chord_path(row.left_breakpoint, row.right_breakpoint, 0.908, genome_length)
        expected_left = circle_point(row.left_breakpoint, 0.908, genome_length)
        expected_right = circle_point(row.right_breakpoint, 0.908, genome_length)
        assert np.allclose(path.vertices[0], expected_left, atol=1e-9)
        assert np.allclose(path.vertices[-1], expected_right, atol=1e-9)


def validate_comparison_calls(calls: pd.DataFrame, genome_length: int) -> None:
    assert calls["exact_deletion_id"].notna().all()
    assert calls["_total_supporting_observations"].ge(1).all()
    assert calls["_absolute_difference"].ge(0).all()
    for row in calls.itertuples():
        path = chord_path(row.left_breakpoint, row.right_breakpoint, 0.908, genome_length)
        assert np.allclose(path.vertices[0], circle_point(row.left_breakpoint, 0.908, genome_length), atol=1e-9)
        assert np.allclose(path.vertices[-1], circle_point(row.right_breakpoint, 0.908, genome_length), atol=1e-9)


def audit_table(calls: pd.DataFrame, genome_length: int) -> pd.DataFrame:
    out = calls.copy()
    out["left_angle_degrees_clockwise_from_coordinate_1"] = (
        (pd.to_numeric(out["left_breakpoint"]) - 1) / genome_length * 360
    )
    out["right_angle_degrees_clockwise_from_coordinate_1"] = (
        (pd.to_numeric(out["right_breakpoint"]) - 1) / genome_length * 360
    )
    columns = [
        "_plot_group",
        "_support_rank",
        "exact_deletion_id",
        "left_breakpoint",
        "right_breakpoint",
        "deleted_size",
        "supporting_reads",
        "_plot_support",
        "replication_arc_context",
        "minor_arc_deleted_bp",
        "major_arc_deleted_bp",
        "left_angle_degrees_clockwise_from_coordinate_1",
        "right_angle_degrees_clockwise_from_coordinate_1",
    ]
    return out[[column for column in columns if column in out.columns]].sort_values(
        ["_plot_group", "_support_rank"], kind="mergesort"
    )


def empty_plot_page(title: str, message: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(10.6, 8.7), constrained_layout=True)
    ax.set_axis_off()
    ax.set_title(title, fontsize=13, pad=10)
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color="#52606d")
    return fig


def clear_sidecars(aggregate_path: Path) -> None:
    for sidecar in aggregate_path.parent.glob(f"{aggregate_path.stem}__*"):
        if sidecar.is_file() or sidecar.is_symlink():
            sidecar.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--burden", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--group-column", default="")
    parser.add_argument("--genome-length", required=True, type=int)
    parser.add_argument("--rainfall-min-support-per-million", type=float, default=0.0)
    parser.add_argument("--rainfall-max-points-per-group", type=int, default=0)
    parser.add_argument("--out-location", required=True, type=Path)
    parser.add_argument("--out-comparison", required=True, type=Path)
    parser.add_argument("--out-displayed-table", required=True, type=Path)
    parser.add_argument("--out-interactive-table", required=True, type=Path)
    parser.add_argument("--out-comparison-table", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    args = parser.parse_args()

    config = read_yaml_safe(str(args.config))
    samples = pd.read_csv(args.burden, sep="\t")
    raw_features = pd.read_csv(args.features, sep="\t")
    observations = deduplicate_evidence_reads(normalize_deletion_ids(pd.read_csv(args.observations, sep="\t")))
    clusters = normalize_deletion_ids(pd.read_csv(args.clusters, sep="\t"))
    comparison_calls = prepare_comparison_calls(
        normalize_deletion_ids(read_tsv_safe(str(args.comparison))), config, args.genome_length
    )

    displayed, groups, support_label = prepare_location_plot_data(
        observations,
        samples,
        args.group_column,
        clusters=clusters,
        mt_length=args.genome_length,
        min_support_per_million=args.rainfall_min_support_per_million,
        max_points_per_group=args.rainfall_max_points_per_group,
    )
    interactive_calls, interactive_groups, interactive_support_label = prepare_location_plot_data(
        observations,
        samples,
        args.group_column,
        clusters=clusters,
        mt_length=args.genome_length,
        min_support_per_million=args.rainfall_min_support_per_million,
        max_points_per_group=0,
    )
    if interactive_groups != groups or interactive_support_label != support_label:
        raise ValueError("static and interactive chord datasets resolved inconsistent groups or support units")

    displayed = annotate_replication_arcs(displayed, config, args.genome_length)
    interactive_calls = annotate_replication_arcs(interactive_calls, config, args.genome_length)
    features = circular_feature_rows(location_features(raw_features, config))
    _, _, support_norm, _ = location_support_scale(displayed)

    validate_coordinates(args.genome_length)
    validate_calls(displayed, args.genome_length)
    validate_calls(interactive_calls, args.genome_length)
    validate_comparison_calls(comparison_calls, args.genome_length)

    for path in [
        args.out_location,
        args.out_comparison,
        args.out_displayed_table,
        args.out_interactive_table,
        args.out_comparison_table,
        args.out_summary,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)
    clear_sidecars(args.out_location)
    clear_sidecars(args.out_comparison)
    audit_table(displayed, args.genome_length).to_csv(args.out_displayed_table, sep="\t", index=False)
    audit_table(interactive_calls, args.genome_length).to_csv(args.out_interactive_table, sep="\t", index=False)
    comparison_calls.to_csv(args.out_comparison_table, sep="\t", index=False)

    summary_rows: list[dict[str, object]] = []
    with PdfPages(args.out_location) as location_pdf:
        if not groups:
            figure = empty_plot_page(
                "Circular Breakpoint Chords", "No exact deletions meet the rainfall display threshold."
            )
            location_pdf.savefig(figure, bbox_inches="tight")
            figure.savefig(args.out_location.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
            plt.close(figure)
        for group in groups:
            safe_group = safe_token(group)
            baseline_group = displayed[displayed["_plot_group"] == group].copy()
            interactive_group = interactive_calls[interactive_calls["_plot_group"] == group].copy()
            baseline_ids = set(baseline_group["exact_deletion_id"].astype(str))
            static_stem = args.out_location.parent / f"{args.out_location.stem}__{safe_group}"
            interactive_stem = args.out_location.parent / f"{args.out_location.stem}__{safe_group}__interactive"
            plot_chord_page(
                baseline_group,
                features,
                args.genome_length,
                support_norm,
                support_label,
                group,
                "All rainfall-displayed deletions",
                static_stem,
                location_pdf,
            )
            plot_chord_page(
                interactive_group,
                features,
                args.genome_length,
                support_norm,
                support_label,
                group,
                "All rainfall-displayed deletions",
                interactive_stem,
                None,
                interactive=True,
                baseline_ids=baseline_ids,
            )
            summary_rows.append(
                {
                    "group": group,
                    "rainfall_displayed_deletions": len(baseline_group),
                    "interactive_loaded_deletions": len(interactive_group),
                    "rainfall_min_support_per_million": args.rainfall_min_support_per_million,
                    "rainfall_max_points_per_group": args.rainfall_max_points_per_group,
                    "support_min": float(interactive_group["_plot_support"].min())
                    if not interactive_group.empty
                    else np.nan,
                    "support_max": float(interactive_group["_plot_support"].max())
                    if not interactive_group.empty
                    else np.nan,
                }
            )
    if groups:
        first_static = args.out_location.parent / f"{args.out_location.stem}__{safe_token(groups[0])}.svg"
        shutil.copyfile(first_static, args.out_location.with_suffix(".svg"))

    comparison_norm = comparison_color_scale(comparison_calls)
    comparison_pairs = list(
        comparison_calls[["left_group", "right_group"]].drop_duplicates().itertuples(index=False, name=None)
    )
    with PdfPages(args.out_comparison) as comparison_pdf:
        if not comparison_pairs:
            figure = empty_plot_page(
                "Exact Deletion Group Comparison Chords", "No exact deletion group comparisons are available."
            )
            comparison_pdf.savefig(figure, bbox_inches="tight")
            figure.savefig(args.out_comparison.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
            plt.close(figure)
        for left_group, right_group in comparison_pairs:
            pair_calls = comparison_calls[
                (comparison_calls["left_group"] == left_group)
                & (comparison_calls["right_group"] == right_group)
            ].copy()
            pair_token = f"{safe_token(left_group)}_vs_{safe_token(right_group)}"
            output_stem = args.out_comparison.parent / f"{args.out_comparison.stem}__{pair_token}"
            plot_comparison_chord_page(
                pair_calls,
                features,
                args.genome_length,
                comparison_norm,
                str(left_group),
                str(right_group),
                "All deletion intervals",
                output_stem,
                comparison_pdf,
            )
    if comparison_pairs:
        first_left, first_right = comparison_pairs[0]
        first_svg = args.out_comparison.parent / (
            f"{args.out_comparison.stem}__{safe_token(first_left)}_vs_{safe_token(first_right)}.svg"
        )
        shutil.copyfile(first_svg, args.out_comparison.with_suffix(".svg"))

    pd.DataFrame(
        summary_rows,
        columns=[
            "group",
            "rainfall_displayed_deletions",
            "interactive_loaded_deletions",
            "rainfall_min_support_per_million",
            "rainfall_max_points_per_group",
            "support_min",
            "support_max",
        ],
    ).to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
