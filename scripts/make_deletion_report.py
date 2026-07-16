#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
from pathlib import Path

import pandas as pd
import yaml

from common import deep_update, ensure_parent
from circular_deletions import (
    circular_closed_interval_pieces,
    circular_distance,
    configured_deletion_targets,
    configured_replication_arcs,
    interval_length,
    pos_within_circular_window,
)


READ_LIST_COLUMNS = [
    "sample",
    "read_id",
    "exact_deletion_id",
    "junction_id",
    "breakpoint_pair_id",
    "left_breakpoint",
    "right_breakpoint",
    "deleted_size",
    "wraps_origin",
    "complement_deleted_size",
    "arc_assignment_method",
    "direction_status",
    "rotation_agreement",
    "read_left_breakpoint",
    "read_right_breakpoint",
    "read_deleted_size",
    "directed_left_breakpoint",
    "directed_right_breakpoint",
    "directed_deleted_size",
    "reported_left_breakpoint",
    "reported_right_breakpoint",
    "deleted_interval",
    "rotation_name",
    "rotation_start",
    "strand",
    "source",
    "left_anchor_length",
    "right_anchor_length",
    "min_anchor_length",
    "left_mapq",
    "right_mapq",
    "min_mapq",
    "query_overlap_bp",
    "query_gap_bp",
    "query_first_start",
    "query_first_end",
    "query_second_start",
    "query_second_end",
    "left_cigar",
    "right_cigar",
    "left_is_primary",
    "right_is_primary",
    "left_is_secondary",
    "right_is_secondary",
    "left_is_supplementary",
    "right_is_supplementary",
    "left_alignment_score",
    "right_alignment_score",
    "left_edit_distance",
    "right_edit_distance",
]

CONFIGURED_SEARCH_READ_LIST_COLUMNS = [
    "sample",
    "read_id",
    "deletion_id",
    "deletion_name",
    "search_strategy",
    "mate",
    "matched_sequence_ids",
    "matched_orientation",
]

REPORT_COLUMN_LABELS = {
    "known_deletion_label": "configured_deletion_target_label",
    "known_deletion_match_reason": "configured_deletion_target_match_reason",
    "known_deletion": "configured_deletion_target",
    "matched_known_deletion_label": "matched_configured_target_label",
    "left_mean_per_million_mt_reads": "left_mean_per_million_normalized_reads",
    "right_mean_per_million_mt_reads": "right_mean_per_million_normalized_reads",
    "difference_per_million_mt_reads": "difference_per_million_normalized_reads",
    "support_per_million_mt_reads": "support_per_million_normalized_reads",
    "total_support_per_million_mt_reads": "total_support_per_million_normalized_reads",
    "deletion_support_per_million_mt_reads": "deletion_support_per_million_normalized_reads",
    "mean_deletion_support_per_million_mt_reads": "mean_deletion_support_per_million_normalized_reads",
    "large_deletion_support_per_million_mt_reads": "large_deletion_support_per_million_normalized_reads",
    "small_lt_1kb_support_per_million_mt_reads": "small_lt_1kb_support_per_million_normalized_reads",
    "medium_1kb_to_4999bp_support_per_million_mt_reads": "medium_1kb_to_4999bp_support_per_million_normalized_reads",
    "large_ge_5kb_support_per_million_mt_reads": "large_ge_5kb_support_per_million_normalized_reads",
    "read_depth_fisher_p": "denominator_depth_fisher_p",
    "read_depth_fisher_q_value_bh": "denominator_depth_fisher_q_value_bh",
}


def normalization_mode(burden: pd.DataFrame, config: dict | None = None) -> str:
    if not burden.empty and "normalization_denominator" in burden.columns:
        values = burden["normalization_denominator"].dropna().astype(str).unique().tolist()
        if values:
            return values[0]
    return str((config or {}).get("analysis", {}).get("normalization_denominator", "total_usable_reads"))


def normalization_phrase(burden: pd.DataFrame, config: dict | None = None) -> str:
    if normalization_mode(burden, config) == "mt_evidence_reads":
        return "per million mitochondrial-evidence reads"
    return "per million usable reads"


def normalization_definition(burden: pd.DataFrame, config: dict | None = None) -> str:
    if normalization_mode(burden, config) == "mt_evidence_reads":
        return "The main normalized plots and burden tables divide deletion-supporting reads by the number of first-pass retained mitochondrial-evidence reads in each sample, then scale to one million reads."
    return "The main normalized plots and burden tables divide deletion-supporting reads by the total usable reads after read preparation in each sample, then scale to one million reads."


def rainfall_display_definition(config: dict, burden: pd.DataFrame) -> str:
    plots = config.get("plots", {}) or {}
    min_support = float(plots.get("rainfall_min_support_per_million", 0) or 0)
    max_points = int(plots.get("rainfall_max_points_per_group", 0) or 0)
    if min_support > 0:
        rule = (
            "Each full-size figure shows exact deletions whose normalized support is at least "
            f"{min_support:g} {normalization_phrase(burden, config)} in that plotted group."
        )
    else:
        rule = "Each full-size figure shows all exact deletions available for that plotted group."
    if max_points > 0:
        rule += f" If more than {max_points:,} deletions pass that threshold, only the {max_points:,} highest-support deletions are drawn."
    return rule


def read_table(path: str) -> pd.DataFrame:
    if not path or not Path(path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def load_report_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        supplied = yaml.safe_load(handle) or {}
    defaults_path = Path("config/defaults.yaml")
    if defaults_path.exists():
        with defaults_path.open("r", encoding="utf-8") as handle:
            defaults = yaml.safe_load(handle) or {}
        return deep_update(defaults, supplied)
    return supplied


def safe_sidecar_name(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip()).strip("._")
    return text or "unnamed"


def read_key(value: object) -> str:
    text = str(value).strip().split()[0] if str(value).strip() else ""
    return text.removesuffix("/1").removesuffix("/2")


def write_exact_deletion_read_lists(junction_reads: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.tsv"):
        old.unlink()
    manifest_rows = []
    if junction_reads.empty or "exact_deletion_id" not in junction_reads.columns or "read_id" not in junction_reads.columns:
        manifest = pd.DataFrame(columns=["exact_deletion_id", "read_count", "read_list_file"])
        manifest.to_csv(out_dir / "manifest.tsv", sep="\t", index=False)
        return manifest
    available_columns = [col for col in READ_LIST_COLUMNS if col in junction_reads.columns]
    for exact_deletion_id, group in junction_reads.groupby("exact_deletion_id", dropna=False, sort=False):
        exact_deletion_id = str(exact_deletion_id)
        filename = f"{safe_sidecar_name(exact_deletion_id)}.read_names.tsv"
        read_list = group[available_columns].copy()
        if {"sample", "read_id"}.issubset(read_list.columns):
            read_list = read_list.drop_duplicates(["sample", "read_id"], keep="first")
        else:
            read_list = read_list.drop_duplicates(keep="first")
        read_list = read_list.sort_values([col for col in ["sample", "read_id"] if col in read_list.columns], kind="mergesort")
        read_list.to_csv(out_dir / filename, sep="\t", index=False)
        manifest_rows.append(
            {
                "exact_deletion_id": exact_deletion_id,
                "read_count": len(read_list),
                "read_list_file": filename,
            }
        )
    manifest = pd.DataFrame(manifest_rows).sort_values("exact_deletion_id", kind="mergesort")
    manifest.to_csv(out_dir / "manifest.tsv", sep="\t", index=False)
    return manifest


def write_configured_sequence_read_lists(summary: pd.DataFrame, hits: pd.DataFrame, out_dir: Path) -> dict[tuple[int, str], str]:
    if summary.empty or hits.empty or "matching_reads" not in summary.columns:
        return {}
    if not {"sample", "deletion_id"}.issubset(summary.columns) or not {"sample", "deletion_id", "read_id"}.issubset(hits.columns):
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)
    valid_hits = hits.copy()
    valid_hits["read_id"] = valid_hits["read_id"].fillna("").astype(str).str.strip()
    valid_hits = valid_hits[
        (valid_hits["read_id"] != "")
        & (valid_hits["read_id"] != "not_recorded_literal_search_count")
    ].copy()
    if valid_hits.empty:
        return {}
    available_columns = [col for col in CONFIGURED_SEARCH_READ_LIST_COLUMNS if col in valid_hits.columns]
    html_cells: dict[tuple[int, str], str] = {}
    for idx, row in summary.iterrows():
        sample = str(row.get("sample", ""))
        deletion_id = str(row.get("deletion_id", ""))
        group = valid_hits[
            (valid_hits["sample"].fillna("").astype(str) == sample)
            & (valid_hits["deletion_id"].fillna("").astype(str) == deletion_id)
        ].copy()
        if group.empty:
            continue
        read_list = group[available_columns].drop_duplicates(keep="first")
        sort_cols = [col for col in ["sample", "read_id", "mate"] if col in read_list.columns]
        if sort_cols:
            read_list = read_list.sort_values(sort_cols, kind="mergesort")
        filename = f"configured_sequence__{safe_sidecar_name(sample)}__{safe_sidecar_name(deletion_id)}.read_names.tsv"
        read_list.to_csv(out_dir / filename, sep="\t", index=False)
        count = html.escape(str(row.get("matching_reads", len(read_list))))
        href = html.escape(f"read_lists/{filename}")
        title = html.escape(f"Open read names matching configured search {deletion_id} in {sample}")
        html_cells[(idx, "matching_reads")] = f'<a class="read-list-link" href="{href}" title="{title}">{count}</a>'
    return html_cells


def configured_search_targets(config: dict) -> dict[str, dict[str, object]]:
    targets: dict[str, dict[str, object]] = {}
    known = config.get("analysis", {}).get("known_deletions", []) or []
    analysis = config.get("analysis", {}) or {}
    species = config.get("dataset", {}).get("species", "")
    mt_length = int(config.get("references", {}).get(species, {}).get("mt_length", 0) or 0)
    known_tolerance = []
    for item in known:
        try:
            left = int(item.get("left_breakpoint"))
            right = int(item.get("right_breakpoint"))
            size = int(item.get("deleted_size", abs(right - left) - 1))
        except (TypeError, ValueError):
            continue
        known_tolerance.append(
            {
                "name": item.get("name", ""),
                "left": left,
                "right": right,
                "size": size,
                "breakpoint_tolerance_bp": int(item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", 100))),
                "size_tolerance_bp": int(item.get("size_tolerance_bp", 150)),
            }
        )
    for item in config.get("analysis", {}).get("known_sequence_searches", []) or []:
        deletion_id = str(item.get("id", ""))
        try:
            if item.get("left_breakpoint") not in {None, ""} and item.get("right_breakpoint") not in {None, ""}:
                left = int(item["left_breakpoint"])
                right = int(item["right_breakpoint"])
            else:
                text = " ".join([deletion_id, str(item.get("name", "")), str(item.get("description", ""))])
                match = re.search(r"(\d{2,6})\D+(\d{2,6})", text)
                if not match:
                    continue
                left = int(match.group(1))
                right = int(match.group(2))
            size = int(item.get("deleted_size", circular_distance(left, right, mt_length) if mt_length else abs(right - left) - 1))
        except (TypeError, ValueError):
            continue
        breakpoint_tolerance = int(item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", analysis.get("sequence_search_breakpoint_tolerance_bp", 100))))
        size_tolerance = int(item.get("size_tolerance_bp", analysis.get("sequence_search_size_tolerance_bp", 150)))
        matched_target = ""
        for known_item in known_tolerance:
            if mt_length:
                left_ok = pos_within_circular_window(left, known_item["left"], known_item["breakpoint_tolerance_bp"], mt_length)
                right_ok = pos_within_circular_window(right, known_item["right"], known_item["breakpoint_tolerance_bp"], mt_length)
            else:
                left_ok = abs(left - known_item["left"]) <= known_item["breakpoint_tolerance_bp"]
                right_ok = abs(right - known_item["right"]) <= known_item["breakpoint_tolerance_bp"]
            if left_ok and right_ok:
                breakpoint_tolerance = int(known_item["breakpoint_tolerance_bp"])
                size_tolerance = int(known_item["size_tolerance_bp"])
                matched_target = str(known_item["name"])
                break
        targets[deletion_id] = {
            "deletion_id": deletion_id,
            "deletion_name": item.get("name", ""),
            "target_left_breakpoint": left,
            "target_right_breakpoint": right,
            "target_deleted_size": size,
            "breakpoint_tolerance_bp": breakpoint_tolerance,
            "size_tolerance_bp": size_tolerance,
            "matched_configured_target_label": matched_target,
        }
    return targets


def deletion_matches_target(row: pd.Series, target: dict[str, object]) -> bool:
    try:
        left = int(float(row.get("left_breakpoint")))
        right = int(float(row.get("right_breakpoint")))
        size = int(float(row.get("deleted_size")))
        target_left = int(target["target_left_breakpoint"])
        target_right = int(target["target_right_breakpoint"])
        target_size = int(target["target_deleted_size"])
        breakpoint_tolerance = int(target["breakpoint_tolerance_bp"])
        size_tolerance = int(target["size_tolerance_bp"])
    except (TypeError, ValueError, KeyError):
        return False
    direct = abs(left - target_left) <= breakpoint_tolerance and abs(right - target_right) <= breakpoint_tolerance
    reciprocal = abs(left - target_right) <= breakpoint_tolerance and abs(right - target_left) <= breakpoint_tolerance
    return (direct or reciprocal) and abs(size - target_size) <= size_tolerance


def write_read_set(path: Path, rows: list[dict[str, object]]) -> None:
    ensure_parent(path)
    columns = ["read_key", "read_id", "sample", "source", "deletion_id", "exact_deletion_id", "junction_id"]
    pd.DataFrame(rows, columns=columns).drop_duplicates().sort_values(
        [col for col in ["sample", "read_key", "source"] if col in columns],
        kind="mergesort",
    ).to_csv(path, sep="\t", index=False)


def sequence_remap_overlap_table(config: dict, known_hits: pd.DataFrame, junction_reads: pd.DataFrame, out_dir: Path) -> tuple[pd.DataFrame, dict[tuple[int, str], str]]:
    targets = configured_search_targets(config)
    if not targets or known_hits.empty:
        return pd.DataFrame(), {}
    out_dir.mkdir(parents=True, exist_ok=True)
    sequence_hits = known_hits.copy()
    if not {"sample", "deletion_id", "read_id"}.issubset(sequence_hits.columns):
        return pd.DataFrame(), {}
    if junction_reads.empty:
        junction_reads = pd.DataFrame(columns=["sample", "read_id", "left_breakpoint", "right_breakpoint", "deleted_size", "exact_deletion_id", "junction_id"])
    elif not {"sample", "read_id", "left_breakpoint", "right_breakpoint", "deleted_size"}.issubset(junction_reads.columns):
        return pd.DataFrame(), {}
    sequence_hits = sequence_hits[sequence_hits["read_id"].fillna("").astype(str).ne("not_recorded_literal_search_count")].copy()
    if sequence_hits.empty:
        return pd.DataFrame(), {}
    sequence_hits["_read_key"] = sequence_hits["read_id"].map(read_key)
    junction = junction_reads.copy()
    junction["_read_key"] = junction["read_id"].map(read_key) if "read_id" in junction.columns else ""
    rows = []
    html_cells: dict[tuple[int, str], str] = {}
    for (sample, deletion_id), seq_group in sequence_hits.groupby(["sample", "deletion_id"], dropna=False, sort=False):
        deletion_id = str(deletion_id)
        target = targets.get(deletion_id)
        if not target:
            continue
        sample = str(sample)
        remap_group = junction[junction["sample"].astype(str).eq(sample)].copy() if "sample" in junction.columns else junction.copy()
        if not remap_group.empty:
            remap_group = remap_group[remap_group.apply(lambda row: deletion_matches_target(row, target), axis=1)].copy()
        seq_keys = {key for key in seq_group["_read_key"].dropna().astype(str) if key}
        remap_keys = {key for key in remap_group["_read_key"].dropna().astype(str) if key}
        shared = seq_keys & remap_keys
        seq_only = seq_keys - remap_keys
        remap_only = remap_keys - seq_keys
        base = f"configured_vs_remap__{safe_sidecar_name(sample)}__{safe_sidecar_name(deletion_id)}"

        def rows_for(keys: set[str], source_name: str) -> list[dict[str, object]]:
            read_rows = []
            if source_name == "sequence_search":
                source_df = seq_group[seq_group["_read_key"].isin(keys)]
            elif source_name == "remap_near_target":
                source_df = remap_group[remap_group["_read_key"].isin(keys)]
            else:
                source_df = pd.concat(
                    [
                        seq_group[seq_group["_read_key"].isin(keys)],
                        remap_group[remap_group["_read_key"].isin(keys)],
                    ],
                    ignore_index=True,
                    sort=False,
                )
            for _, read_row in source_df.iterrows():
                read_rows.append(
                    {
                        "read_key": read_row.get("_read_key", ""),
                        "read_id": read_row.get("read_id", ""),
                        "sample": read_row.get("sample", sample),
                        "source": source_name,
                        "deletion_id": deletion_id,
                        "exact_deletion_id": read_row.get("exact_deletion_id", ""),
                        "junction_id": read_row.get("junction_id", ""),
                    }
                )
            return read_rows

        files = {
            "sequence_search_reads": (seq_keys, "sequence_search"),
            "remap_nearby_reads": (remap_keys, "remap_near_target"),
            "shared_reads": (shared, "both"),
            "sequence_only_reads": (seq_only, "sequence_search"),
            "remap_only_reads": (remap_only, "remap_near_target"),
        }
        file_names = {}
        for label, (keys, source_name) in files.items():
            filename = f"{base}__{label}.tsv"
            write_read_set(out_dir / filename, rows_for(keys, source_name))
            file_names[label] = filename
        exact_ids = sorted(set(remap_group.get("exact_deletion_id", pd.Series(dtype=str)).dropna().astype(str))) if not remap_group.empty else []
        rows.append(
            {
                "sample": sample,
                "configured_sequence_search": deletion_id,
                "deletion_name": target.get("deletion_name", ""),
                "target_left_breakpoint": target.get("target_left_breakpoint", ""),
                "target_right_breakpoint": target.get("target_right_breakpoint", ""),
                "target_deleted_size": target.get("target_deleted_size", ""),
                "breakpoint_tolerance_bp": target.get("breakpoint_tolerance_bp", ""),
                "size_tolerance_bp": target.get("size_tolerance_bp", ""),
                "sequence_search_reads": len(seq_keys),
                "remap_nearby_reads": len(remap_keys),
                "shared_reads": len(shared),
                "sequence_only_reads": len(seq_only),
                "remap_only_reads": len(remap_only),
                "overlap_fraction_of_sequence_search": len(shared) / len(seq_keys) if seq_keys else 0,
                "overlap_fraction_of_remap_nearby": len(shared) / len(remap_keys) if remap_keys else 0,
                "matched_remap_exact_deletions": ";".join(exact_ids[:20]),
                "_files": file_names,
            }
        )
    if not rows:
        return pd.DataFrame(), {}
    table = pd.DataFrame(rows)
    for idx, row in table.iterrows():
        for col in ["sequence_search_reads", "remap_nearby_reads", "shared_reads", "sequence_only_reads", "remap_only_reads"]:
            filename = row["_files"].get(col, "")
            if filename:
                html_cells[(idx, col)] = f'<a class="read-list-link" href="read_lists/{html.escape(filename)}">{html.escape(str(row[col]))}</a>'
    table = table.drop(columns=["_files"])
    return table, html_cells


def table_html(df: pd.DataFrame, rows: int | None = 200, html_cells: dict[tuple[int, str], str] | None = None) -> str:
    if df.empty:
        return '<p class="empty">No rows.</p>'
    view = df.copy() if rows is None else df.head(rows).copy()
    view = view.drop(columns=[col for col in ["fastq_1", "fastq_2"] if col in view.columns], errors="ignore")
    view = view.dropna(axis=1, how="all")
    view = view.where(pd.notna(view), "")
    view = view.replace("nan", "")
    replacements = {}
    if html_cells:
        for token_index, ((row_index, col), cell_html) in enumerate(html_cells.items()):
            if row_index in view.index and col in view.columns:
                token = f"__HTML_CELL_{token_index}__"
                view[col] = view[col].astype(object)
                view.loc[row_index, col] = token
                replacements[token] = cell_html
    table = view.to_html(index=False, classes="data-table", escape=True, border=0)
    for original, display in REPORT_COLUMN_LABELS.items():
        table = table.replace(f"<th>{html.escape(original)}</th>", f"<th>{html.escape(display)}</th>")
    table = table.replace("<table ", f'<table data-row-count="{len(view)}" ', 1)
    for token, cell_html in replacements.items():
        table = table.replace(token, cell_html)
    return f'<div class="table-wrap">{table}</div>'


def add_sort_columns(work: pd.DataFrame) -> pd.DataFrame:
    out = work.copy()
    if "known_deletion_label" in out.columns:
        out["_known_deletion_rank"] = out["known_deletion_label"].fillna("").astype(str).eq("")
    if "abundance_change" in out.columns:
        out["_presence_change_rank"] = ~out["abundance_change"].fillna("").astype(str).str.startswith("present_only")
    for col in [
        "read_depth_fisher_q_value_bh",
        "read_depth_fisher_p",
        "q_value_bh",
        "p_value",
        "fisher_presence_p",
        "mann_whitney_p",
    ]:
        if col in out.columns:
            out[f"_{col}_sort"] = pd.to_numeric(out[col], errors="coerce").fillna(float("inf"))
    for col in [
        "difference_per_million_mt_reads",
        "log2_fold_change_right_over_left",
        "total_supporting_reads",
        "right_total_supporting_reads",
        "left_total_supporting_reads",
        "support_per_million_mt_reads",
        "deletion_support_per_million_mt_reads",
        "unique_exact_deletions",
        "samples_with_signal",
        "deleted_size",
    ]:
        if col in out.columns:
            values = pd.to_numeric(out[col], errors="coerce")
            if col in {"difference_per_million_mt_reads", "log2_fold_change_right_over_left"}:
                values = values.abs()
            out[f"_{col}_sort"] = values.fillna(float("-inf"))
    return out


def default_sort_table(df: pd.DataFrame, title: str) -> pd.DataFrame:
    if df.empty:
        return df
    work = add_sort_columns(df)
    sort_cols: list[str] = []
    ascending: list[bool] = []
    candidates = [
        ("_known_deletion_rank", True),
        ("_presence_change_rank", True),
        ("_read_depth_fisher_q_value_bh_sort", True),
        ("_read_depth_fisher_p_sort", True),
        ("_q_value_bh_sort", True),
        ("_p_value_sort", True),
        ("_fisher_presence_p_sort", True),
        ("_difference_per_million_mt_reads_sort", False),
        ("_log2_fold_change_right_over_left_sort", False),
        ("_total_supporting_reads_sort", False),
        ("_right_total_supporting_reads_sort", False),
        ("_left_total_supporting_reads_sort", False),
        ("_support_per_million_mt_reads_sort", False),
        ("_deletion_support_per_million_mt_reads_sort", False),
        ("_unique_exact_deletions_sort", False),
        ("_samples_with_signal_sort", False),
        ("_deleted_size_sort", False),
    ]
    for col, asc in candidates:
        if col in work.columns:
            sort_cols.append(col)
            ascending.append(asc)
    if not sort_cols:
        return df
    sorted_work = work.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    return sorted_work[df.columns].reset_index(drop=True)


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fmt_int(value) -> str:
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError):
        return ""


def fmt_float(value, digits: int = 3) -> str:
    try:
        return f"{float(value):,.{digits}g}"
    except (TypeError, ValueError):
        return ""


def compact_sample_table(samples: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    preferred = [
        "sample",
        "sample_name",
        "biological_replicate",
        "dataset",
        "species",
        "layout",
        "condition",
        "deletion_status",
        "age",
        "treatment",
        "tissue",
        "run_accession",
    ]
    columns = []
    for col in ["sample", *group_columns, *preferred]:
        if col in samples.columns and col not in columns:
            columns.append(col)
    return samples[columns].copy() if columns else samples.copy()


def group_count_table(samples: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if not group_col or group_col not in samples.columns:
        return pd.DataFrame()
    table = (
        samples[group_col]
        .fillna("missing")
        .astype(str)
        .value_counts()
        .rename_axis(group_col)
        .reset_index(name="sample_count")
    )
    if {"age", "treatment"}.issubset(samples.columns):
        reps = samples[[group_col, "age", "treatment"]].drop_duplicates(group_col)
        table = table.merge(reps, on=group_col, how="left")
        table["_age_num"] = pd.to_numeric(table["age"].astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce").fillna(1e9)
        table["_treatment_rank"] = table["treatment"].astype(str).str.lower().map(lambda value: 0 if "control" in value or value in {"ctrl", "vehicle", "untreated"} else 1)
        table = table.sort_values(["_age_num", "_treatment_rank", group_col]).drop(columns=["_age_num", "_treatment_rank"])
    else:
        table = table.sort_values(group_col)
    return table


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (json.JSONDecodeError, OSError):
        return {}


def read_preparation_table(samples: pd.DataFrame, results_dir: Path) -> pd.DataFrame:
    if samples.empty or "sample" not in samples.columns:
        return pd.DataFrame()
    rows = []
    for _, sample_row in samples.iterrows():
        sample = str(sample_row.get("sample", ""))
        qc_dir = results_dir / "qc" / sample
        fastp = read_json(qc_dir / "fastp.json")
        decision = read_json(qc_dir / "qc_decision.json")
        read_input = read_json(qc_dir / "read_input.json")
        summary = fastp.get("summary", {}) if isinstance(fastp, dict) else {}
        before = summary.get("before_filtering", {}) if isinstance(summary, dict) else {}
        after = summary.get("after_filtering", {}) if isinstance(summary, dict) else {}
        r1_before = fastp.get("read1_before_filtering", {}) if isinstance(fastp, dict) else {}
        r2_before = fastp.get("read2_before_filtering", {}) if isinstance(fastp, dict) else {}
        filtering = fastp.get("filtering_result", {}) if isinstance(fastp, dict) else {}
        adapter = fastp.get("adapter_cutting", {}) if isinstance(fastp, dict) else {}
        trimmed = decision.get("trimmed", "") if decision else ""
        if trimmed == "":
            trimmed = "not recorded"
        rows.append(
            {
                "sample": sample,
                "layout": sample_row.get("layout", read_input.get("declared_layout", "")),
                "fastp_trimming_run": trimmed,
                "minimum_length_after_trimming": decision.get("minimum_length", ""),
                "sequencing_cycles": summary.get("sequencing", ""),
                "read1_cycles_raw": r1_before.get("total_cycles", ""),
                "read2_cycles_raw": r2_before.get("total_cycles", ""),
                "read1_mean_length_raw": before.get("read1_mean_length", ""),
                "read2_mean_length_raw": before.get("read2_mean_length", ""),
                "read1_mean_length_after_fastp": after.get("read1_mean_length", ""),
                "read2_mean_length_after_fastp": after.get("read2_mean_length", ""),
                "total_reads_raw": before.get("total_reads", ""),
                "total_reads_after_fastp": after.get("total_reads", ""),
                "q30_rate_raw": before.get("q30_rate", ""),
                "q30_rate_after_fastp": after.get("q30_rate", ""),
                "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads", ""),
                "too_short_reads_removed": filtering.get("too_short_reads", ""),
            }
        )
    return pd.DataFrame(rows)


def feature_class(row: pd.Series) -> str:
    name = str(row.get("gene_name", row.get("feature", ""))).lower()
    feature_type = str(row.get("feature_type", "")).lower()
    if feature_type == "region":
        return "configured region"
    if "trna" in feature_type or name.startswith(("mt-t", "trn")):
        return "tRNA"
    if "rrna" in feature_type or name.startswith(("mt-r", "rrn")):
        return "rRNA"
    if "protein_coding" in feature_type or feature_type in {"cds", "gene"} and name.startswith(("mt-co", "mt-cy", "mt-nd", "mt-atp")):
        return "protein_coding"
    if feature_type in {"cds", "gene"} and name:
        return "gene"
    return feature_type or "feature"


def configured_region_table(config: dict) -> pd.DataFrame:
    rows = []
    for item in config.get("analysis", {}).get("mt_regions", []) or []:
        rows.append(
            {
                "region_name": item.get("name", ""),
                "start": item.get("start", ""),
                "end": item.get("end", ""),
                "description": item.get("reason", ""),
            }
        )
    return pd.DataFrame(rows)


def configured_replication_arc_table(config: dict) -> pd.DataFrame:
    species = str((config.get("dataset", {}) or {}).get("species", ""))
    reference = (config.get("references", {}) or {}).get(species, {}) or {}
    mt_length = int(reference.get("mt_length", 0) or 0)
    rows = []
    for item in configured_replication_arcs(config):
        try:
            start = int(item["start"])
            end = int(item["end"])
            length = interval_length(circular_closed_interval_pieces(start, end, mt_length))
        except (KeyError, TypeError, ValueError):
            start = item.get("start", "")
            end = item.get("end", "")
            length = ""
        rows.append(
            {
                "arc_name": item.get("display_name", item.get("name", "")),
                "start": start,
                "end": end,
                "wraps_coordinate_origin": "yes" if isinstance(start, int) and isinstance(end, int) and start > end else "no",
                "length_bp": length,
                "boundary_definition": item.get("boundary_definition", ""),
            }
        )
    return pd.DataFrame(rows)


def apply_display_aliases_to_features(features: pd.DataFrame, config: dict) -> pd.DataFrame:
    if features.empty:
        return features.copy()
    out = features.copy()
    out["display_name"] = out.apply(lambda row: str(row.get("gene_name") or row.get("gene_id") or row.get("product") or "").strip(), axis=1)
    for item in config.get("annotations", {}).get("feature_aliases", []) or []:
        display = str(item.get("display_name", "")).strip()
        if not display:
            continue
        mask = pd.Series(True, index=out.index)
        if item.get("raw_name"):
            mask &= out.get("gene_name", "").astype(str).eq(str(item["raw_name"]))
        if item.get("gene_id"):
            mask &= out.get("gene_id", "").astype(str).eq(str(item["gene_id"]))
        if item.get("start") is not None:
            mask &= pd.to_numeric(out.get("start"), errors="coerce").eq(int(item["start"]))
        if item.get("end") is not None:
            mask &= pd.to_numeric(out.get("end"), errors="coerce").eq(int(item["end"]))
        out.loc[mask, "display_name"] = display
    return out


def configured_deletion_target_table(config: dict) -> pd.DataFrame:
    species = config.get("dataset", {}).get("species", "")
    mt_length = config.get("references", {}).get(species, {}).get("mt_length", 0)
    try:
        targets = configured_deletion_targets(config, int(mt_length))
    except (TypeError, ValueError):
        targets = []
    rows = []
    for item in targets:
        rows.append(
            {
                "name": item.get("name", ""),
                "source": item.get("source", ""),
                "left_breakpoint": item.get("left_breakpoint", ""),
                "right_breakpoint": item.get("right_breakpoint", ""),
                "deleted_size": item.get("deleted_size", ""),
                "breakpoint_tolerance_bp": item.get("breakpoint_tolerance_bp", item.get("tolerance_bp", "")),
                "size_tolerance_bp": item.get("size_tolerance_bp", ""),
            }
        )
    return pd.DataFrame(rows)


def known_sequence_search_table(config: dict) -> pd.DataFrame:
    rows = []
    for item in config.get("analysis", {}).get("known_sequence_searches", []) or []:
        strategy = item.get("search_strategy", {}) or {}
        sequences = item.get("search_sequences", []) or []
        rows.append(
            {
                "deletion_id": item.get("id", ""),
                "name": item.get("name", ""),
                "strategy": strategy.get("type", ""),
                "search_sequence_count": len(sequences),
                "description": item.get("description", ""),
            }
        )
    return pd.DataFrame(rows)


def simplified_feature_table(features: pd.DataFrame, config: dict) -> pd.DataFrame:
    if features.empty:
        base = pd.DataFrame()
    else:
        work = apply_display_aliases_to_features(features, config)
        work["feature_name"] = work.apply(lambda row: str(row.get("display_name") or row.get("gene_name") or row.get("gene_id") or row.get("product") or "").strip(), axis=1)
        work["start"] = pd.to_numeric(work["start"], errors="coerce")
        work["end"] = pd.to_numeric(work["end"], errors="coerce")
        work = work.dropna(subset=["start", "end"])
        work = work[work["feature_name"] != ""].copy()
        priority = {"gene": 0, "tRNA": 1, "rRNA": 1, "CDS": 2, "transcript": 3, "exon": 4}
        work["_priority"] = work["feature_type"].map(priority).fillna(99).astype(int)
        work["_length"] = work["end"] - work["start"] + 1
        work = work.sort_values(["feature_name", "_priority", "_length", "start", "end"])
        work = work.drop_duplicates(subset=["feature_name"], keep="first")
        base = pd.DataFrame(
            {
                "feature_name": work["feature_name"],
                "feature_class": work.apply(feature_class, axis=1),
                "start": work["start"].astype(int),
                "end": work["end"].astype(int),
                "strand": work.get("strand", ""),
                "length_bp": (work["end"] - work["start"] + 1).astype(int),
            }
        )
    region_rows = []
    for item in config.get("analysis", {}).get("mt_regions", []) or []:
        start = item.get("start", "")
        end = item.get("end", "")
        try:
            length = abs(int(end) - int(start)) + 1
        except (TypeError, ValueError):
            length = ""
        region_rows.append(
            {
                "feature_name": item.get("name", ""),
                "feature_class": "configured region",
                "start": start,
                "end": end,
                "strand": ".",
                "length_bp": length,
            }
        )
    if region_rows:
        base = pd.concat([base, pd.DataFrame(region_rows)], ignore_index=True)
    if base.empty:
        return base
    return base.sort_values(["start", "end", "feature_name"]).reset_index(drop=True)


def card(label: str, value: str, help_text: str) -> str:
    return f"""
    <div class="metric-card">
      <div class="metric-label">{html.escape(label)}</div>
      <div class="metric-value">{html.escape(str(value))}</div>
      <div class="metric-help">{html.escape(help_text)}</div>
    </div>
    """


def section(title: str, text: str, body: str) -> str:
    return f"""
    <section>
      <div class="section-heading">
        <h2>{html.escape(title)}</h2>
        <p>{html.escape(text)}</p>
      </div>
      {body}
    </section>
    """


def svg_data_attribute(svg: str, name: str, fallback: str = "") -> str:
    match = re.search(rf"<svg\b[^>]*\bdata-{re.escape(name)}=(['\"])(.*?)\1", svg, flags=re.DOTALL)
    return html.unescape(match.group(2)) if match else fallback


def chord_target_id(prefix: str, value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return f"{prefix}__{token or 'all'}"


def circular_location_plot_panel(path: str, title: str, caption: str, link_prefix: str) -> str:
    aggregate = Path(path)
    interactive_svgs = sorted(aggregate.parent.glob(f"{aggregate.stem}__*__interactive.svg"))
    if not interactive_svgs:
        return ""
    subpanels = []
    for svg_path in interactive_svgs:
        svg = svg_path.read_text(encoding="utf-8", errors="ignore")
        fallback_group = svg_path.stem.removeprefix(f"{aggregate.stem}__").removesuffix("__interactive")
        group = svg_data_attribute(svg, "group", fallback_group.replace("_", " "))
        target_id = chord_target_id(f"interactive__{aggregate.stem}", group)
        static_pdf = svg_path.with_name(svg_path.name.replace("__interactive.svg", ".pdf"))
        subpanels.append(
            f"""
            <div class="plot-subpanel">
              <div class="plot-title-row">
                <h4>{html.escape(group)}</h4>
                <a class="plot-link" href="{html.escape(link_prefix)}/{html.escape(static_pdf.name)}">Open baseline PDF</a>
              </div>
              <div class="chord-controls" data-chord-controls data-target="{html.escape(target_id)}">
                <label class="slider-control">
                  <span>Minimum normalized support</span>
                  <input type="range" min="0" max="1000" step="1" value="0" data-support-slider>
                  <output data-support-output>All loaded calls</output>
                </label>
                <label class="observation-control">
                  <span>Minimum supporting observations</span>
                  <select data-observation-filter>
                    <option value="linked" data-linked-option>Auto</option>
                    <option value="1">1</option>
                    <option value="2">2</option>
                    <option value="5">5</option>
                    <option value="10">10</option>
                  </select>
                </label>
                <button type="button" data-reset-controls>Reset to PDF view</button>
                <div class="filter-status" data-filter-status></div>
              </div>
              <div class="control-note">The support slider uses a logarithmic scale. In Auto mode, the observation value reports the lowest raw count among calls passing the support filter; choose a number to enforce an additional raw-evidence cutoff. Moving the support slider returns this setting to Auto. Controls affect the HTML view only.</div>
              <div id="{html.escape(target_id)}" class="plot-svg">{svg}</div>
            </div>
            """
        )
    return f"""
    <article class="plot-panel circular_breakpoint_chords_all">
      <div class="plot-title-row">
        <h3>{html.escape(title)}</h3>
        <a class="plot-link" href="{html.escape(link_prefix)}/{html.escape(aggregate.name)}">Open multipage baseline PDF</a>
      </div>
      <p>{html.escape(caption)}</p>
      <div class="control-guidance"><strong>Using location controls</strong><p>Use normalized support to compare evidence across samples or datasets with different usable-read depths. Auto reports the lowest raw observation count among calls retained by the support slider; choose a fixed observation count only when an additional absolute-evidence requirement is intended.</p></div>
      {''.join(subpanels)}
    </article>
    """


def circular_comparison_plot_panel(path: str, title: str, caption: str, link_prefix: str) -> str:
    aggregate = Path(path)
    sidecar_svgs = sorted(aggregate.parent.glob(f"{aggregate.stem}__*.svg"))
    if not sidecar_svgs:
        return ""
    subpanels = []
    for svg_path in sidecar_svgs:
        svg = svg_path.read_text(encoding="utf-8", errors="ignore")
        left_group = svg_data_attribute(svg, "left-group", "left group")
        right_group = svg_data_attribute(svg, "right-group", "right group")
        target_id = chord_target_id(f"interactive__{aggregate.stem}", svg_path.stem)
        subpanels.append(
            f"""
            <div class="plot-subpanel comparison-plot-panel">
              <div class="plot-title-row">
                <h4>{html.escape(right_group)} compared with {html.escape(left_group)}</h4>
                <a class="plot-link" href="{html.escape(link_prefix)}/{html.escape(svg_path.with_suffix('.pdf').name)}">Open baseline PDF</a>
              </div>
              <p>Each chord represents one exact deletion comparison. Color shows the normalized mean-support difference ({html.escape(right_group)} minus {html.escape(left_group)}); the controls use statistics and support fields from the exact-deletion comparison table.</p>
              <div class="comparison-filter-block" data-comparison-controls data-target="{html.escape(target_id)}">
                <div class="comparison-primary-controls">
                  <label class="select-control">
                    <span>Comparison view</span>
                    <select data-comparison-preset>
                      <option value="all">All comparisons</option>
                      <option value="replicate-significant">Replicate-significant (BH q &le; 0.05)</option>
                      <option value="replicate-suggestive">Replicate p &le; 0.05 (exploratory)</option>
                      <option value="depth-significant">Read-depth enriched (BH q &le; 0.05; technical)</option>
                    </select>
                  </label>
                  <div class="filter-status" data-comparison-status></div>
                </div>
                <div class="preset-guidance" data-comparison-preset-guidance></div>
                <details class="advanced-comparison-filters">
                  <summary>Optional display refinements</summary>
                  <div class="comparison-controls">
                    <label class="slider-control">
                      <span>Minimum total supporting observations</span>
                      <input type="range" min="0" max="1000" step="1" value="0" data-comparison-observation-slider>
                      <output data-comparison-observation-output>&ge; 1</output>
                    </label>
                    <label class="slider-control">
                      <span>Minimum absolute normalized mean difference</span>
                      <input type="range" min="0" max="1000" step="1" value="0" data-comparison-difference-slider>
                      <output data-comparison-difference-output>All differences</output>
                    </label>
                    <label class="select-control">
                      <span>Direction</span>
                      <select data-comparison-direction>
                        <option value="both">Both directions</option>
                        <option value="right">Higher in {html.escape(right_group)}</option>
                        <option value="left">Higher in {html.escape(left_group)}</option>
                      </select>
                    </label>
                    <button type="button" data-reset-comparison-refinements>Reset refinements</button>
                  </div>
                </details>
              </div>
              <div id="{html.escape(target_id)}" class="plot-svg">{svg}</div>
            </div>
            """
        )
    return f"""
    <article class="plot-panel exact_deletion_comparison_chords">
      <div class="plot-title-row">
        <h3>{html.escape(title)}</h3>
        <a class="plot-link" href="{html.escape(link_prefix)}/{html.escape(aggregate.name)}">Open multipage baseline PDF</a>
      </div>
      <p>{html.escape(caption)}</p>
      <div class="control-guidance"><strong>Choosing a comparison view</strong><p>Use Replicate-significant (BH q &le; 0.05) for biological group conclusions. An empty view means no exact deletion passes that threshold. Replicate p &le; 0.05 is exploratory before multiple-testing correction. Read-depth enriched is technical read-count evidence and must not be reported as biological-replicate significance. Effect-size and observation refinements change only the display; they do not determine statistical significance.</p></div>
      {''.join(subpanels)}
    </article>
    """


def plot_panel(path: str, title: str, caption: str, link_prefix: str = "plots") -> str:
    if Path(path).stem == "circular_breakpoint_chords_all":
        panel = circular_location_plot_panel(path, title, caption, link_prefix)
        if panel:
            return panel
    if Path(path).stem == "exact_deletion_comparison_chords":
        panel = circular_comparison_plot_panel(path, title, caption, link_prefix)
        if panel:
            return panel
    svg_path = Path(path).with_suffix(".svg")
    pdf_name = html.escape(Path(path).name)
    panel_class = "plot-panel " + Path(path).stem.replace("-", "_")
    sidecar_svgs = sorted(Path(path).parent.glob(f"{Path(path).stem}__*.svg"))
    if Path(path).stem in {
        "deletion_rainfall_left_breakpoint",
        "deletion_rainfall_right_breakpoint",
        "deletion_rainfall_midpoint",
        "breakpoint_pair_support_map",
        "pooled_breakpoint_support_density",
        "pooled_breakpoint_support_density_capped",
    } and sidecar_svgs:
        previews = []
        for sidecar_svg in sidecar_svgs:
            group_label = sidecar_svg.stem.split("__", 1)[-1].replace("_", " ")
            sidecar_pdf = sidecar_svg.with_suffix(".pdf")
            sidecar_pdf_name = html.escape(sidecar_pdf.name)
            svg = sidecar_svg.read_text(encoding="utf-8", errors="ignore")
            previews.append(
                f"""
                <div class="plot-subpanel">
                  <div class="plot-title-row">
                    <h4>{html.escape(group_label)}</h4>
                    <a class="plot-link" href="{html.escape(link_prefix)}/{sidecar_pdf_name}">Open PDF</a>
                  </div>
                  <div class="plot-svg">{svg}</div>
                </div>
                """
            )
        return f"""
        <article class="{html.escape(panel_class)}">
          <div class="plot-title-row">
            <h3>{html.escape(title)}</h3>
          </div>
          <p>{html.escape(caption)}</p>
          {''.join(previews)}
        </article>
        """
    if svg_path.exists():
        svg = svg_path.read_text(encoding="utf-8", errors="ignore")
        preview = f'<div class="plot-svg">{svg}</div>'
    else:
        preview = '<div class="plot-missing">Plot preview not available.</div>'
    return f"""
    <article class="{html.escape(panel_class)}">
      <div class="plot-title-row">
        <h3>{html.escape(title)}</h3>
        <a class="plot-link" href="{html.escape(link_prefix)}/{pdf_name}">Open PDF</a>
      </div>
      <p>{html.escape(caption)}</p>
      {preview}
    </article>
    """


def reference_section(config: dict, features: pd.DataFrame) -> str:
    dataset = config.get("dataset", {})
    species = dataset.get("species", "")
    ref = config.get("references", {}).get(species, {})
    ref_dir = Path(config.get("project", {}).get("work_dir", "resources")) / "references" / str(species)
    genome_path = ref_dir / "genome.fa"
    annotation_path = ref_dir / "annotation.gtf"
    mt_path = ref_dir / "mt.fa"
    rows = [
        ("Species", species),
        ("Genome FASTA source", ref.get("genome_url") or ref.get("genome_path") or ""),
        ("Annotation GTF source", ref.get("annotation_url") or ref.get("annotation_path") or ""),
        ("Genome FASTA used", str(genome_path)),
        ("Genome FASTA SHA-256", sha256_file(genome_path)),
        ("Annotation GTF used", str(annotation_path)),
        ("Annotation GTF SHA-256", sha256_file(annotation_path)),
        ("Mitochondrial coordinate standard", ref.get("mt_reference_name", "")),
        ("Mitochondrial reference accession", ref.get("mt_reference_accession", "")),
        ("Extracted mtDNA FASTA used", str(mt_path)),
        ("Extracted mtDNA FASTA SHA-256", sha256_file(mt_path)),
        ("Mitochondrial length", ref.get("mt_length", "")),
        ("Mitochondrial contig names", ", ".join(ref.get("mt_contig_names", []))),
        ("Circular treatment", "normal and rotated mitochondrial references; alignment-directed circular deletion coordinates"),
    ]
    rows = [(field, value) for field, value in rows if value != ""]
    table = pd.DataFrame(rows, columns=["Field", "Value"])
    known = configured_deletion_target_table(config)
    known_sequences = known_sequence_search_table(config)
    regions = configured_region_table(config)
    replication_arcs = configured_replication_arc_table(config)
    extra = ""
    if not replication_arcs.empty:
        boundary_basis = str(ref.get("replication_arc_boundary_basis", "")).strip()
        basis_text = f'<p><strong>Boundary basis:</strong> {html.escape(boundary_basis)}</p>' if boundary_basis else ""
        extra += (
            '<h3>Configured Major And Minor Arcs</h3>'
            '<p>These reference-specific replication-arc labels annotate the alignment-directed deleted interval after calling. '
            'They do not select the deleted circular interval and are kept separate from affected-gene and configured-region annotations.</p>'
            + basis_text
            + table_html(replication_arcs, rows=10)
        )
    if not regions.empty:
        extra += '<h3>Configured Mitochondrial Regions</h3>' + table_html(regions, rows=100)
    if not known.empty:
        extra += '<h3>Configured Deletion Targets</h3><p>These coordinate targets are used to label matching remap calls. Rows can come from explicit known-deletion configuration or be inferred from coordinate-bearing configured sequence searches.</p>' + table_html(known, rows=100)
    if not known_sequences.empty:
        extra += '<h3>Configured Sequence Searches</h3>' + table_html(known_sequences, rows=100)
    return section(
        "Reference And Annotation",
        "Deletion coordinates and affected-feature labels depend directly on the mitochondrial reference, gene annotation, and any configured noncoding regions or deletion target intervals used here.",
        table_html(table, rows=20)
        + extra
        + '<h3>Features Shown In This Report</h3>'
        + table_html(simplified_feature_table(features, config), rows=300),
    )


def first_pass_selection_explanation(config: dict) -> str:
    mapping = config.get("mapping", {}) or {}
    mode = str(mapping.get("first_pass_read_selection", "whole_genome_mt_best"))
    if mode == "whole_genome_mt_best":
        return "Reads were competitively aligned with nuclear and mitochondrial references present together; reads with selected mitochondrial best evidence were retained for mitochondrial remapping."
    if mode == "nuclear_unmapped_reads":
        return "Reads not mapped to the nuclear-only reference were retained for mitochondrial remapping. This sensitivity mode can lose mitochondrial reads with NUMT-like nuclear alignments."
    return "Reads with configured mitochondrial evidence in the first-pass full-genome alignment were retained for mitochondrial remapping."


def assay_limitations(config: dict) -> str:
    dataset = config.get("dataset", {}) or {}
    molecule = str(dataset.get("molecule_type", "unknown")).strip().lower()
    assay = str(dataset.get("assay_type", "unknown")).strip().lower()
    parts = [
        "Reported calls are coordinate-focused deletion-like evidence from accepted split alignments, not automatic proof of biological mtDNA deletions."
    ]
    if molecule == "rna":
        parts.append(
            "RNA read support does not directly measure mtDNA heteroplasmy or genome copy fraction."
        )
    elif molecule == "dna":
        parts.append(
            "DNA-derived reads are closer to genome-molecule evidence, but local split-support fraction is not automatically a heteroplasmy estimate."
        )
    else:
        parts.append("Molecule type is not specified; DNA- versus RNA-specific biological interpretation is therefore limited.")
    if assay == "single_cell_rna_seq":
        parts.append(
            "Single-cell RNA-seq support can be affected by cell-level sparsity, amplification, barcode or UMI processing, and pooling. Unless cell identifiers are retained in the workflow inputs, deletion support is summarized at the configured sample or group level rather than as per-cell prevalence."
        )
    return " ".join(parts)


def potential_alternative_explanations(config: dict) -> pd.DataFrame:
    dataset = config.get("dataset", {}) or {}
    technology = str(dataset.get("read_technology", "unknown")).strip().lower()
    molecule = str(dataset.get("molecule_type", "unknown")).strip().lower()
    assay = str(dataset.get("assay_type", "unknown")).strip().lower()
    rows = [
        {
            "Applies to": "All datasets",
            "Potential alternative explanation": "Alternative or non-unique alignment placement",
            "How deletion-like evidence can arise": "NUMTs, repeats, or equivalent secondary placements can assign split segments to misleading mitochondrial coordinates.",
            "Relevant workflow control or limitation": "Competitive whole-genome read selection, MAPQ fields, and secondary-alignment provenance expose this uncertainty but cannot eliminate every ambiguous placement.",
        },
        {
            "Applies to": "All datasets",
            "Potential alternative explanation": "Chimeric source read or library molecule",
            "How deletion-like evidence can arise": "Fragments joined during library preparation or sequencing can contain an adjacency that was not present in the original mitochondrial molecule.",
            "Relevant workflow control or limitation": "Alignment chains keep paired mates separate and require one physical read sequence, but this does not prove that the source read itself is free of library chimeras.",
        },
        {
            "Applies to": "All datasets",
            "Potential alternative explanation": "Repeated observations from one original molecule",
            "How deletion-like evidence can arise": "PCR or assay amplification can make one source molecule appear as multiple supporting reads.",
            "Relevant workflow control or limitation": "Read and fragment identifiers are deduplicated where available; without validated molecular barcodes, supporting-read counts do not guarantee independent original molecules.",
        },
        {
            "Applies to": "Circular remapping",
            "Potential alternative explanation": "Reciprocal direction or reference-rotation disagreement",
            "How deletion-like evidence can arise": "Alternative split-alignment arrangements can support complementary circular arcs or appear in only one rotated coordinate system.",
            "Relevant workflow control or limitation": "The caller preserves directed adjacency, consolidates rotations, and reports or excludes reciprocal conflicts according to configuration. Crossing the coordinate origin alone is not evidence of an artifact.",
        },
    ]

    if technology == "illumina":
        rows.append(
            {
                "Applies to": "Illumina",
                "Potential alternative explanation": "Short or non-unique split anchors",
                "How deletion-like evidence can arise": "A short segment on either side of a junction can have several plausible placements, particularly near repeats or NUMTs.",
                "Relevant workflow control or limitation": "Anchor-length, aligned-fraction, and MAPQ settings constrain accepted segments. Mate distance alone is not used to call a deletion.",
            }
        )
    elif technology == "nanopore":
        rows.extend(
            [
                {
                    "Applies to": "Nanopore",
                    "Potential alternative explanation": "Basecalling errors and difficult sequence contexts",
                    "How deletion-like evidence can arise": "Errors, including those near homopolymers, can shift breakpoint placement or promote a split alignment.",
                    "Relevant workflow control or limitation": "Long anchors and alignment-quality fields support review, but breakpoint precision and validity still require read-level inspection for prioritized calls.",
                },
                {
                    "Applies to": "Nanopore",
                    "Potential alternative explanation": "Ligation products, internal adapters, or concatemers",
                    "How deletion-like evidence can arise": "Multiple source molecules can be sequenced as one apparent read and resemble a breakpoint-spanning molecule.",
                    "Relevant workflow control or limitation": "Supporting-read counts and alignment provenance are reported, but the caller does not by itself prove that a long read represents one original molecule.",
                },
            ]
        )
    else:
        rows.append(
            {
                "Applies to": "Unknown read technology",
                "Potential alternative explanation": "Technology-specific effects cannot be selected",
                "How deletion-like evidence can arise": "The relevant error profile, anchor-length limitations, and library artifacts depend on the sequencing technology.",
                "Relevant workflow control or limitation": "Set dataset.read_technology so the report can display the applicable interpretation guidance.",
            }
        )

    if molecule == "rna":
        rows.extend(
            [
                {
                    "Applies to": "RNA",
                    "Potential alternative explanation": "Mitochondrial transcript processing",
                    "How deletion-like evidence can arise": "Processed polycistronic transcripts can contain RNA adjacencies that are absent from the mitochondrial genome.",
                    "Relevant workflow control or limitation": "Configured transcript-compatible junctions are annotated and can be excluded, but that model does not cover every possible RNA-processing product.",
                },
                {
                    "Applies to": "RNA",
                    "Potential alternative explanation": "Reverse-transcription template switching",
                    "How deletion-like evidence can arise": "cDNA synthesis can join non-adjacent RNA segments into an apparent deletion junction.",
                    "Relevant workflow control or limitation": "This cannot be ruled out from alignment alone; independent support, replication, controls, and orthogonal DNA validation strengthen interpretation.",
                },
            ]
        )
    elif molecule == "dna":
        rows.append(
            {
                "Applies to": "DNA",
                "Potential alternative explanation": "NUMT-derived reads or DNA-library chimeras",
                "How deletion-like evidence can arise": "Nuclear mitochondrial sequence or PCR/ligation products can resemble a mitochondrial genomic breakpoint.",
                "Relevant workflow control or limitation": "Competitive mapping and alignment provenance reduce this risk, but prioritized calls still require molecule-level or orthogonal validation.",
            }
        )
    else:
        rows.append(
            {
                "Applies to": "Unknown molecule type",
                "Potential alternative explanation": "Molecule-specific effects cannot be selected",
                "How deletion-like evidence can arise": "RNA processing and reverse transcription differ from DNA-library and genome-molecule interpretations.",
                "Relevant workflow control or limitation": "Set dataset.molecule_type so the report can distinguish RNA- and DNA-specific alternatives.",
            }
        )

    if assay == "single_cell_rna_seq":
        rows.append(
            {
                "Applies to": "Single-cell RNA-seq",
                "Potential alternative explanation": "Amplification, barcode/UMI processing, sparsity, or pooling",
                "How deletion-like evidence can arise": "Technical amplification can exaggerate rare junctions, while barcode handling or pooled summaries can obscure whether support comes from distinct cells and molecules.",
                "Relevant workflow control or limitation": "Unless validated cell and molecule identifiers are retained, results are summarized at the configured sample or group level rather than interpreted as per-cell prevalence.",
            }
        )

    return pd.DataFrame(rows)


def assumptions_section(config: dict, clusters: pd.DataFrame, ambiguous_reads: pd.DataFrame, qc: pd.DataFrame) -> str:
    dataset = config.get("dataset", {}) or {}
    mt = config.get("mt_realign", {}) or {}
    junctions = config.get("junctions", {}) or {}
    warnings = []
    if junctions.get("arc_assignment", "alignment_directed") == "legacy_shortest_arc":
        warnings.append("Shortest-arc mode is active. Deleted intervals are not alignment-directed in this result.")
    secondary_calls = int(pd.to_numeric(qc.get("primary_calls_using_secondary_alignments", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
    mapq_zero_calls = int(pd.to_numeric(qc.get("primary_calls_with_min_mapq_zero", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
    if secondary_calls:
        noun = "row" if secondary_calls == 1 else "rows"
        verb = "uses" if secondary_calls == 1 else "use"
        warnings.append(f"{secondary_calls:,} retained deletion-supporting read {noun} {verb} at least one secondary alignment; these calls have increased alternative-placement uncertainty.")
    if mapq_zero_calls:
        noun = "row" if mapq_zero_calls == 1 else "rows"
        verb = "has" if mapq_zero_calls == 1 else "have"
        pronoun = "its" if mapq_zero_calls == 1 else "their"
        warnings.append(f"{mapq_zero_calls:,} retained deletion-supporting read {noun} {verb} minimum MAPQ 0; inspect {pronoun} read-level provenance before biological interpretation.")
    if str(dataset.get("read_technology", "unknown")).lower() == "unknown":
        warnings.append("Read technology is unknown.")
    if str(dataset.get("molecule_type", "unknown")).lower() == "unknown":
        warnings.append("Molecule type is unknown; DNA/RNA-specific interpretation is limited.")
    if not ambiguous_reads.empty:
        warnings.append(f"{len(ambiguous_reads):,} reciprocal-direction evidence rows were excluded from primary summaries.")
    if not clusters.empty and "rotation_agreement" in clusters.columns:
        single = int(clusters["rotation_agreement"].astype(str).eq("single_rotation").sum())
        if single:
            warnings.append(f"{single:,} exact deletions have support recorded from only one reference rotation.")
    warning_body = ""
    if warnings:
        warning_body = '<div class="notice"><strong>Interpretation warnings:</strong><ul>' + "".join(f"<li>{html.escape(item)}</li>" for item in warnings) + "</ul></div>"
    assumptions = pd.DataFrame(
        [
            ("Directed arc", "Split segments are ordered on the query and normalized to a forward-reference retained adjacency L -> R. On the plus strand the earlier query segment supplies L; on the minus strand the later query segment supplies L. The inferred deleted interval is the forward circular arc from retained base L to retained base R; the complementary arc is a different hypothesis."),
            ("Coordinate convention", "Breakpoints are retained flanking bases. Deleted size excludes both breakpoint bases."),
            ("Alignment chain", "Accepted split segments must come from one physical read sequence. SAM read1/read2 flags keep paired mates in separate chains; unpaired and long reads use their query name. This grouping does not establish that a chain is free from alternative placement or library artifacts."),
            ("Reference", "The configured mitochondrial reference, length, contig identity, and coordinate origin are assumed appropriate for the samples."),
            ("Mapping uniqueness", "Retained evidence is assumed not to be better explained by NUMTs, repeats, or equivalent secondary placements."),
            ("Clustering", "Directed breakpoints within the configured circular slop are assumed to represent the same exact coordinate-level event."),
            ("Read support", "Support measures detected split-alignment evidence; it is not automatically molecule frequency, heteroplasmy, viability, or proof of a biological deletion."),
            ("Detection limit", "Failure to detect a junction is not evidence that it is absent; sensitivity depends on coverage, molecule abundance, read length, mapping, and filters."),
        ],
        columns=["Assumption", "Meaning"],
    )
    arc_example = (
        "<h3>How The Deleted Arc Is Assigned</h3>"
        "<p>After query order is normalized by alignment strand, a split read with reference sequence ending at retained base L followed by sequence beginning at retained base R supports adjacency <code>L|R</code>. "
        "The inferred deleted bases are the forward circular interval between L and R. A read supporting <code>R|L</code> represents the complementary deletion model. "
        "Reference rotation changes coordinate origin but must not reverse this directed adjacency. Conflicting directions are retained as ambiguous evidence rather than resolved by interval length.</p>"
    )
    ambiguous_body = ""
    if not ambiguous_reads.empty:
        audit_columns = [
            col
            for col in [
                "sample",
                "read_id",
                "breakpoint_pair_id",
                "left_breakpoint",
                "right_breakpoint",
                "deleted_size",
                "rotation_name",
                "strand",
                "direction_status",
                "min_mapq",
            ]
            if col in ambiguous_reads.columns
        ]
        ambiguous_body = (
            "<h3>Excluded Ambiguous Direction Evidence</h3>"
            "<p>These rows support conflicting reciprocal directions for the same read and unordered breakpoint pair. They are retained for audit but excluded from primary summaries under the configured policy.</p>"
            + table_html(ambiguous_reads[audit_columns], rows=200)
        )
    return section(
        "Analysis Assumptions And Limitations",
        assay_limitations(config),
        warning_body
        + "<h3>Potential Alternative Explanations For Deletion-like Evidence</h3>"
        + "<p>This table lists technical artifacts, alignment ambiguity, and biological RNA phenomena that can resemble an mtDNA deletion in this dataset. These are possible explanations to evaluate, not findings that every reported call is artifactual.</p>"
        + table_html(potential_alternative_explanations(config), rows=30)
        + arc_example
        + "<h3>Core Assumptions</h3>"
        + table_html(assumptions, rows=20)
        + ambiguous_body,
    )


def method_section(config: dict, burden: pd.DataFrame) -> str:
    dataset = config.get("dataset", {}) or {}
    mt = config.get("mt_realign", {}) or {}
    mapping = config.get("mapping", {}) or {}
    qc = config.get("qc", {}) or {}
    junctions = config.get("junctions", {}) or {}
    analysis = config.get("analysis", {}) or {}
    settings = pd.DataFrame(
        [
            ("Result schema version", config.get("project", {}).get("result_schema_version", "unknown")),
            ("Read technology", dataset.get("read_technology", "unknown")),
            ("Molecule type", dataset.get("molecule_type", "unknown")),
            ("Assay type", dataset.get("assay_type", "unknown")),
            ("Library strategy", dataset.get("library_strategy", "unknown")),
            ("Sample source", config.get("samples", {}).get("source", "")),
            ("Read preparation", config.get("downloads", {}).get("method", "")),
            (
                "Adapter/QC trimming",
                "disabled by dataset configuration"
                if qc.get("trim_reads") is False
                else "fastp enabled by dataset configuration",
            ),
            ("First-pass aligner", mapping.get("first_pass_aligner", "star")),
            ("First-pass read selection", mapping.get("first_pass_read_selection", mt.get("input_strategy", "mt_evidence_reads"))),
            ("Short-read nuclear-depletion candidate", mapping.get("short_read_nuclear_depletion_candidate", "")),
            ("Long-read nuclear-depletion candidate", mapping.get("long_read_nuclear_depletion_candidate", "")),
            ("STAR chimSegmentMin", mapping.get("star_chimeric_options", {}).get("chimSegmentMin", "")),
            ("STAR alignIntronMax", mapping.get("star_chimeric_options", {}).get("alignIntronMax", "")),
            ("Full-genome mapQ filter", mapping.get("minimum_mapq_full_genome", "")),
            ("Remap-input selection mode", mapping.get("first_pass_read_selection", mt.get("input_strategy", "mt_evidence_reads"))),
            ("Keep ambiguous mitochondrial/nuclear reads", mapping.get("keep_ambiguous_mt_nuclear_reads", "")),
            ("Mitochondrial remap aligner", "minimap2"),
            ("Mitochondrial remap preset", mt.get("minimap2_preset", "sr")),
            ("Mitochondrial remap index options", mt.get("minimap2_index_extra", "")),
            ("Mitochondrial remap alignment options", mt.get("minimap2_extra", "")),
            ("Circular-coordinate handling", "normal plus rotated mitochondrial references, converted to standard coordinates while preserving directed junctions"),
            ("Reference rotations", ", ".join(str(item.get("name", "")) for item in mt.get("rotations", []))),
            ("Arc assignment", junctions.get("arc_assignment", "alignment_directed")),
            ("Alignment pairing mode", junctions.get("alignment_pairing_mode", "all_compatible")),
            (
                "Alignment-chain identity",
                "SAM read1/read2 flags separate mates when present; otherwise each query name identifies one physical read sequence; fragment support is deduplicated by base query name",
            ),
            ("Ambiguous direction policy", junctions.get("ambiguous_direction_policy", "exclude")),
            ("Include secondary remap alignments", mt.get("minimap2_include_secondary", True)),
            ("Include supplementary remap alignments", mt.get("minimap2_include_supplementary", True)),
            ("Minimum remap MAPQ", mt.get("minimap2_min_mapq", 0)),
            ("Minimum segment aligned fraction", mt.get("min_segment_aligned_fraction", "")),
            ("Maximum soft-clip fraction", mt.get("max_soft_clip_fraction", "")),
            ("Maximum query overlap", mt.get("max_query_overlap_bp", "")),
            ("Maximum query gap", mt.get("max_query_gap_bp", "")),
            ("Minimum deletion anchor length", junctions.get("min_anchor_length", "")),
            ("Minimum deletion size", junctions.get("min_deletion_size", "")),
            ("Maximum deletion size", junctions.get("max_deletion_size", "")),
            ("Breakpoint clustering slop", junctions.get("breakpoint_slop_bp", "")),
            ("Minimum split-read support", junctions.get("min_split_read_support", "")),
            ("Normalization denominator", analysis.get("normalization_denominator", normalization_mode(burden, config))),
            (
                "Expected transcript-compatible junctions",
                "excluded from deletion summaries; counts retained in QC"
                if junctions.get("exclude_expected_transcript_junctions", True)
                else "annotated and retained in deletion summaries",
            ),
        ],
        columns=["Setting", "Value"],
    )
    text = (
        "The workflow resolves sample metadata, prepares FASTQ inputs, applies read QC/trimming when configured, and uses a first-pass genome alignment to select mitochondrial-evidence reads for remapping. "
        + first_pass_selection_explanation(config)
        + " "
        "The retained reads are then remapped to mitochondrial references with minimap2 using normal and rotated coordinate systems so deletion breakpoints near the artificial linear boundary can be recovered. "
        "Split alignments are converted back to the original mitochondrial coordinate system, then query order is normalized by alignment strand into a forward-reference retained adjacency. Reciprocal directions remain distinct, and direction conflicts are handled by the configured ambiguity policy. "
        "Directed calls are consolidated across rotations, filtered, annotated against mitochondrial features, and summarized in group comparisons. "
        "Configured adjacent mitochondrial transcript pairs are labeled as transcript-compatible; when the exclusion setting is enabled, those reads are removed from deletion summaries but remain visible in QC counts. "
        "Configured deletion targets are used only for labeling and targeted summaries. They can come from explicit known-deletion entries or from coordinate-bearing configured sequence searches. "
        "Configured sequence searches are supplementary literal-read checks for named breakpoint sequences in the retained remap-input FASTQs; they do not replace the remapped deletion-calling analysis."
    )
    denominators = pd.DataFrame(
        [
            {
                "term": "Normalization denominator used here",
                "definition": normalization_definition(burden, config),
            },
            {
                "term": "Mitochondrial-evidence reads",
                "definition": "Reads retained after the first-pass genome assignment because their best or selected alignment evidence is mitochondrial. These reads are written to the remap-input FASTQs for mitochondrial remapping. They can be selected as the normalization denominator, but are not necessarily the denominator in this report.",
            },
            {
                "term": "Deletion-supporting reads",
                "definition": "Retained reads whose mitochondrial remap contains accepted split/supplementary alignment evidence for a directed retained adjacency and its inferred circular deletion interval after conversion, deduplication, filtering, and annotation.",
            },
            {
                "term": "Local reference-spanning reads",
                "definition": "Primary mitochondrial remap alignments that span the undeleted reference sequence near a reported breakpoint. These are used only for the local breakpoint reference-support columns, not for the main per-million normalization.",
            },
        ]
    )
    body = "<h3>Key Denominators</h3>" + table_html(denominators, rows=10) + "<h3>Run Settings</h3>" + table_html(settings, rows=40)
    return section("Workflow Method", text, body)


def known_deletion_rank_table(clusters: pd.DataFrame) -> pd.DataFrame:
    if clusters.empty or "known_deletion_label" not in clusters.columns:
        return pd.DataFrame()
    work = clusters.copy()
    work["known_deletion_label"] = work["known_deletion_label"].fillna("").astype(str)
    work = work[work["known_deletion_label"] != ""].copy()
    if work.empty:
        return pd.DataFrame()
    support = pd.to_numeric(work.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0)
    all_clusters = clusters.copy()
    all_clusters["_support"] = pd.to_numeric(all_clusters.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0)
    all_clusters = all_clusters.sort_values("_support", ascending=False, kind="mergesort").reset_index(drop=True)
    rank_by_id = {
        str(row.get("junction_id") or row.get("exact_deletion_id")): idx + 1
        for idx, row in all_clusters.iterrows()
    }
    work["_support"] = support
    rows = []
    for label, group in work.groupby("known_deletion_label", dropna=False):
        top = group.sort_values("_support", ascending=False, kind="mergesort").iloc[0]
        top_id = str(top.get("junction_id") or top.get("exact_deletion_id") or "")
        rows.append(
            {
                "configured_deletion_target": label,
                "matching_exact_deletion_calls": len(group),
                "total_supporting_reads_across_matches": int(group["_support"].sum()),
                "top_matching_exact_deletion": top_id,
                "top_matching_rank_among_all_exact_deletions": rank_by_id.get(top_id, ""),
                "top_matching_supporting_reads": int(top["_support"]),
                "top_left_breakpoint": top.get("left_breakpoint", ""),
                "top_right_breakpoint": top.get("right_breakpoint", ""),
                "top_deleted_size": top.get("deleted_size", ""),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["total_supporting_reads_across_matches", "top_matching_supporting_reads"],
        ascending=[False, False],
        kind="mergesort",
    )


def exact_deletion_support_read_links(clusters: pd.DataFrame, read_list_manifest: pd.DataFrame) -> dict[tuple[int, str], str]:
    if clusters.empty or read_list_manifest.empty or "total_supporting_reads" not in clusters.columns:
        return {}
    id_col = "exact_deletion_id" if "exact_deletion_id" in clusters.columns else "junction_id" if "junction_id" in clusters.columns else ""
    if not id_col or "exact_deletion_id" not in read_list_manifest.columns or "read_list_file" not in read_list_manifest.columns:
        return {}
    manifest = {
        str(row.get("exact_deletion_id", "")): str(row.get("read_list_file", ""))
        for _, row in read_list_manifest.iterrows()
        if str(row.get("read_list_file", "")).strip()
    }
    html_cells: dict[tuple[int, str], str] = {}
    for idx, row in clusters.iterrows():
        exact_id = str(row.get(id_col, ""))
        filename = manifest.get(exact_id)
        if not filename:
            continue
        count = html.escape(str(row.get("total_supporting_reads", "")))
        href = html.escape(f"read_lists/{filename}")
        title = html.escape(f"Open supporting read list for {exact_id}")
        html_cells[(idx, "total_supporting_reads")] = f'<a class="read-list-link" href="{href}" title="{title}">{count}</a>'
    return html_cells


def configured_target_read_links(
    rank_table: pd.DataFrame,
    clusters: pd.DataFrame,
    junction_reads: pd.DataFrame,
    out_dir: Path,
) -> dict[tuple[int, str], str]:
    if (
        rank_table.empty
        or clusters.empty
        or junction_reads.empty
        or "configured_deletion_target" not in rank_table.columns
        or "known_deletion_label" not in clusters.columns
        or "read_id" not in junction_reads.columns
    ):
        return {}
    cluster_id_col = "exact_deletion_id" if "exact_deletion_id" in clusters.columns else "junction_id" if "junction_id" in clusters.columns else ""
    read_id_col = "exact_deletion_id" if "exact_deletion_id" in junction_reads.columns else "junction_id" if "junction_id" in junction_reads.columns else ""
    if not cluster_id_col or not read_id_col:
        return {}

    out_dir.mkdir(parents=True, exist_ok=True)
    cluster_labels = clusters.copy()
    cluster_labels["known_deletion_label"] = cluster_labels["known_deletion_label"].fillna("").astype(str)
    cluster_labels[cluster_id_col] = cluster_labels[cluster_id_col].fillna("").astype(str)
    read_work = junction_reads.copy()
    read_work[read_id_col] = read_work[read_id_col].fillna("").astype(str)
    available_columns = [col for col in READ_LIST_COLUMNS if col in read_work.columns]
    html_cells: dict[tuple[int, str], str] = {}

    for idx, row in rank_table.iterrows():
        label = str(row.get("configured_deletion_target", "")).strip()
        if not label:
            continue
        exact_ids = set(
            cluster_labels.loc[cluster_labels["known_deletion_label"].eq(label), cluster_id_col]
            .dropna()
            .astype(str)
        )
        exact_ids.discard("")
        if not exact_ids:
            continue
        read_list = read_work[read_work[read_id_col].isin(exact_ids)][available_columns].copy()
        if read_list.empty:
            continue
        dedup_cols = [col for col in ["sample", "read_id", read_id_col] if col in read_list.columns]
        read_list = read_list.drop_duplicates(dedup_cols if dedup_cols else None, keep="first")
        sort_cols = [col for col in ["sample", read_id_col, "read_id"] if col in read_list.columns]
        if sort_cols:
            read_list = read_list.sort_values(sort_cols, kind="mergesort")
        filename = f"configured_target__{safe_sidecar_name(label)}.read_names.tsv"
        read_list.to_csv(out_dir / filename, sep="\t", index=False)
        count = html.escape(str(row.get("total_supporting_reads_across_matches", len(read_list))))
        href = html.escape(f"read_lists/{filename}")
        title = html.escape(f"Open supporting read rows for configured target {label}")
        html_cells[(idx, "total_supporting_reads_across_matches")] = f'<a class="read-list-link" href="{href}" title="{title}">{count}</a>'
    return html_cells


def evidence_streams_section(
    qc: pd.DataFrame,
    known_sequence_summary: pd.DataFrame,
    known_sequence_hits: pd.DataFrame,
    read_list_dir: Path,
    overlap: pd.DataFrame,
    overlap_html_cells: dict[tuple[int, str], str],
    config: dict | None = None,
    report_profile: str = "",
) -> str:
    dual_caller = bool((config or {}).get("quality", {}).get("short_read_rna_dual_caller", {}).get("enabled", False)) and bool(report_profile)
    rows = []
    if dual_caller:
        rows.extend(
            [
                {
                    "approach": "Combined canonical deletion evidence",
                    "role": "Primary profile result stream",
                    "input reads": "Unique physical observations retained by the selected quality profile",
                    "question answered": "Which circularly canonicalized deletion models pass this profile after cross-caller deduplication?",
                    "main outputs": "Exact deletions, caller-specific and combined support, gene-pair summaries, plots, matrices, and group comparisons",
                },
                {
                    "approach": "STAR chimeric alignment evidence",
                    "role": "Short-read RNA evidence source",
                    "input reads": "All reads aligned competitively to the full nuclear-plus-mitochondrial reference",
                    "question answered": "Which reads contain mitochondrial-to-mitochondrial chimeric alignments under the configured STAR geometry and quality filters?",
                    "main outputs": "STAR-supported canonical observations; full-genome nuclear competition remains available for quality annotation",
                },
                {
                    "approach": "Minimap2 mitochondrial remap evidence",
                    "role": "Targeted circular-remap evidence source",
                    "input reads": "Reads retained by the configured first-pass mitochondrial-evidence selection mode",
                    "question answered": "Which retained reads support deletion-like splits against normal and rotated mitochondrial references?",
                    "main outputs": "Minimap2-supported canonical observations and local breakpoint reference support",
                },
            ]
        )
    else:
        rows.append(
            {
                "approach": "Main remapping-based deletion caller",
                "role": "Primary discovery and reporting stream",
                "input reads": "Reads retained by the configured first-pass mitochondrial-evidence selection mode",
                "question answered": "Which reads support deletion-like split alignments after mitochondrial remapping, circular coordinate consolidation, filtering, and annotation?",
                "main outputs": "Exact deletion calls, affected-feature categories, size distributions, recurrence plots, group comparisons, and normalized burden tables",
            }
        )
    rows.extend([
        {
            "approach": "Configured literal sequence search",
            "role": "Supplementary targeted check",
            "input reads": "The same retained remap-input FASTQs",
            "question answered": "Do reads contain configured breakpoint-spanning motifs for named deletions?",
            "main outputs": "Read counts for configured sequence motifs only; this does not discover unconfigured deletions",
        },
    ])
    body = table_html(pd.DataFrame(rows), rows=10)

    if not known_sequence_summary.empty:
        body += "<h3>Configured Sequence Search Counts</h3>"
        body += (
            "<p>These counts come from literal sequence matching in the retained remap-input FASTQs. "
            "They are independent targeted checks for configured motifs, so their counts need not equal the remap caller counts. "
            "When read-level hits are available, click a value in the matching_reads column to open a TSV file with the matching read names.</p>"
        )
        linked_summary = default_sort_table(known_sequence_summary, "Configured Sequence Checks")
        html_cells = write_configured_sequence_read_lists(linked_summary, known_sequence_hits, read_list_dir)
        body += table_html(linked_summary, rows=300, html_cells=html_cells)

    if not overlap.empty:
        body += '<h3 id="sequence-remap-overlap">Configured Sequence Search Versus Remap Read Overlap</h3>'
        body += (
            "<p>This table compares the reads found by configured literal sequence searches with reads supporting remap-based exact deletions near the same configured target sites. "
            "Nearby remap calls are selected using the breakpoint and size tolerances shown in the table, because alignment placement can shift by a few bases even when reads support the same deletion model. "
            "The sequence_search_reads, remap_nearby_reads, shared_reads, sequence_only_reads, and remap_only_reads values link to TSV files containing the corresponding read names.</p>"
        )
        body += table_html(overlap, rows=300, html_cells=overlap_html_cells)

    if not qc.empty:
        before = pd.to_numeric(qc.get("clustered_split_read_alignments_before_transcript_filter", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()
        removed = pd.to_numeric(qc.get("expected_transcript_compatible_split_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()
        after = pd.to_numeric(qc.get("deletion_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()
        if before or removed or after:
            pct = (removed / before * 100.0) if before else 0.0
            body += "<h3>RNA-Processing Artifact Filter</h3>"
            body += (
                f"<p>The workflow labels configured expected mitochondrial transcript-compatible split reads before summarizing deletions. "
                f"In this run, {fmt_int(removed)} of {fmt_int(before)} split-read alignments ({fmt_float(pct)}%) were labeled transcript-compatible and excluded from deletion summaries, leaving {fmt_int(after)} deletion-supporting reads. "
                "This filter removes configured transcript patterns; it does not prove that every remaining split alignment is a biological mtDNA deletion.</p>"
            )
    return section(
        "How The Evidence Streams Differ",
        "This report keeps caller provenance and supplementary targeted sequence checks separate. When multiple callers support the same physical read and canonical deletion, combined support counts that observation once.",
        body,
    )


def quality_profile_section(config: dict, report_profile: str, clusters: pd.DataFrame, observations: pd.DataFrame) -> str:
    if not report_profile:
        return ""
    profiles = config.get("quality", {}).get("report_profiles", {}) or {}
    profile_config = profiles.get(report_profile, {})
    tiers = profile_config.get("include_tiers", []) if isinstance(profile_config, dict) else profile_config
    criteria = pd.DataFrame(
        [
            {"Field": "Report profile", "Value": report_profile},
            {"Field": "Included quality tiers", "Value": ", ".join(str(value) for value in tiers)},
            {"Field": "Exact deletion clusters retained", "Value": len(clusters)},
            {"Field": "Distinct physical observations retained", "Value": len(observations)},
            {"Field": "PCA behavior", "Value": "Matrices and PCA are rebuilt from this profile only; axes are not assumed equivalent across profiles"},
        ]
    )
    source_rows = []
    if not clusters.empty and "evidence_status" in clusters.columns:
        counts = clusters["evidence_status"].fillna("unknown").astype(str).value_counts()
        source_rows = [{"Evidence status": name, "Deletion clusters": int(value)} for name, value in counts.items()]
    return section(
        "Evidence Quality Profile",
        "This is one filtered view of shared canonical evidence. Stable deletion IDs and physical observations are not duplicated between report profiles.",
        table_html(criteria, rows=20)
        + ("<h3>Caller Evidence Composition</h3>" + table_html(pd.DataFrame(source_rows), rows=20) if source_rows else ""),
    )


def method_concordance_section(
    config: dict,
    report_profile: str,
    source_candidates: pd.DataFrame,
    clusters: pd.DataFrame,
    observations: pd.DataFrame,
) -> str:
    if not report_profile:
        return ""
    pieces = []
    if not observations.empty and "evidence_sources" in observations.columns:
        counts = observations["evidence_sources"].fillna("unknown").astype(str).value_counts().reset_index()
        counts.columns = ["retained evidence source(s)", "distinct physical observations"]
        pieces.append("<h3>Retained Observation Evidence</h3>" + table_html(counts, rows=20))
    if not clusters.empty and "evidence_status" in clusters.columns:
        counts = clusters["evidence_status"].fillna("unknown").astype(str).value_counts().reset_index()
        counts.columns = ["cluster evidence status", "exact deletion clusters"]
        pieces.append("<h3>Retained Cluster Evidence</h3>" + table_html(counts, rows=20))
    if not source_candidates.empty and {"evidence_source", "filter_status", "filter_reason"}.issubset(source_candidates.columns):
        source = source_candidates.copy()
        source["filter_reason"] = source["filter_reason"].fillna("").replace("", "passed source filters")
        disposition = (
            source.groupby(["evidence_source", "filter_status", "filter_reason"], dropna=False)
            .size()
            .reset_index(name="source records")
            .sort_values(["evidence_source", "filter_status", "source records"], ascending=[True, True, False])
        )
        pieces.append("<h3>Source-Specific Candidate Disposition</h3>" + table_html(disposition, rows=300))
    dual = bool(config.get("quality", {}).get("short_read_rna_dual_caller", {}).get("enabled", False))
    if dual:
        pieces.append(
            "<div class=\"notice\"><strong>STAR evidence boundary.</strong> "
            "The workflow parses mitochondrial-to-mitochondrial records from STAR Chimeric.out.junction, "
            "requires the configured query geometry and anchors, and by default requires distinct annotated mitochondrial gene anchors. "
            "Same-gene alignments, substantially overlapping gene anchors, configured expected transcript pairs, and configured size failures are excluded. "
            "This is direct STAR chimeric-record processing; STAR-Fusion is not run. Gene-pair labels are an RNA fusion-context aggregation and do not define exact deletion identity.</div>"
        )
    return section(
        "Caller Concordance And Source Filters",
        "Caller agreement is evidence provenance, not an additional molecule. A physical read or paired-end fragment supporting the same canonical deletion through more than one caller is counted once in combined support.",
        "".join(pieces),
    )


def configured_target_matches_panel(clusters: pd.DataFrame, junction_reads: pd.DataFrame, read_list_dir: Path) -> str:
    known_ranks = known_deletion_rank_table(clusters)
    if known_ranks.empty:
        return ""
    top_support = pd.to_numeric(clusters.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).max()
    rank_rows = []
    for _, row in known_ranks.iterrows():
        rank = row.get("top_matching_rank_among_all_exact_deletions", "")
        if rank == 1:
            interpretation = "This configured target is the strongest exact deletion signal in the remap caller."
        else:
            interpretation = (
                "This configured target is detected by the remap caller, but stronger exact deletion signals are also present."
            )
        rank_rows.append(
            {
                **row.to_dict(),
                "interpretation": interpretation,
                "top_overall_exact_deletion_supporting_reads": int(top_support),
            }
        )
    rank_display = pd.DataFrame(rank_rows).reset_index(drop=True)
    html_cells = configured_target_read_links(rank_display, clusters, junction_reads, read_list_dir)
    return result_table_panel(
        "Configured Target Matches",
        "Remap-called exact deletions are matched to configured deletion targets when the dataset configuration defines or implies target coordinates with breakpoint and size tolerances. A target can come from an explicit known-deletion entry or from a coordinate-bearing configured sequence search. Click total_supporting_reads_across_matches to open the read rows supporting all remap calls assigned to that target.",
        rank_display,
        rows=100,
        html_cells=html_cells,
        presorted=True,
    )


def sequence_remap_overlap_section(overlap: pd.DataFrame, html_cells: dict[tuple[int, str], str]) -> str:
    if overlap.empty:
        return ""
    return section(
        "Configured Sequence Search Versus Remap Read Overlap",
        "This cross-check compares reads found by configured literal sequence searches with reads supporting remap-based exact deletions near the configured target coordinates. Nearby remap calls are selected using the breakpoint and size tolerances shown in the table, because alignment ambiguity can shift a breakpoint without changing the underlying read-level evidence. Count values link to TSV files containing the corresponding read names.",
        table_html(overlap, rows=300, html_cells=html_cells),
    )


def experimental_design_section(
    samples: pd.DataFrame,
    burden: pd.DataFrame,
    group_col: str,
    config: dict | None = None,
    report_profile: str = "",
) -> str:
    if samples.empty:
        return ""
    pieces = []
    layouts = sorted(samples.get("layout", pd.Series(dtype=str)).dropna().astype(str).unique())
    if layouts:
        dual_caller = bool((config or {}).get("quality", {}).get("short_read_rna_dual_caller", {}).get("enabled", False)) and bool(report_profile)
        evidence_path = (
            "from individual STAR chimeric or mitochondrial-remap alignments"
            if dual_caller
            else "from individual split-read alignments after mitochondrial remapping"
        )
        if layouts == ["single"]:
            note = (
                "Read layout: single-end. "
                f"Each deletion-supporting observation is evaluated {evidence_path} within one read. "
                "No mate exists, so mate evidence is neither required nor used."
            )
        elif "paired" in layouts:
            note = "Read layout: " + ", ".join(layouts) + ". "
            note += f"Deletion evidence is evaluated {evidence_path}. Mates are not joined across reads to create a deletion call, and R1/R2 identifiers are collapsed to one fragment-level observation when they support the same canonical event. Mate-placement context is marked unavailable when the required alignment records are not present in the retained intermediates."
        else:
            note = "Read layout: " + ", ".join(layouts) + ". "
            note += f"Deletion evidence is evaluated {evidence_path}."
        pieces.append(f'<div class="notice">{html.escape(note)}</div>')
    if {"age", "treatment"}.issubset(samples.columns):
        counts = samples.groupby(["age", "treatment"], dropna=False).size().reset_index(name="sample_count")
        pieces.append("<h3>Factorial Sample Counts</h3>" + table_html(counts, rows=100))
        if not burden.empty and "deletion_support_per_million_mt_reads" in burden.columns:
            summary = (
                burden.groupby(["age", "treatment"], dropna=False)
                .agg(
                    sample_count=("sample", "count"),
                    mean_deletion_support_per_million_mt_reads=("deletion_support_per_million_mt_reads", "mean"),
                    mean_distinct_exact_deletions=("unique_exact_deletions", "mean"),
                )
                .reset_index()
            )
            pieces.append("<h3>Factorial Burden Summary</h3>" + table_html(summary, rows=100))
        text = (
            "This dataset has both age and treatment metadata, so the flat condition labels should be read as cells in a factorial design. "
            "The burden summary helps show whether the highest deletion burden is concentrated in one age-by-treatment cell, but interpretation should consider age, treatment, and their interaction rather than only pairwise condition contrasts."
        )
    elif group_col and group_col in samples.columns:
        text = f"The primary report grouping is `{group_col}`. Group counts are shown below so small-sample contrasts are visible before interpreting p-values."
        pieces.append("<h3>Primary Group Counts</h3>" + table_html(group_count_table(samples, group_col), rows=100))
    else:
        text = "No primary group column was configured for this dataset."
    return section("Experimental Design", text, "".join(pieces))


def circular_validation_section(config: dict, clusters: pd.DataFrame) -> str:
    mt = config.get("mt_realign", {}) or {}
    rotations = mt.get("rotations", []) or []
    rotation_names = ", ".join(str(item.get("name", "")) for item in rotations)
    rows = [
        {
            "check": "Rotated references configured",
            "result": rotation_names or "none",
            "interpretation": "Multiple coordinate starts reduce dependence on a single artificial linear mtDNA origin.",
        },
        {
            "check": "Coordinate consolidation",
            "result": "standard directed circular coordinates",
            "interpretation": "Calls from all rotations are reported on the original coordinate system without swapping reciprocal directions.",
        },
    ]
    if not clusters.empty:
        wraps = int((clusters.get("wraps_origin", pd.Series(dtype=str)).astype(str).str.lower() == "yes").sum())
        known = int(clusters.get("known_deletion_label", pd.Series(dtype=str)).fillna("").astype(str).ne("").sum())
        rows.extend(
            [
                {
                    "check": "Origin-spanning exact deletions",
                    "result": wraps,
                    "interpretation": "These alignment-directed deleted intervals pass through the configured coordinate origin.",
                },
                {
                    "check": "Configured deletion target matches",
                    "result": known,
                    "interpretation": "Matches are assigned only from dataset configuration, using breakpoint and size tolerances recorded in the reference section.",
                },
            ]
        )
    rows.append(
        {
            "check": "Rotation deduplication and directed reciprocal tests",
            "result": "covered by repository unit tests",
            "interpretation": "The test suite checks strand-normalized direction, rotated-coordinate conversion, separate reciprocal models, origin wrapping, and duplicate rotation support from the same read.",
        }
    )
    dual_caller = bool((config.get("quality", {}) or {}).get("short_read_rna_dual_caller", {}).get("enabled", False))
    source_text = (
        "Minimap2 mitochondrial-remap observations use normal and rotated references; STAR chimeric observations are converted directly from their full-genome mitochondrial coordinates. "
        if dual_caller
        else "The mitochondrial remap step uses normal and rotated references. "
    )
    return section(
        "Circular Coordinate Checks",
        source_text
        + "All observations are represented in one circular coordinate system while preserving directed retained adjacency. These checks summarize the circular-coordinate behavior visible in this run.",
        table_html(pd.DataFrame(rows), rows=20),
    )


def plot_sections(plots: dict[str, tuple[str, str]], link_prefix: str = "plots") -> str:
    return "".join(plot_panel(path, title, caption, link_prefix=link_prefix) for path, (title, caption) in plots.items())


def select_plots(plots: dict[str, tuple[str, str]], names: list[str]) -> dict[str, tuple[str, str]]:
    wanted = set(names)
    return {path: meta for path, meta in plots.items() if Path(path).name in wanted}


def stream_summary_cards(label: str, clusters: pd.DataFrame, burden: pd.DataFrame) -> str:
    n_deletions = len(clusters)
    total_support = int(pd.to_numeric(clusters.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) if not clusters.empty else 0
    mean_burden = ""
    if not burden.empty and "deletion_support_per_million_mt_reads" in burden.columns:
        mean_burden = f"{pd.to_numeric(burden['deletion_support_per_million_mt_reads'], errors='coerce').fillna(0).mean():.3g}"
    return "".join(
        [
            card(f"{label}: exact deletions", n_deletions, "Number of alignment-directed deletion intervals passing the clustering and support filters."),
            card(f"{label}: supporting reads", total_support, "Total split-read support assigned to reportable deletion intervals after rotation deduplication."),
            card(f"{label}: mean burden", mean_burden or "n/a", f"Average normalized deletion support across samples, {normalization_phrase(burden)}."),
        ]
    )


def reference_support_explanation(clusters: pd.DataFrame) -> str:
    required = {"left_reference_spanning_reads", "right_reference_spanning_reads", "reference_spanning_reads_min", "local_split_support_fraction"}
    if clusters.empty or not required.issubset(clusters.columns):
        return ""
    evidence_status = clusters.get("evidence_status", pd.Series(dtype=str)).fillna("").astype(str)
    evidence_sources = clusters.get("evidence_sources", pd.Series(dtype=str)).fillna("").astype(str)
    has_star_evidence = evidence_status.str.contains("star", case=False).any() or evidence_sources.str.contains("star", case=False).any()
    method_meaning = "How the local reference-spanning denominator was obtained. Remap-supported calls use primary-alignment depth and the maximum count across normal and rotated mitochondrial remaps."
    star_unavailable = ""
    if has_star_evidence:
        method_meaning += " STAR-only calls are marked unavailable because they have no comparable minimap2 split-support numerator."
        star_unavailable = "STAR-only exact deletions are marked unavailable for this metric because a STAR chimeric numerator and a minimap2 remap denominator would not be a like-for-like comparison. "
    definitions = pd.DataFrame(
        [
            {
                "column": "left_reference_spanning_reads",
                "meaning": "Primary-alignment spanning depth across the configured local window around the left breakpoint.",
            },
            {
                "column": "right_reference_spanning_reads",
                "meaning": "Primary-alignment spanning depth across the configured local window around the right breakpoint.",
            },
            {
                "column": "reference_spanning_reads_min",
                "meaning": "The smaller of the left and right reference-spanning counts. This conservative value prevents one high-coverage side from dominating the denominator.",
            },
            {
                "column": "local_split_support_fraction",
                "meaning": "split_supporting_reads / (split_supporting_reads + reference_spanning_reads_min). This is a local read-support fraction, not a DNA heteroplasmy estimate.",
            },
            {
                "column": "reference_support_method",
                "meaning": method_meaning,
            },
        ]
    )
    return (
        "<h3>Local Breakpoint Reference Support</h3>"
        "<p>For exact deletions with minimap2 remap evidence, the table includes local reference-spanning support at each breakpoint. "
        "These columns compare minimap2 deletion-supporting split reads with primary alignments that span the undeleted reference sequence near each breakpoint in the same mitochondrial remap stream. "
        "For each breakpoint, the workflow counts local spanning depth in the normal and rotated mitochondrial remaps and uses the larger value, rather than summing the two rotations. "
        f"{star_unavailable}"
        "This is a breakpoint-local reference-support denominator, separate from the configured per-million normalization denominator used in the main plots. "
        "The metric is most interpretable when reads are long enough and coverage exists at both breakpoint neighborhoods. For RNA data, it should not be interpreted as mtDNA heteroplasmy; for DNA data, it remains a local breakpoint-support summary unless the dataset and coverage assumptions justify a heteroplasmy interpretation.</p>"
        + table_html(definitions, rows=10)
    )


def exact_deletion_table_settings(config: dict) -> dict[str, object]:
    report = config.get("report", {}) or {}
    settings = report.get("exact_deletion_table", {}) or {}
    return {
        "min_total_supporting_reads": int(settings.get("min_total_supporting_reads", 0) or 0),
        "always_include_configured_targets": bool(settings.get("always_include_configured_targets", True)),
        "max_rows": int(settings.get("max_rows", 500) or 0),
    }


def configured_target_mask(df: pd.DataFrame) -> pd.Series:
    if "known_deletion_label" not in df.columns:
        return pd.Series(False, index=df.index)
    return df["known_deletion_label"].fillna("").astype(str).str.strip().ne("")


def exact_deletion_display_table(
    sorted_clusters: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, str]:
    if sorted_clusters.empty:
        return sorted_clusters, "No exact deletion rows are available."
    settings = exact_deletion_table_settings(config)
    min_support = int(settings["min_total_supporting_reads"])
    max_rows = int(settings["max_rows"])
    support = pd.to_numeric(sorted_clusters.get("total_supporting_reads", pd.Series(0, index=sorted_clusters.index)), errors="coerce").fillna(0)
    total_rows = len(sorted_clusters)
    total_support = int(support.sum())
    if max_rows <= 0 or total_rows <= max_rows:
        shown_support = int(support.sum())
        return (
            sorted_clusters.copy(),
            f"The embedded report table shows all {total_rows} exact deletions"
            f" ({shown_support:,} supporting reads) because the complete call set fits within the display cap."
            " The complete unfiltered exact-deletions table is delivered as tables/exact_deletions.tsv.",
        )
    include = support >= min_support
    target_rows = configured_target_mask(sorted_clusters)
    if settings["always_include_configured_targets"]:
        include = include | target_rows
    filtered = sorted_clusters.loc[include].copy()
    filtered_before_cap = len(filtered)
    if max_rows > 0 and len(filtered) > max_rows:
        filtered = filtered.head(max_rows).copy()
    target_retained = int(configured_target_mask(filtered).sum())
    target_total = int(target_rows.sum())
    shown_support = int(pd.to_numeric(filtered.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
    if min_support > 0:
        parts = [
            f"The embedded report table is filtered for readability: it shows exact deletions with at least {min_support} supporting reads",
        ]
    else:
        parts = [
            "The embedded report table is limited for readability: it prioritizes exact deletions by configured-target status and supporting-read count",
        ]
    if settings["always_include_configured_targets"]:
        parts.append(" and always includes configured deletion-target matches")
    cap_text = f", then caps the display at {max_rows} rows" if max_rows > 0 else ""
    parts.append(f"{cap_text}.")
    parts.append(
        f" Showing {len(filtered)} of {total_rows} exact deletions"
        f" ({shown_support:,} of {total_support:,} supporting reads)."
    )
    if max_rows > 0 and filtered_before_cap > max_rows:
        parts.append(f" {filtered_before_cap} rows passed the filter before the display cap.")
    if target_total:
        parts.append(f" Configured-target rows shown: {target_retained} of {target_total}.")
    parts.append(" The complete unfiltered exact-deletions table is delivered as tables/exact_deletions.tsv.")
    return filtered, "".join(parts)


def stream_result_section(
    section_id: str,
    title: str,
    text: str,
    config: dict,
    plots_body: str,
    diagnostics_body: str,
    clusters: pd.DataFrame,
    burden: pd.DataFrame,
    exact_comp: pd.DataFrame,
    affected_comp: pd.DataFrame,
    impact_comp: pd.DataFrame,
    size_tests: pd.DataFrame,
    size_bin_summary: pd.DataFrame,
    factorial_model_summary: pd.DataFrame,
    metadata_assoc: pd.DataFrame,
    per_gene: pd.DataFrame,
    junction_reads: pd.DataFrame,
    read_list_dir: Path,
    read_list_manifest: pd.DataFrame,
) -> str:
    sorted_clusters = default_sort_table(clusters, "Exact Deletion Calls")
    display_clusters, display_note = exact_deletion_display_table(sorted_clusters, config)
    exact_read_links = exact_deletion_support_read_links(display_clusters, read_list_manifest)
    return f"""
    <section id="{html.escape(section_id)}">
      <div class="section-heading">
        <h2>{html.escape(title)}</h2>
        <p>{html.escape(text)}</p>
      </div>
      <div class="metric-grid">{stream_summary_cards(title, clusters, burden)}</div>
      <h3>Core Plots</h3>
      {plots_body}
      <h3>Exploratory Ordination Plots</h3>
      {diagnostics_body}
      {reference_support_explanation(sorted_clusters)}
      <h3>Result Tables</h3>
      {configured_target_matches_panel(clusters, junction_reads, read_list_dir)}
      {result_table_panel("Deletion Burden By Sample", "One row per sample with normalized deletion burden, raw support, and sample-level QC fields.", burden)}
      {result_table_panel("Exact Deletion Group Comparisons", "Exact breakpoint intervals compared across the configured primary groups.", exact_comp)}
      {result_table_panel("Affected-Feature Group Comparisons", "Exact deletions collapsed to the mitochondrial features affected by the deleted interval.", affected_comp)}
      {result_table_panel("Exact Deletion Calls", "Alignment-directed exact deletion intervals with coordinates, complement diagnostics, direction and rotation status, support, sample recurrence, feature annotations, configured major/minor replication-arc context, and links from total_supporting_reads to read-level evidence files. " + display_note, display_clusters, html_cells=exact_read_links, presorted=True)}
      {result_table_panel("Collapsed Feature-Impact Comparisons", "Broad structural classes derived from affected-feature annotations.", impact_comp)}
      {result_table_panel("Deletion Size Distribution Tests", "Group-level tests comparing deletion-size distributions.", size_tests)}
      {result_table_panel("Deletion Size Bin Summary", "Group summaries for small, medium, and large deletion support.", size_bin_summary)}
      {optional_result_table_panel("Factorial Model Summary", "When age and treatment form a two-factor design, this table estimates age, treatment, and age-by-treatment terms for sample-level deletion outcomes.", factorial_model_summary)}
      {result_table_panel("Metadata Associations", "Associations between deletion burden and available metadata columns.", metadata_assoc)}
      {result_table_panel("Per-Gene Affected Burden", "Per-feature burden where each deletion can contribute to every mitochondrial feature it overlaps.", per_gene)}
    </section>
    """


def result_table_panel(
    title: str,
    description: str,
    df: pd.DataFrame,
    rows: int | None = None,
    html_cells: dict[tuple[int, str], str] | None = None,
    presorted: bool = False,
) -> str:
    sorted_df = df if presorted else default_sort_table(df, title)
    return f"""
    <article class="result-table-panel">
      <h3>{html.escape(title)}</h3>
      <p>{html.escape(description)}</p>
      {table_html(sorted_df, rows, html_cells=html_cells)}
    </article>
    """


def optional_result_table_panel(title: str, description: str, df: pd.DataFrame, rows: int | None = None) -> str:
    if df.empty:
        return ""
    return result_table_panel(title, description, df, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--title", required=True)
    parser.add_argument("--report-profile", default="")
    parser.add_argument("--run-results-dir", default="")
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--qc-summary", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--junction-reads", required=True)
    parser.add_argument("--source-candidates", default="")
    parser.add_argument("--ambiguous-reads", required=True)
    parser.add_argument("--burden", required=True)
    parser.add_argument("--exact-comparison", required=True)
    parser.add_argument("--affected-comparison", required=True)
    parser.add_argument("--impact-class-comparison", required=True)
    parser.add_argument("--size-tests", required=True)
    parser.add_argument("--size-bin-summary", required=True)
    parser.add_argument("--factorial-model-summary", required=True)
    parser.add_argument("--metadata-associations", required=True)
    parser.add_argument("--per-gene-burden", required=True)
    parser.add_argument("--plots", nargs="+", required=True)
    parser.add_argument("--known-sequence-summary", required=True)
    parser.add_argument("--known-sequence-hits", required=True)
    parser.add_argument("--read-list-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_report_config(args.config)
    samples = read_table(args.samples)
    features = read_table(args.features)
    qc = read_table(args.qc_summary)
    clusters = read_table(args.clusters)
    junction_reads = read_table(args.junction_reads)
    source_candidates = read_table(args.source_candidates) if args.source_candidates else pd.DataFrame()
    ambiguous_reads = read_table(args.ambiguous_reads)
    burden = read_table(args.burden)
    exact_comp = read_table(args.exact_comparison)
    affected_comp = read_table(args.affected_comparison)
    impact_comp = read_table(args.impact_class_comparison)
    size_tests = read_table(args.size_tests)
    size_bin_summary = read_table(args.size_bin_summary)
    factorial_model_summary = read_table(args.factorial_model_summary)
    metadata_assoc = read_table(args.metadata_associations)
    per_gene = read_table(args.per_gene_burden)
    known_sequence_summary = read_table(args.known_sequence_summary)
    known_sequence_hits = read_table(args.known_sequence_hits)
    read_list_dir = Path(args.read_list_manifest).parent
    read_list_manifest = write_exact_deletion_read_lists(junction_reads, read_list_dir)
    overlap_table, overlap_html_cells = sequence_remap_overlap_table(config, known_sequence_hits, junction_reads, read_list_dir)

    n_samples = len(samples)
    n_deletions = len(clusters)
    total_support = int(pd.to_numeric(clusters.get("total_supporting_reads", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) if not clusters.empty else 0
    group_col = config.get("dataset", {}).get("primary_group_column", "")
    group_columns = config.get("dataset", {}).get("group_columns", [])
    groups = ", ".join(sorted(samples[group_col].dropna().astype(str).unique())) if group_col in samples.columns else ""
    compact_samples = compact_sample_table(samples, group_columns)
    group_counts = group_count_table(samples, group_col)
    results_dir = Path(args.run_results_dir).resolve() if args.run_results_dir else Path(args.output).resolve().parent.parent
    read_prep = read_preparation_table(samples, results_dir)
    sample_note = ""
    if group_col in samples.columns and not group_counts.empty and group_counts["sample_count"].min() < 2:
        sample_note = (
            f'<div class="notice">At least one {html.escape(group_col)} group has fewer than two samples. '
            "Group-difference statistics and ordinations are descriptive for that group.</div>"
        )

    plot_meta = {
        "deletion_burden_by_sample.pdf": ("Total Deletion Burden", "Normalized deletion-supporting reads by sample and group. Colored dots are individual samples; the black diamond is the group mean. Confidence intervals are omitted when any group has fewer than three samples."),
        "unique_exact_deletions_by_sample.pdf": ("Distinct Exact Deletions", "Number of distinct alignment-directed exact deletions detected in each sample."),
        "deletion_burden_factorial_interaction.pdf": ("Deletion Burden: Age By Treatment", "Factorial interaction view for datasets with both age and treatment metadata. Points are samples; connected diamonds are treatment means at each age. This directly shows age effects, treatment effects, and possible age-by-treatment interaction patterns."),
        "unique_exact_deletions_factorial_interaction.pdf": ("Distinct Exact Deletions: Age By Treatment", "The same age-by-treatment interaction view applied to the number of distinct exact deletions per sample."),
        "deletion_size_distribution_unweighted.pdf": ("Deletion Size Distribution, Unweighted", "Distribution of deletion sizes where each deletion-supporting read contributes one count. This shows raw support and can be influenced by samples with more usable reads or more retained mitochondrial-evidence reads."),
        "deletion_size_distribution_support_weighted.pdf": ("Deletion Size Distribution, Normalized", f"Distribution of deletion sizes after each read is scaled to support {normalization_phrase(burden)} for its sample. This is better for comparing groups with different sequencing depth."),
        "deletion_size_distribution_support_weighted_log_y.pdf": ("Deletion Size Distribution, Normalized Log Scale", "The same normalized size distribution with a log y-axis, which keeps high-count small deletions visible while allowing lower-abundance large-deletion peaks to be seen."),
        "deletion_size_distribution_small.pdf": ("Small Deletions (<1 kb)", "Restricted normalized size distribution for small deletions. This separates the dense short-deletion range from larger events."),
        "deletion_size_distribution_medium.pdf": ("Medium Deletions (1-5 kb)", "Restricted normalized size distribution for medium deletions, where group-specific peaks can be hidden in the full-range plot."),
        "deletion_size_distribution_large.pdf": ("Large Deletions (>=5 kb)", "Restricted normalized size distribution for large deletions. This is useful for common-deletion-sized or paper-sized events."),
        "deletion_rainfall_left_breakpoint.pdf": ("Deletion Rainfall: Left Breakpoint", f"{rainfall_display_definition(config, burden)} Each point is one displayed exact deletion, placed by alignment-directed left breakpoint and deleted size on a log y-axis. Larger and brighter points have more normalized support. Marker numbers give unique support ranks within each group and match the right-breakpoint, circular-midpoint, and breakpoint-pair views; numbers are omitted when they do not fit. A cyan outline marks directed intervals spanning the coordinate origin. This is a display filter for the plot only; the exact-deletion table remains the complete call list subject to its own table-display settings."),
        "deletion_rainfall_right_breakpoint.pdf": ("Deletion Rainfall: Right Breakpoint", f"{rainfall_display_definition(config, burden)} This companion view uses the same displayed exact deletions, group-specific support ranks, and support scale as the left-breakpoint rainfall plot, but places each point by alignment-directed right breakpoint. Comparing left and right views helps identify fixed-endpoint patterns."),
        "deletion_rainfall_midpoint.pdf": ("Deletion Rainfall: Circular Midpoint", f"{rainfall_display_definition(config, burden)} This companion view uses the same displayed exact deletions and group-specific support ranks, placed by the circular midpoint of the deleted interval. Origin-spanning deletions are positioned by the midpoint along the deleted circular path rather than by a simple linear average."),
        "circular_breakpoint_chords_all.pdf": ("Circular Breakpoint Chords", f"{rainfall_display_definition(config, burden)} Each chord joins the alignment-directed left and right breakpoints of one exact deletion. The baseline PDF uses the rainfall display threshold and count cap; the HTML view loads every threshold-eligible call and provides normalized-support and raw-observation controls. Chord color uses the same normalized-support scale as the rainfall plots."),
        "exact_deletion_comparison_chords.pdf": ("Exact Deletion Group Comparison Chords", "Each chord represents one exact deletion in a configured group comparison. Chord color is the normalized mean-support difference between groups. The HTML views provide explicit replicate-level, exploratory, and technical read-depth presets plus optional display refinements."),
        "breakpoint_pair_support_map.pdf": ("Breakpoint-Pair Support Map", f"{rainfall_display_definition(config, burden)} Each point is one unique left/right breakpoint pair after applying the same display threshold and optional per-group cap as the rainfall plots. Marker numbers use the same group-specific support ranks as the three rainfall views. The x-axis is the left breakpoint; the y-axis is the right breakpoint, with origin-crossing right breakpoints shown above the horizontal genome-end line."),
        "pooled_breakpoint_support_density.pdf": ("Pooled Breakpoint Support Density", f"{rainfall_display_definition(config, burden)} This group-split view summarizes where deletion endpoints accumulate along the mitochondrial genome after pooling left and right breakpoints within each group. Stacked bars show binned support split by left versus right endpoint; the line is circular-smoothed total endpoint support."),
        "pooled_breakpoint_support_density_capped.pdf": ("Pooled Breakpoint Support Density: Capped Scale", f"{rainfall_display_definition(config, burden)} This is the same group-split endpoint-density view with the y-axis capped so smaller secondary breakpoint hotspots remain visible when one region dominates."),
        "affected_feature_support.pdf": ("Affected Features: Normalized Abundance", f"This bar chart compares affected-feature categories after normalizing each sample {normalization_phrase(burden)}. Use this as the main abundance view when groups have different sequencing depth or mitochondrial read recovery."),
        "affected_feature_counts.pdf": ("Affected Features: Raw Supporting Reads", "This uses the same affected-feature categories as the normalized plot, but shows raw supporting read counts. It can look similar when sample depths are similar; disagreement between this and the normalized plot suggests depth or recovery differences."),
        "affected_feature_proportions.pdf": ("Affected Features: Within-Group Percent", "This uses the same affected-feature categories again, but rescales each group to 100 percent. It asks whether the mix of affected features changes, independent of total deletion burden."),
        "feature_impact_classes.pdf": ("Collapsed Feature-Impact Classes", "This collapses detailed affected-feature labels into broad structural classes such as single-feature, two-feature, or mixed multi-feature deletions. It is less specific but easier to compare across datasets."),
        "per_gene_affected_burden.pdf": ("Per-Gene Affected Burden", "This abandons feature-pair categories and counts every gene or feature overlapped by deleted intervals. A single deletion can contribute to several genes, so this answers which mitochondrial features are most often touched overall."),
        "exact_deletion_recurrence.pdf": ("Exact Deletion Recurrence", "Top alignment-directed exact deletions ranked by supporting reads and shown as separate group bars. Row labels show directed left and right breakpoints, deleted size, and a shortened affected-feature label; full identifiers and annotations are in the exact-deletions table."),
        "exact_deletion_pca.pdf": ("Exact Deletion PCA", "Sample ordination using normalized exact-deletion support. Static labels and centroids are omitted so crowded datasets remain readable."),
        "exact_deletion_bray_curtis_mds.pdf": ("Exact Deletion Bray-Curtis MDS", "Distance-based sample ordination using exact-deletion profiles. Static labels and centroids are omitted so crowded datasets remain readable."),
        "affected_feature_pca.pdf": ("Affected-Feature PCA", "Sample ordination after exact deletions are summed to affected-feature categories. Static labels and centroids are omitted so crowded datasets remain readable."),
        "affected_feature_bray_curtis_mds.pdf": ("Affected-Feature Bray-Curtis MDS", "Distance-based sample ordination using affected-feature profiles. Static labels and centroids are omitted so crowded datasets remain readable."),
        "gene_pair_pca.pdf": ("Mitochondrial Gene-Pair PCA", "Short-read RNA sample ordination after retained canonical observations are aggregated by their STAR annotated mitochondrial gene pair when available, with breakpoint-flanking features used for other evidence sources. Gene-pair aggregation does not change exact deletion identity. The matrix is rebuilt for this report profile, so axes are not assumed equivalent to another profile."),
    }
    plot_lookup = {Path(path).name: path for path in args.plots}
    plots = {plot_lookup[name]: meta for name, meta in plot_meta.items() if name in plot_lookup}
    primary_plot_names = [
        "deletion_burden_by_sample.pdf",
        "unique_exact_deletions_by_sample.pdf",
        "deletion_size_distribution_support_weighted_log_y.pdf",
        "deletion_size_distribution_small.pdf",
        "deletion_size_distribution_medium.pdf",
        "deletion_size_distribution_large.pdf",
        "deletion_rainfall_left_breakpoint.pdf",
        "deletion_rainfall_right_breakpoint.pdf",
        "deletion_rainfall_midpoint.pdf",
        "circular_breakpoint_chords_all.pdf",
        "exact_deletion_comparison_chords.pdf",
        "breakpoint_pair_support_map.pdf",
        "pooled_breakpoint_support_density.pdf",
        "pooled_breakpoint_support_density_capped.pdf",
        "affected_feature_support.pdf",
        "affected_feature_proportions.pdf",
        "feature_impact_classes.pdf",
        "per_gene_affected_burden.pdf",
        "exact_deletion_recurrence.pdf",
    ]
    if {"age", "treatment"}.issubset(samples.columns):
        primary_plot_names[2:2] = [
            "deletion_burden_factorial_interaction.pdf",
            "unique_exact_deletions_factorial_interaction.pdf",
        ]
    secondary_plot_names = [
        "exact_deletion_pca.pdf",
        "exact_deletion_bray_curtis_mds.pdf",
        "affected_feature_pca.pdf",
        "affected_feature_bray_curtis_mds.pdf",
    ]
    dual_caller = bool((config.get("quality", {}) or {}).get("short_read_rna_dual_caller", {}).get("enabled", False))
    if dual_caller and "gene_pair_pca.pdf" in plot_lookup:
        primary_plot_names.insert(2, "gene_pair_pca.pdf")
    primary_plots = select_plots(plots, primary_plot_names)
    secondary_plots = select_plots(plots, secondary_plot_names)

    css = """
    body { margin: 0; font-family: Inter, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #1f2933; background: #f6f7f9; }
    header { background: #243447; color: white; padding: 34px 42px; }
    header h1 { margin: 0 0 8px; font-size: 30px; }
    header p { margin: 0; max-width: 980px; line-height: 1.5; color: #dce5ef; }
    nav { padding: 12px 42px; background: white; border-bottom: 1px solid #d8dee8; position: sticky; top: 0; z-index: 5; }
    nav a { margin-right: 18px; color: #285f8f; text-decoration: none; font-weight: 600; }
    main { box-sizing: border-box; padding: 24px 36px 48px; width: 100%; }
    section { background: white; border: 1px solid #d8dee8; border-radius: 8px; box-sizing: border-box; margin: 0 0 24px; overflow: hidden; padding: 24px; width: 100%; }
    .section-heading h2 { margin: 0 0 8px; font-size: 22px; }
    .section-heading p { margin: 0 0 18px; max-width: 1100px; line-height: 1.55; color: #52606d; }
    .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 14px; margin-bottom: 24px; }
    .metric-card { border-left: 5px solid #3b82a0; background: #f3f7fb; padding: 14px; border-radius: 6px; }
    .metric-label { font-size: 12px; text-transform: uppercase; letter-spacing: .04em; color: #52606d; font-weight: 700; }
    .metric-value { font-size: 24px; margin-top: 4px; font-weight: 760; }
    .metric-help { margin-top: 6px; color: #66788a; font-size: 13px; line-height: 1.35; }
    .plot-panel { border-top: 1px solid #e1e6ef; padding-top: 18px; margin-top: 18px; }
    .plot-title-row { display: flex; align-items: baseline; justify-content: space-between; gap: 14px; }
    .plot-title-row h3 { margin: 0; font-size: 18px; }
    .plot-panel p { color: #52606d; line-height: 1.45; }
    .plot-subpanel { border-top: 1px solid #edf1f6; margin-top: 18px; padding-top: 14px; }
    .plot-subpanel h4 { font-size: 16px; margin: 0; }
    .plot-link { background: #285f8f; color: white; border-radius: 5px; padding: 7px 10px; text-decoration: none; font-size: 13px; white-space: nowrap; }
    .plot-svg { max-width: 1100px; }
    .plot-svg svg { width: 100%; height: auto; max-height: 430px; }
    .deletion_rainfall_left_breakpoint .plot-svg,
    .deletion_rainfall_right_breakpoint .plot-svg,
    .deletion_rainfall_midpoint .plot-svg,
    .breakpoint_pair_support_map .plot-svg,
    .pooled_breakpoint_support_density .plot-svg,
    .pooled_breakpoint_support_density_capped .plot-svg { max-width: 1180px; }
    .deletion_rainfall_left_breakpoint .plot-svg svg,
    .deletion_rainfall_right_breakpoint .plot-svg svg,
    .deletion_rainfall_midpoint .plot-svg svg,
    .breakpoint_pair_support_map .plot-svg svg,
    .pooled_breakpoint_support_density .plot-svg svg,
    .pooled_breakpoint_support_density_capped .plot-svg svg { max-height: 900px; }
    .plot-missing, .empty, .notice { color: #8a4b16; background: #fff7ed; border: 1px solid #fed7aa; border-radius: 6px; padding: 12px; margin: 12px 0; }
    .table-wrap { background: white; border: 1px solid #d8dee8; border-radius: 6px; box-sizing: border-box; margin-top: 12px; max-height: 560px; max-width: 100%; overflow: auto; }
    .table-controls { align-items: center; background: #ffffff; border-bottom: 1px solid #d8dee8; display: flex; flex-wrap: wrap; gap: 10px; left: 0; margin: 0; padding: 8px; position: sticky; top: 0; z-index: 4; }
    .table-controls input { border: 1px solid #cbd5e1; border-radius: 5px; min-width: 260px; padding: 7px 9px; }
    .table-controls button { background: #edf2f7; border: 1px solid #cbd5e1; border-radius: 5px; color: #243447; cursor: pointer; padding: 7px 10px; }
    .table-controls button:disabled { color: #94a3b8; cursor: default; }
    .table-status { color: #52606d; font-size: 13px; }
    table.data-table { border-collapse: collapse; min-width: 100%; width: max-content; font-size: 12px; }
    table.data-table th { background: #eef2f7; color: #1f2933; cursor: pointer; position: sticky; top: 0; user-select: none; z-index: 2; }
    table.data-table th::after { background: #dbeafe; border: 1px solid #bfdbfe; border-radius: 999px; color: #285f8f; content: "sort"; display: inline-block; font-size: 9px; font-weight: 700; line-height: 1; margin-left: 6px; padding: 2px 4px; text-transform: uppercase; vertical-align: middle; }
    table.data-table th.sort-asc::after { background: #285f8f; border-color: #285f8f; color: white; content: "asc"; }
    table.data-table th.sort-desc::after { background: #285f8f; border-color: #285f8f; color: white; content: "desc"; }
    table.data-table th, table.data-table td { border: 1px solid #d8dee8; max-width: 300px; min-width: 92px; overflow-wrap: anywhere; padding: 5px 7px; text-align: left; vertical-align: top; white-space: normal; }
    table.data-table tr:nth-child(even) { background: #fafbfc; }
    .result-table-panel { border-top: 1px solid #e1e6ef; margin-top: 18px; padding-top: 16px; }
    .result-table-panel h3 { margin: 0 0 6px; }
    .result-table-panel p { color: #52606d; line-height: 1.45; margin: 0 0 8px; max-width: 1080px; }
    .result-table-panel table.data-table { font-size: 10.5px; }
    .result-table-panel table.data-table th, .result-table-panel table.data-table td { max-width: 260px; min-width: 108px; padding: 4px 6px; }
    .read-id-details summary { color: #285f8f; cursor: pointer; font-weight: 700; text-decoration: underline; text-underline-offset: 2px; }
    .read-id-details pre { background: #f8fafc; border: 1px solid #d8dee8; border-radius: 5px; color: #1f2933; font-size: 10px; line-height: 1.35; margin: 8px 0 0; max-height: 260px; min-width: 360px; overflow: auto; padding: 8px; white-space: pre; }
    .read-list-link { color: #285f8f; font-weight: 700; text-decoration: underline; text-underline-offset: 2px; }
    """
    js = """
    document.querySelectorAll('table.data-table').forEach((table) => {
      let page = 0;
      const pageSize = 40;
      const allRows = Array.from(table.tBodies[0]?.rows || []);
      let filteredRows = allRows.slice();
      const useControls = allRows.length > pageSize;
      let search = null, prev = null, next = null, status = null;
      if (useControls) {
        const controls = document.createElement('div');
        controls.className = 'table-controls';
        controls.innerHTML = '<input type="search" placeholder="Search this table"><button type="button" data-action="prev">Previous</button><button type="button" data-action="next">Next</button><span class="table-status"></span>';
        table.parentNode.insertBefore(controls, table);
        search = controls.querySelector('input');
        prev = controls.querySelector('[data-action="prev"]');
        next = controls.querySelector('[data-action="next"]');
        status = controls.querySelector('.table-status');
      }

      function renderTable() {
        if (!useControls) {
          allRows.forEach(row => row.style.display = '');
          return;
        }
        const totalPages = Math.max(1, Math.ceil(filteredRows.length / pageSize));
        page = Math.min(page, totalPages - 1);
        allRows.forEach(row => row.style.display = 'none');
        filteredRows.slice(page * pageSize, (page + 1) * pageSize).forEach(row => row.style.display = '');
        prev.disabled = page === 0;
        next.disabled = page >= totalPages - 1;
        const start = filteredRows.length ? page * pageSize + 1 : 0;
        const end = Math.min(filteredRows.length, (page + 1) * pageSize);
        status.textContent = `${start}-${end} of ${filteredRows.length} rows`;
      }

      if (useControls) {
        search.addEventListener('input', () => {
          const q = search.value.trim().toLowerCase();
          filteredRows = q ? allRows.filter(row => row.textContent.toLowerCase().includes(q)) : allRows.slice();
          page = 0;
          renderTable();
        });
        prev.addEventListener('click', () => { page -= 1; renderTable(); });
        next.addEventListener('click', () => { page += 1; renderTable(); });
      }

      table.querySelectorAll('th').forEach((th, index) => {
        th.title = 'Click to sort';
        th.addEventListener('click', () => {
          const tbody = table.tBodies[0];
          const asc = th.dataset.asc !== 'true';
          table.querySelectorAll('th').forEach(other => {
            if (other !== th) {
              other.dataset.asc = '';
              other.classList.remove('sort-asc', 'sort-desc');
            }
          });
          filteredRows.sort((a, b) => {
            const av = a.cells[index]?.textContent.trim() || '';
            const bv = b.cells[index]?.textContent.trim() || '';
            const an = Number(av), bn = Number(bv);
            const cmp = Number.isFinite(an) && Number.isFinite(bn) ? an - bn : av.localeCompare(bv);
            return asc ? cmp : -cmp;
          });
          th.dataset.asc = asc;
          th.classList.toggle('sort-asc', asc);
          th.classList.toggle('sort-desc', !asc);
          filteredRows.forEach(row => tbody.appendChild(row));
          allRows.filter(row => !filteredRows.includes(row)).forEach(row => tbody.appendChild(row));
          page = 0;
          renderTable();
        });
      });
      renderTable();
    });
    """

    report_asset_dir = Path(__file__).resolve().parents[1] / "report_assets"
    css += "\n" + (report_asset_dir / "circular_chords.css").read_text(encoding="utf-8")
    js += "\n" + (report_asset_dir / "circular_chords.js").read_text(encoding="utf-8")

    summary_cards = "".join(
        [
            card("Samples", n_samples, "Samples included after metadata resolution."),
            card("Groups", groups or "not configured", "Primary comparison labels from the dataset metadata."),
            card("Exact deletions", n_deletions, "Alignment-directed circular deletion intervals retained by this report view."),
            card("Distinct supporting observations", total_support, "Unique physical read or fragment observations assigned to retained exact deletions."),
        ]
    )

    html_out = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(args.title)} deletion report</title>
  <style>{css}</style>
</head>
<body>
  <header>
    <h1>{html.escape(args.title)}</h1>
    <p>Workflow report for coordinate-focused mitochondrial deletion evidence in sequencing data. The report identifies its evidence sources and quality profile; configured literal sequence searches remain supplementary checks for named breakpoint motifs.</p>
  </header>
  <nav>
    <a href="#samples">Samples</a>
    <a href="#read-prep">Read Prep</a>
    <a href="#design">Design</a>
    <a href="#reference">Reference</a>
    <a href="#method">Method</a>
    <a href="#assumptions">Assumptions</a>
    <a href="#evidence">Evidence Streams</a>
    <a href="#quality-profile">Quality Profile</a>
    <a href="#caller-concordance">Caller Concordance</a>
    <a href="#circular-checks">Circular Checks</a>
    <a href="#qc">QC</a>
    <a href="#remap-stream">Deletion Results</a>
    <a href="#sequence-remap-overlap">Sequence/Remap Overlap</a>
  </nav>
  <main>
    <div class="metric-grid">{summary_cards}</div>
    <section id="samples"><div class="section-heading"><h2>Sample Metadata</h2><p>These labels define the group and continuous-variable comparisons in the report. This compact table hides long FASTQ path columns so every included sample is visible; full input paths remain in the delivered configuration and generated metadata files.</p></div>{sample_note}<h3>Included Samples</h3>{table_html(compact_samples, 300)}<h3>Group Counts</h3>{table_html(group_counts, 100)}</section>
    <section id="read-prep"><div class="section-heading"><h2>Read Preparation</h2><p>This table reports the observed read layout, read cycles, mean read length before and after fastp, and whether the workflow ran trimming for each sample. Values are read from the per-sample QC files generated by this run.</p></div>{table_html(read_prep, None)}</section>
    <div id="design">{experimental_design_section(samples, burden, group_col, config, args.report_profile)}</div>
    <div id="reference">{reference_section(config, features)}</div>
    <div id="method">{method_section(config, burden)}</div>
    <div id="assumptions">{assumptions_section(config, clusters, ambiguous_reads, qc)}</div>
    <div id="evidence">{evidence_streams_section(qc, known_sequence_summary, known_sequence_hits, read_list_dir, overlap_table, overlap_html_cells, config, args.report_profile)}</div>
    <div id="quality-profile">{quality_profile_section(config, args.report_profile, clusters, junction_reads)}</div>
    <div id="caller-concordance">{method_concordance_section(config, args.report_profile, source_candidates, clusters, junction_reads)}</div>
    <div id="circular-checks">{circular_validation_section(config, clusters)}</div>
    <section id="qc"><div class="section-heading"><h2>Processing QC</h2><p>This table summarizes the first-pass read selection, mitochondrial remapping, and deletion-call denominators used by the report.</p></div>{table_html(qc, 300)}</section>
    {stream_result_section("remap-stream", "Canonical Deletion Evidence Results", f"These results contain the canonical observations retained by the {args.report_profile or 'configured'} report profile and normalize support {normalization_phrase(burden, config)}. Caller-specific support remains visible in the exact-deletion table. These are coordinate-focused deletion-like evidence; review the earlier Analysis Assumptions And Limitations section before biological interpretation.", config, plot_sections(primary_plots), plot_sections(secondary_plots), clusters, burden, exact_comp, affected_comp, impact_comp, size_tests, size_bin_summary, factorial_model_summary, metadata_assoc, per_gene, junction_reads, read_list_dir, read_list_manifest)}
  </main>
  <script>{js}</script>
</body>
</html>
"""
    ensure_parent(args.output)
    Path(args.output).write_text(html_out, encoding="utf-8")


if __name__ == "__main__":
    main()
