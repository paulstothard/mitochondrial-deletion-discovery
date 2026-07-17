import sys
import unittest
from argparse import Namespace
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from classify_mt_reads import classify_alignment, reverse_complement
from annotate_junctions import append_configured_regions, apply_feature_aliases, biological_features
from analyze_deletions import deduplicate_evidence_reads, expected_transcript_mask
from circular_deletions import (
    affected_feature_impact,
    configured_deletion_targets,
    directed_breakpoints,
    known_deletion_match,
    normalize_pos as normalize_rotated_pos,
    replication_arc_annotation,
)
from call_minimap2_deletions import alignment_chain_key, candidate_rows_from_chains, deletion_from_segments
from cluster_junctions import canonical_junction
from consolidate_deletions import cluster_rows, split_direction_conflicts
from make_rotated_mt_reference import rotate_sequence
from parse_split_alignments import aligned_bases_from_cigar, deletion_size, normalize_pos, parse_star_junction_line
from prepare_reads import normalized_layout
from plot_deletion_results import (
    MITOCHONDRIAL_FEATURE_COLORS,
    apply_cluster_coordinates,
    assign_group_support_ranks,
    draw_feature_track_axis,
    draw_location_feature_track,
    endpoint_density_figure,
    endpoint_density_hotspots,
    endpoint_density_pages,
    breakpoint_pair_support_map,
    burden_plot,
    category_bar,
    exact_recurrence,
    location_rainfall,
    location_features,
    mitochondrial_axis_bounds,
    pooled_endpoint_density,
    rainfall_point_sizes,
    rank_label_font_size,
    support_scale_limits,
    rainfall_support_limits,
    rainfall_y_axis_min,
    size_distribution,
    per_gene_plot,
    ordination,
    factorial_interaction_plot,
    prepare_location_plot_data,
    support_legend_values,
    support_size_legend_values,
    value_columns,
)
from estimate_breakpoint_reference_support import circular_window, window_covered
from make_deletion_report import (
    assay_limitations,
    assumptions_section,
    breakpoint_pair_plot_panel,
    circular_comparison_plot_panel,
    circular_location_plot_panel,
    endpoint_density_plot_panel,
    plot_panel,
    rainfall_location_plot_panel,
    configured_replication_arc_table,
    exact_deletion_display_table,
    exact_deletion_table_settings,
    exact_deletion_support_read_links,
    method_section,
    namespace_inline_svg,
    potential_alternative_explanations,
    sequence_remap_overlap_table,
    table_html,
    write_configured_sequence_read_lists,
    write_exact_deletion_read_lists,
)
from plot_circular_chords import (
    add_feature_ring,
    add_svg_chord_metadata,
    chord_dom_id,
    chord_path,
    circle_point,
    prepare_comparison_calls,
)
from resolve_samples import derive_age, derive_replicate, derive_treatment, make_sample_id, validate_dataset_inputs
from search_known_sequences import compiled_searches, match_multi_required, match_single, sample_from_fastq, scan_fastq_for_searches
from select_whole_genome_mt_from_sam import classify_group


class FakeRead:
    def __init__(self, reference_name="MT", mapq=60, unmapped=False, nh=1):
        self.reference_name = reference_name
        self.mapping_quality = mapq
        self.is_unmapped = unmapped
        self._nh = nh

    def has_tag(self, tag):
        return tag == "NH" and self._nh is not None

    def get_tag(self, tag):
        if tag != "NH" or self._nh is None:
            raise KeyError(tag)
        return self._nh


class FakeWholeGenomeRead:
    def __init__(
        self,
        reference_name="MT",
        mapq=60,
        unmapped=False,
        secondary=False,
        supplementary=False,
        query_length=100,
        query_alignment_length=100,
    ):
        self.reference_name = reference_name
        self.mapping_quality = mapq
        self.is_unmapped = unmapped
        self.is_secondary = secondary
        self.is_supplementary = supplementary
        self.query_length = query_length
        self.query_alignment_length = query_alignment_length
        self.query_sequence = "A" * query_length


class FakeChainRead:
    def __init__(self, query_name, read1=False, read2=False):
        self.query_name = query_name
        self.is_read1 = read1
        self.is_read2 = read2


def caller_segment(query_start, query_end, ref_start, ref_end, strand="+", read_id="readA", rotation="normal"):
    return {
        "read_id": read_id,
        "query_name": read_id,
        "alignment_chain_id": read_id,
        "query_start": query_start,
        "query_end": query_end,
        "anchor_length": query_end - query_start,
        "ref_start_raw": ref_start,
        "ref_end_raw": ref_end,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "strand": strand,
        "mapq": 60,
        "is_primary": "yes" if query_start == 0 else "no",
        "is_secondary": "no",
        "is_supplementary": "no" if query_start == 0 else "yes",
        "aligned_fraction": 0.5,
        "soft_clip_fraction": 0.5,
        "cigar": "50M50S" if query_start == 0 else "50S50M",
        "flag": 0,
        "alignment_score": 50,
        "edit_distance": 0,
        "sa_tag": "",
        "minimap2_type": "P" if query_start == 0 else "",
        "rotation_name": rotation,
    }


def caller_args(mt_length=1000):
    return Namespace(
        sample="s1",
        species="human",
        mt_length=mt_length,
        rotation_start=1,
        rotation_name="normal",
        min_deletion_size=1,
        max_deletion_size=mt_length - 2,
        max_query_overlap_bp=5,
        max_query_gap_bp=5,
        arc_assignment="alignment_directed",
        pairing_mode="all_compatible",
    )


class CoreTests(unittest.TestCase):
    def test_sample_resolution_checkpoint_has_one_tracked_output(self):
        snakefile = (Path(__file__).resolve().parents[1] / "Snakefile").read_text(encoding="utf-8")
        checkpoint = snakefile.split("checkpoint resolve_samples:", 1)[1].split("\n\nrule download_genome:", 1)[0]

        self.assertIn("sample_source=lambda wildcards", checkpoint)
        self.assertIn("    output:\n        samples=RESOLVED_SAMPLES,\n    params:", checkpoint)
        self.assertNotIn("\n        config=RESOLVED_CONFIG,", checkpoint)
        self.assertNotIn("\n        run_table=f", checkpoint)

    def test_nanopore_splice_preset_configures_index_seed_parameters(self):
        config_path = Path(__file__).resolve().parents[1] / "config" / "datasets" / "human_nanopore.yaml"
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

        self.assertEqual(config["mapping"]["first_pass_minimap2_preset"], "splice")
        self.assertEqual(config["mapping"]["first_pass_minimap2_index_extra"], "-k15 -w5")
        self.assertEqual(config["mt_realign"]["minimap2_preset"], "map-ont")
        self.assertEqual(config["mt_realign"]["minimap2_index_extra"], "-k15 -w10")

    def test_minimap2_index_paths_include_seed_profile(self):
        root = Path(__file__).resolve().parents[1]
        snakefile = (root / "Snakefile").read_text(encoding="utf-8")

        self.assertIn("minimap2_full_{FIRST_PASS_MINIMAP2_INDEX_TAG}.mmi", snakefile)
        self.assertIn("minimap2_mt_{MT_MINIMAP2_INDEX_TAG}_{{rotation}}.mmi", snakefile)
        self.assertIn("minimap2 -x {params.preset} {params.extra} -d", snakefile)

        for name in ("human_common_deletion", "human_bulkseq_matched_nanopore"):
            config = yaml.safe_load((root / "config" / "datasets" / f"{name}.yaml").read_text(encoding="utf-8"))
            self.assertEqual(config["mt_realign"]["minimap2_index_extra"], "-k21 -w11")

    def test_value_columns_excludes_normalization_denominators(self):
        matrix = pd.DataFrame(
            {
                "sample": ["s1"],
                "dataset": ["example"],
                "condition": ["treated"],
                "total_usable_reads": [1000],
                "reads_passed_to_minimap2": [900],
                "normalization_denominator": ["total_usable_reads"],
                "normalization_reads": [1000],
                "MT-CO1+MT-CYB": [5.0],
            }
        )
        samples = pd.DataFrame({"sample": ["s1"], "condition": ["treated"]})
        self.assertEqual(value_columns(matrix, samples), ["MT-CO1+MT-CYB"])
        self.assertEqual(value_columns(matrix, None), ["MT-CO1+MT-CYB"])

    def test_normalize_pos_wraps_padded_reference(self):
        self.assertEqual(normalize_pos(1001, 16313, 1000), 1)
        self.assertEqual(normalize_pos(1000, 16313, 1000), 16313)

    def test_normalize_pos_rotated_reference(self):
        self.assertEqual(normalize_pos(1, 100, rotation_start=51), 51)
        self.assertEqual(normalize_pos(50, 100, rotation_start=51), 100)
        self.assertEqual(normalize_pos(51, 100, rotation_start=51), 1)
        self.assertEqual(normalize_rotated_pos(51, 100, rotation_start=51), 1)

    def test_rotate_sequence(self):
        self.assertEqual(rotate_sequence("ABCDEFGHIJ", 1), "ABCDEFGHIJ")
        self.assertEqual(rotate_sequence("ABCDEFGHIJ", 6), "FGHIJABCDE")

    def test_deletion_size_linear_and_wrapped(self):
        self.assertEqual(deletion_size(100, 500, 16313), 399)
        self.assertEqual(deletion_size(16000, 200, 16313), 512)

    def test_directed_junction_preserves_reciprocal_breakpoint_models(self):
        mt_length = 16569
        direct = canonical_junction(
            {"left_breakpoint": 8472, "right_breakpoint": 13448, "deleted_size": 4975},
            mt_length,
        )
        reciprocal = canonical_junction(
            {"left_breakpoint": 13448, "right_breakpoint": 8472, "deleted_size": 11592},
            mt_length,
        )
        self.assertEqual((direct["left_breakpoint"], direct["right_breakpoint"], direct["deleted_size"]), (8472, 13448, 4975))
        self.assertEqual((reciprocal["left_breakpoint"], reciprocal["right_breakpoint"], reciprocal["deleted_size"]), (13448, 8472, 11592))
        self.assertEqual(reciprocal["canonical_orientation"], "alignment_directed")
        self.assertEqual(direct["breakpoint_pair_id"], reciprocal["breakpoint_pair_id"])

    def test_directed_breakpoints_do_not_choose_shorter_arc(self):
        directed = directed_breakpoints(100, 800, 1000)
        self.assertEqual(directed["deleted_size"], 699)
        self.assertEqual(directed["complement_deleted_size"], 299)
        self.assertEqual(directed["wraps_origin"], "no")

    def test_minimap2_plus_strand_query_order_assigns_directed_wrapping_arc(self):
        first = caller_segment(0, 50, 951, 980, "+")
        second = caller_segment(50, 100, 100, 149, "+")
        row = deletion_from_segments(first, second, caller_args())
        self.assertEqual((row["left_breakpoint"], row["right_breakpoint"]), (980, 100))
        self.assertEqual(row["deleted_size"], 119)
        self.assertEqual(row["wraps_origin"], "yes")
        self.assertEqual(row["arc_assignment_method"], "alignment_directed")

    def test_minimap2_minus_strand_query_order_is_normalized_to_forward_adjacency(self):
        first = caller_segment(0, 50, 100, 149, "-")
        second = caller_segment(50, 100, 951, 980, "-")
        row = deletion_from_segments(first, second, caller_args())
        self.assertEqual((row["left_breakpoint"], row["right_breakpoint"]), (980, 100))
        self.assertEqual(row["deleted_size"], 119)
        self.assertEqual(row["left_anchor_length"], 50)
        self.assertEqual(row["right_anchor_length"], 50)

    def test_minimap2_reverse_strand_discovers_unconfigured_linear_arc(self):
        first = caller_segment(0, 50, 800, 849, "-")
        second = caller_segment(50, 100, 100, 149, "-")
        row = deletion_from_segments(first, second, caller_args())
        self.assertEqual((row["left_breakpoint"], row["right_breakpoint"]), (149, 800))
        self.assertEqual(row["deleted_size"], 650)
        self.assertEqual(row["wraps_origin"], "no")

    def test_minimap2_reverse_strand_rat_cigar_layout_selects_small_linear_arc(self):
        first = caller_segment(0, 27, 4716, 4742, "-", read_id="SRR17380112.235192")
        second = caller_segment(20, 49, 4339, 4367, "-", read_id="SRR17380112.235192")
        args = caller_args(mt_length=16313)
        args.max_query_overlap_bp = 12
        row = deletion_from_segments(first, second, args)
        self.assertEqual((row["left_breakpoint"], row["right_breakpoint"]), (4367, 4716))
        self.assertEqual(row["deleted_size"], 348)
        self.assertEqual(row["complement_deleted_size"], 15963)

    def test_common_deletion_style_overlapping_split_chain_is_retained(self):
        first = caller_segment(0, 64, 8419, 8482, "+", read_id="fragmentA")
        second = caller_segment(54, 150, 13450, 13545, "+", read_id="fragmentA")
        args = caller_args(mt_length=16569)
        args.max_query_overlap_bp = 10
        row = deletion_from_segments(first, second, args)
        self.assertEqual((row["left_breakpoint"], row["right_breakpoint"]), (8482, 13450))
        self.assertEqual(row["deleted_size"], 4967)

    def test_alignment_chain_key_separates_mates_with_shared_bam_query_name(self):
        read1 = FakeChainRead("fragmentA", read1=True)
        read2 = FakeChainRead("fragmentA", read2=True)
        self.assertEqual(alignment_chain_key(read1), "fragmentA/1")
        self.assertEqual(alignment_chain_key(read2), "fragmentA/2")
        self.assertNotEqual(alignment_chain_key(read1), alignment_chain_key(read2))

    def test_alignment_chain_key_preserves_unpaired_read_identity(self):
        read = FakeChainRead("nanopore_read_42")
        self.assertEqual(alignment_chain_key(read), "nanopore_read_42")

    def test_candidate_generation_keeps_mates_separate_and_deduplicates_fragment_support(self):
        mate1_first = caller_segment(0, 64, 8419, 8482, "+", read_id="fragmentA")
        mate1_second = caller_segment(54, 150, 13450, 13545, "+", read_id="fragmentA")
        mate2_first = caller_segment(0, 96, 13450, 13545, "-", read_id="fragmentA")
        mate2_second = caller_segment(86, 150, 8419, 8482, "-", read_id="fragmentA")
        for segment in (mate1_first, mate1_second):
            segment["alignment_chain_id"] = "fragmentA/1"
        for segment in (mate2_first, mate2_second):
            segment["alignment_chain_id"] = "fragmentA/2"
        args = caller_args(mt_length=16569)
        args.max_query_overlap_bp = 10
        rows = candidate_rows_from_chains(
            {"fragmentA/1": [mate1_first, mate1_second], "fragmentA/2": [mate2_first, mate2_second]},
            args,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["read_id"], "fragmentA")
        self.assertEqual((rows[0]["left_breakpoint"], rows[0]["right_breakpoint"]), (8482, 13450))

    def test_candidate_generation_never_joins_single_segments_from_opposite_mates(self):
        mate1 = caller_segment(0, 50, 1200, 1249, "+", read_id="fragmentA")
        mate2 = caller_segment(50, 100, 9200, 9249, "+", read_id="fragmentA")
        mate1["alignment_chain_id"] = "fragmentA/1"
        mate2["alignment_chain_id"] = "fragmentA/2"
        rows = candidate_rows_from_chains(
            {"fragmentA/1": [mate1], "fragmentA/2": [mate2]},
            caller_args(mt_length=16569),
        )
        self.assertEqual(rows, [])

    def test_candidate_generation_discovers_arbitrary_non_target_junction(self):
        first = caller_segment(0, 50, 1200, 1249, "+", read_id="discovery_read")
        second = caller_segment(50, 100, 9200, 9249, "+", read_id="discovery_read")
        rows = candidate_rows_from_chains({"discovery_read": [first, second]}, caller_args(mt_length=16569))
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["left_breakpoint"], rows[0]["right_breakpoint"]), (1249, 9200))
        self.assertEqual(rows[0]["deleted_size"], 7950)

    def test_reciprocal_directed_calls_do_not_cluster_together(self):
        rows = [
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "100", "right_breakpoint": "800", "deleted_size": "699", "rotation_name": "normal"},
            {"sample": "s1", "species": "human", "read_id": "readB", "left_breakpoint": "800", "right_breakpoint": "100", "deleted_size": "299", "rotation_name": "normal"},
        ]
        _, clusters, _ = cluster_rows(rows, slop=5, min_support=1, mt_length=1000)
        self.assertEqual(len(clusters), 2)
        self.assertEqual({(row["left_breakpoint"], row["right_breakpoint"]) for row in clusters}, {(100, 800), (800, 100)})

    def test_same_read_reciprocal_rotations_are_ambiguous(self):
        rows = [
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "100", "right_breakpoint": "800", "deleted_size": "699", "rotation_name": "normal"},
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "800", "right_breakpoint": "100", "deleted_size": "299", "rotation_name": "half"},
        ]
        accepted, ambiguous = split_direction_conflicts(rows)
        self.assertEqual(accepted, [])
        self.assertEqual(len(ambiguous), 2)
        self.assertTrue(all(row["direction_status"] == "ambiguous_across_rotations" for row in ambiguous))

    def test_same_read_reciprocal_rotations_are_ambiguous_with_breakpoint_slop(self):
        rows = [
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "100", "right_breakpoint": "800", "deleted_size": "699", "rotation_name": "normal"},
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "803", "right_breakpoint": "98", "deleted_size": "294", "rotation_name": "half"},
        ]
        accepted, ambiguous = split_direction_conflicts(rows, mt_length=1000, slop=5)
        self.assertEqual(accepted, [])
        self.assertEqual(len(ambiguous), 2)

    def test_circular_breakpoint_slop_clusters_across_coordinate_origin(self):
        rows = [
            {"sample": "s1", "species": "human", "read_id": "readA", "left_breakpoint": "999", "right_breakpoint": "400", "deleted_size": "400", "rotation_name": "normal"},
            {"sample": "s1", "species": "human", "read_id": "readB", "left_breakpoint": "1", "right_breakpoint": "401", "deleted_size": "399", "rotation_name": "half"},
        ]
        _, clusters, _ = cluster_rows(rows, slop=3, min_support=1, mt_length=1000)
        self.assertEqual(len(clusters), 1)

    def test_known_sequence_single_search_matches_reverse_complement(self):
        deletion = {
            "search_sequences": [
                {
                    "id": "junction_sequence",
                    "sequence": "AAAAACCCCC",
                    "reverse_complement": "GGGGGTTTTT",
                }
            ]
        }
        matched, ids, orientations = match_single("TTTGGGGGTTTTTAAA", deletion)
        self.assertTrue(matched)
        self.assertEqual(ids, ["junction_sequence"])
        self.assertEqual(orientations, ["reverse_complement"])

    def test_known_sequence_multi_search_requires_all_sequences(self):
        deletion = {
            "search_sequences": [
                {"id": "primary", "sequence": "AAA", "reverse_complement": "TTT"},
                {"id": "support", "sequence": "CCC", "reverse_complement": "GGG"},
            ]
        }
        self.assertTrue(match_multi_required("NNNAAANNNGGGNNN", deletion)[0])
        self.assertFalse(match_multi_required("NNNAAANNN", deletion)[0])

    def test_known_sequence_sample_name_from_mt_evidence_fastq(self):
        self.assertEqual(sample_from_fastq("results/d/mt_reads/KSS-95.mt_evidence.fastq.gz", "_R1.fastq.gz"), "KSS-95")

    def test_known_sequence_scan_records_matching_read_ids(self):
        import gzip
        import tempfile

        deletion = {
            "id": "del1",
            "name": "configured deletion",
            "search_strategy": {"type": "single_sequence"},
            "search_sequences": [{"id": "junction", "sequence": "AACCGG", "reverse_complement": "CCGGTT"}],
        }
        with tempfile.NamedTemporaryFile(suffix=".fastq.gz") as handle:
            with gzip.open(handle.name, "wt", encoding="utf-8") as fastq:
                fastq.write("@readA cell=cell1\nTTTAACCGGTTT\n+\nFFFFFFFFFFFF\n")
                fastq.write("@readB cell=cell2\nTTTTTTTTTTTT\n+\nFFFFFFFFFFFF\n")
            counts, examined, hits = scan_fastq_for_searches(handle.name, compiled_searches([deletion]), sample="sample1", mate="R1")
        self.assertEqual(examined, 2)
        self.assertEqual(counts["del1"], 1)
        self.assertEqual(hits[0]["read_id"], "readA cell=cell1")
        self.assertEqual(hits[0]["sample"], "sample1")

    def test_consolidate_deletions_deduplicates_rotations_by_read(self):
        rows = [
            {
                "sample": "s1",
                "species": "human",
                "read_id": "readA/1",
                "left_breakpoint": "8472",
                "right_breakpoint": "13448",
                "deleted_size": "4975",
                "rotation_name": "normal",
            },
            {
                "sample": "s1",
                "species": "human",
                "read_id": "readA/1",
                "left_breakpoint": "8473",
                "right_breakpoint": "13449",
                "deleted_size": "4975",
                "rotation_name": "half",
            },
        ]
        all_rows, clusters, id_rows = cluster_rows(rows, slop=10, min_support=1, mt_length=16569)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["total_supporting_reads"], 1)
        self.assertEqual(len(all_rows), 1)
        self.assertEqual(len(id_rows), 1)

    def test_consolidate_deletions_recomputes_size_from_representative_breakpoints(self):
        rows = [
            {
                "sample": "s1",
                "species": "human",
                "read_id": "readA",
                "left_breakpoint": "100",
                "right_breakpoint": "201",
                "deleted_size": "100",
                "rotation_name": "normal",
            },
            {
                "sample": "s1",
                "species": "human",
                "read_id": "readB",
                "left_breakpoint": "101",
                "right_breakpoint": "202",
                "deleted_size": "100",
                "rotation_name": "half",
            },
        ]
        _, clusters, _ = cluster_rows(rows, slop=2, min_support=1, mt_length=1000)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["left_breakpoint"], 100)
        self.assertEqual(clusters[0]["right_breakpoint"], 202)
        self.assertEqual(clusters[0]["deleted_size"], 101)

    def test_location_plot_data_uses_cluster_representative_coordinates(self):
        reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "junction_id": "j1", "left_breakpoint": 100, "right_breakpoint": 200, "deleted_size": 99},
                {"sample": "s1", "read_id": "readB", "junction_id": "j1", "left_breakpoint": 110, "right_breakpoint": 210, "deleted_size": 99},
            ]
        )
        samples = pd.DataFrame({"sample": ["s1"], "condition": ["treated"], "normalization_reads": [1000]})
        clusters = pd.DataFrame({"junction_id": ["j1"], "left_breakpoint": [105], "right_breakpoint": [205], "deleted_size": [99]})
        grouped, groups, _ = prepare_location_plot_data(reads, samples, "condition", clusters=clusters, mt_length=1000)
        self.assertEqual(groups, ["treated"])
        self.assertEqual(len(grouped), 1)
        self.assertEqual(int(grouped.loc[0, "left_breakpoint"]), 105)
        self.assertEqual(int(grouped.loc[0, "right_breakpoint"]), 205)
        self.assertEqual(int(grouped.loc[0, "deleted_size"]), 99)
        self.assertEqual(int(grouped.loc[0, "supporting_reads"]), 2)
        self.assertEqual(int(grouped.loc[0, "_support_rank"]), 1)

    def test_location_support_ranks_are_unique_stable_and_group_specific(self):
        data = pd.DataFrame(
            [
                {"_plot_group": "a", "_plot_support": 10.0, "exact_deletion_id": "mtDel_b", "left_breakpoint": 20, "right_breakpoint": 40, "deleted_size": 19},
                {"_plot_group": "a", "_plot_support": 5.0, "exact_deletion_id": "mtDel_c", "left_breakpoint": 30, "right_breakpoint": 60, "deleted_size": 29},
                {"_plot_group": "a", "_plot_support": 10.0, "exact_deletion_id": "mtDel_a", "left_breakpoint": 10, "right_breakpoint": 30, "deleted_size": 19},
                {"_plot_group": "b", "_plot_support": 2.0, "exact_deletion_id": "mtDel_d", "left_breakpoint": 40, "right_breakpoint": 80, "deleted_size": 39},
            ]
        )
        ranked = assign_group_support_ranks(data).set_index("exact_deletion_id")
        self.assertEqual(int(ranked.loc["mtDel_a", "_support_rank"]), 1)
        self.assertEqual(int(ranked.loc["mtDel_b", "_support_rank"]), 2)
        self.assertEqual(int(ranked.loc["mtDel_c", "_support_rank"]), 3)
        self.assertEqual(int(ranked.loc["mtDel_d", "_support_rank"]), 1)

    def test_rank_label_font_size_uses_small_range_and_omits_tiny_markers(self):
        self.assertEqual(rank_label_font_size(600.0, 1), 6.0)
        fitted = rank_label_font_size(100.0, 12)
        self.assertIsNotNone(fitted)
        self.assertGreaterEqual(fitted, 4.0)
        self.assertLessEqual(fitted, 6.0)
        self.assertIsNone(rank_label_font_size(20.0, 1))

    def test_cluster_coordinate_adapter_recomputes_size_from_representative_breakpoints(self):
        reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "junction_id": "j1", "left_breakpoint": 100, "right_breakpoint": 200, "deleted_size": 99},
                {"sample": "s1", "read_id": "readB", "junction_id": "j1", "left_breakpoint": 110, "right_breakpoint": 210, "deleted_size": 99},
            ]
        )
        clusters = pd.DataFrame(
            {
                "junction_id": ["j1"],
                "left_breakpoint": [105],
                "right_breakpoint": [207],
                "deleted_size": [999],
                "affected_feature_label": ["MT-ND1+MT-ND2"],
                "replication_arc_context": ["major_arc_only"],
            }
        )
        corrected = apply_cluster_coordinates(reads, clusters, mt_length=1000)
        self.assertEqual(corrected["left_breakpoint"].tolist(), [105, 105])
        self.assertEqual(corrected["right_breakpoint"].tolist(), [207, 207])
        self.assertEqual(corrected["deleted_size"].tolist(), [101, 101])
        self.assertEqual(corrected["affected_feature_label"].tolist(), ["MT-ND1+MT-ND2", "MT-ND1+MT-ND2"])
        self.assertEqual(corrected["replication_arc_context"].tolist(), ["major_arc_only", "major_arc_only"])

    def test_wrapping_deletion_annotates_features_on_both_sides_of_origin(self):
        import pandas as pd

        features = pd.DataFrame(
            [
                {"gene_name": "D-loop", "feature_type": "regulatory", "start": 1, "end": 100},
                {"gene_name": "MT-CYB", "feature_type": "gene", "start": 900, "end": 1000},
            ]
        )
        impact = affected_feature_impact(features, left=950, right=50, mt_length=1000)
        self.assertEqual(impact.affected_feature_label, "D-loop+MT-CYB")
        self.assertIn("d_loop", impact.feature_impact_class)

    def test_mt_t_and_mt_r_aliases_are_rna_feature_impact(self):
        import pandas as pd

        features = pd.DataFrame(
            [
                {"gene_name": "Mt-tf", "feature_type": "tRNA", "start": 10, "end": 50},
                {"gene_name": "Mt-rnr1", "feature_type": "rRNA", "start": 70, "end": 180},
                {"gene_name": "Mt-nd1", "feature_type": "protein_coding", "start": 200, "end": 300},
            ]
        )
        impact = affected_feature_impact(features, left=1, right=190, mt_length=1000)
        self.assertEqual(impact.feature_impact_class, "rrna_trna_involved")

    def test_annotation_features_are_deduplicated_to_biological_features(self):
        import pandas as pd

        raw = pd.DataFrame(
            [
                {"gene_name": "MT-ND1", "feature_type": "exon", "start": 1, "end": 100},
                {"gene_name": "MT-ND1", "feature_type": "gene", "start": 1, "end": 100},
                {"gene_name": "MT-ND1", "feature_type": "CDS", "start": 2, "end": 99},
                {"gene_name": "MT-TF", "feature_type": "gene", "start": 120, "end": 180},
            ]
        )
        deduped = biological_features(raw)
        self.assertEqual(deduped["gene_name"].tolist(), ["MT-ND1", "MT-TF"])
        self.assertEqual(deduped.loc[deduped["gene_name"] == "MT-ND1", "feature_type"].iloc[0], "gene")

    def test_configured_wrapping_region_is_split_for_annotation(self):
        import pandas as pd

        features = pd.DataFrame(columns=["contig", "start", "end", "strand", "feature_type", "gene_id", "gene_name", "transcript_id", "product"])
        config = {"analysis": {"mt_regions": [{"name": "control_region", "start": 900, "end": 100, "reason": "wraps origin"}]}}
        with_regions = append_configured_regions(features, config, mt_length=1000)
        self.assertEqual(len(with_regions), 2)
        self.assertEqual(with_regions[["start", "end"]].astype(int).values.tolist(), [[900, 1000], [1, 100]])

    def test_replication_arc_annotations_follow_directed_deleted_interval(self):
        config = {
            "dataset": {"species": "human"},
            "references": {
                "human": {
                    "mt_length": 16569,
                    "replication_arcs": [
                        {"name": "minor_arc", "start": 408, "end": 5746},
                        {"name": "major_arc", "start": 5747, "end": 407},
                    ],
                }
            },
        }
        major_only = replication_arc_annotation(config, left=8469, right=13447, mt_length=16569)
        self.assertEqual(major_only["replication_arc_context"], "major_arc_only")
        self.assertEqual(major_only["major_arc_deleted_bp"], 4977)
        self.assertEqual(major_only["minor_arc_deleted_bp"], 0)

        reciprocal = replication_arc_annotation(config, left=13447, right=8469, mt_length=16569)
        self.assertEqual(reciprocal["replication_arc_context"], "major_and_minor_arcs")
        self.assertEqual(reciprocal["major_arc_deleted_bp"] + reciprocal["minor_arc_deleted_bp"], 11590)
        self.assertEqual(reciprocal["minor_arc_deleted_bp"], 5339)

        minor_only = replication_arc_annotation(config, left=1000, right=2000, mt_length=16569)
        self.assertEqual(minor_only["replication_arc_context"], "minor_arc_only")
        self.assertEqual(minor_only["minor_arc_deleted_bp"], 999)
        self.assertEqual(minor_only["major_arc_deleted_bp"], 0)

    def test_replication_arcs_are_not_affected_features(self):
        import pandas as pd

        features = pd.DataFrame(columns=["contig", "start", "end", "strand", "feature_type", "gene_id", "gene_name", "transcript_id", "product"])
        config = {
            "dataset": {"species": "human"},
            "references": {
                "human": {
                    "replication_arcs": [
                        {"name": "minor_arc", "start": 408, "end": 5746},
                        {"name": "major_arc", "start": 5747, "end": 407},
                    ]
                }
            },
        }
        self.assertTrue(append_configured_regions(features, config, mt_length=16569).empty)

    def test_replication_arc_configuration_must_partition_reference(self):
        config = {
            "dataset": {"species": "human"},
            "references": {
                "human": {
                    "replication_arcs": [
                        {"name": "minor_arc", "start": 408, "end": 5746},
                        {"name": "major_arc", "start": 5748, "end": 407},
                    ]
                }
            },
        }
        with self.assertRaisesRegex(ValueError, "partition"):
            replication_arc_annotation(config, left=1000, right=2000, mt_length=16569)

    def test_report_replication_arc_table_handles_wrapping_reference_arcs(self):
        config = {
            "dataset": {"species": "rat"},
            "references": {
                "rat": {
                    "mt_length": 16313,
                    "replication_arcs": [
                        {"name": "minor_arc", "display_name": "Minor arc", "start": 16160, "end": 5170},
                        {"name": "major_arc", "display_name": "Major arc", "start": 5171, "end": 16159},
                    ],
                }
            },
        }
        table = configured_replication_arc_table(config)
        self.assertEqual(table["arc_name"].tolist(), ["Minor arc", "Major arc"])
        self.assertEqual(table["wraps_coordinate_origin"].tolist(), ["yes", "no"])
        self.assertEqual(table["length_bp"].sum(), 16313)

    def test_dataset_replication_arcs_partition_each_mitochondrial_reference(self):
        repo_root = Path(__file__).resolve().parents[1]
        for config_path in sorted((repo_root / "config" / "datasets").glob("*.yaml")):
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            table = configured_replication_arc_table(config)
            self.assertEqual(set(table["arc_name"]), {"Minor arc", "Major arc"}, config_path.name)
            species = config["dataset"]["species"]
            reference = config["references"][species]
            self.assertTrue(str(reference.get("replication_arc_boundary_basis", "")).strip(), config_path.name)
            mt_length = int(reference["mt_length"])
            self.assertEqual(int(table["length_bp"].sum()), mt_length, config_path.name)

    def test_reference_support_windows_handle_origin_wrapping(self):
        self.assertEqual(circular_window(5, 10, 100), [(95, 100), (1, 15)])
        self.assertTrue(window_covered(5, 10, 100, [(90, 100), (1, 20)]))
        self.assertFalse(window_covered(5, 10, 100, [(1, 20)]))

    def test_configured_feature_aliases_replace_display_names(self):
        import pandas as pd

        features = pd.DataFrame(
            [
                {"gene_name": "AY172581.13", "gene_id": "ENSRNOG1", "feature_type": "gene", "start": 1, "end": 67},
                {"gene_name": "Mt-nd1", "gene_id": "ENSRNOG2", "feature_type": "gene", "start": 100, "end": 200},
            ]
        )
        config = {"annotations": {"feature_aliases": [{"raw_name": "AY172581.13", "display_name": "Mt-tf"}]}}
        aliased = apply_feature_aliases(features, config)
        self.assertEqual(aliased["gene_name"].tolist(), ["Mt-tf", "Mt-nd1"])
        self.assertEqual(aliased["raw_gene_name"].tolist(), ["AY172581.13", "Mt-nd1"])

    def test_exact_deletion_display_table_shows_small_call_sets(self):
        clusters = pd.DataFrame(
            [
                {"exact_deletion_id": "mtDel_1", "total_supporting_reads": 41},
                {"exact_deletion_id": "mtDel_2", "total_supporting_reads": 2},
            ]
        )
        config = {"report": {"exact_deletion_table": {"min_total_supporting_reads": 50, "max_rows": 500}}}
        display, note = exact_deletion_display_table(clusters, config)
        self.assertEqual(len(display), 2)
        self.assertIn("shows all 2 exact deletions", note)

    def test_exact_deletion_display_default_has_no_support_threshold(self):
        self.assertEqual(exact_deletion_table_settings({})["min_total_supporting_reads"], 0)
        clusters = pd.DataFrame(
            [
                {
                    "exact_deletion_id": f"mtDel_{index:04d}",
                    "total_supporting_reads": 1,
                    "known_deletion_label": "configured" if index == 0 else "",
                }
                for index in range(600)
            ]
        )
        config = {
            "report": {
                "exact_deletion_table": {
                    "min_total_supporting_reads": 0,
                    "always_include_configured_targets": True,
                    "max_rows": 500,
                }
            }
        }
        display, note = exact_deletion_display_table(clusters, config)
        self.assertEqual(len(display), 500)
        self.assertTrue((display["total_supporting_reads"] == 1).all())
        self.assertIn("supporting-read count", note)
        self.assertNotIn("at least 0", note)

    def test_rainfall_location_features_apply_configured_aliases(self):
        features = pd.DataFrame(
            [
                {"gene_name": "AY172581.9", "gene_id": "ENSRNOG1", "feature_type": "gene", "start": 68, "end": 1025},
                {"gene_name": "AY172581.3", "gene_id": "ENSRNOG2", "feature_type": "gene", "start": 1026, "end": 1093},
                {"gene_name": "Mt-nd1", "gene_id": "ENSRNOG3", "feature_type": "gene", "start": 2740, "end": 3694},
            ]
        )
        config = {
            "annotations": {
                "feature_aliases": [
                    {"raw_name": "AY172581.9", "display_name": "Mt-rnr1"},
                    {"raw_name": "AY172581.3", "display_name": "Mt-tv"},
                ]
            }
        }
        plotted = location_features(features, config)
        classes = dict(zip(plotted["name"], plotted["class"]))
        self.assertEqual(classes["Mt-rnr1"], "rRNA")
        self.assertEqual(classes["Mt-tv"], "tRNA")
        self.assertEqual(classes["Mt-nd1"], "protein_coding")

    def test_rainfall_feature_track_uses_full_mt_axis_and_clipped_labels(self):
        import matplotlib.pyplot as plt
        import pandas as pd

        features = pd.DataFrame(
            [
                {"gene_name": "MT-ND1", "feature_type": "protein_coding", "start": 3307, "end": 4262},
                {"gene_name": "MT-CYB", "feature_type": "protein_coding", "start": 14747, "end": 15887},
                {"gene_name": "MT-TP", "feature_type": "tRNA", "start": 15956, "end": 16023},
            ]
        )
        self.assertEqual(mitochondrial_axis_bounds(features), (1.0, 16023.0))

        fig, ax = plt.subplots()
        try:
            x_min, x_max = draw_feature_track_axis(ax, features)
            self.assertEqual((x_min, x_max), (1.0, 16023.0))
            self.assertEqual(ax.get_xlim(), (1.0, 16023.0))
            feature_labels = [text for text in ax.texts if text.get_text() in {"MT-ND1", "MT-CYB"}]
            self.assertEqual(len(feature_labels), 2)
            self.assertTrue(all(text.get_clip_on() for text in feature_labels))
        finally:
            plt.close(fig)

    def test_mitochondrial_feature_palette_is_shared_by_linear_and_circular_tracks(self):
        from matplotlib.colors import to_hex

        self.assertEqual(
            MITOCHONDRIAL_FEATURE_COLORS,
            {
                "protein_coding": "#7CAE00",
                "rRNA": "#00BFC4",
                "tRNA": "#C77CFF",
                "region": "#F8766D",
                "other": "#9AA3AF",
            },
        )
        features = pd.DataFrame(
            [
                {"name": "MT-ND1", "class": "protein_coding", "start": 100, "end": 300},
                {"name": "MT-RNR1", "class": "rRNA", "start": 350, "end": 550},
                {"name": "MT-TF", "class": "tRNA", "start": 600, "end": 680},
                {"name": "D-loop/control", "class": "region", "start": 700, "end": 900},
            ]
        )
        expected = {value.lower() for key, value in MITOCHONDRIAL_FEATURE_COLORS.items() if key != "other"}

        linear_fig, linear_ax = plt.subplots()
        circular_fig, circular_ax = plt.subplots()
        try:
            draw_location_feature_track(linear_ax, features, 1000)
            add_feature_ring(circular_ax, features, 1000)
            linear_colors = {to_hex(patch.get_facecolor()).lower() for patch in linear_ax.patches}
            circular_colors = {
                to_hex(patch.get_facecolor()).lower()
                for patch in circular_ax.patches
                if patch.get_fill()
            }
            self.assertTrue(expected.issubset(linear_colors))
            self.assertTrue(expected.issubset(circular_colors))
        finally:
            plt.close(linear_fig)
            plt.close(circular_fig)

    def test_circular_chord_endpoints_use_clockwise_mt_coordinates(self):
        genome_length = 1000
        self.assertAlmostEqual(circle_point(1, 1.0, genome_length)[0], 0.0, places=10)
        self.assertAlmostEqual(circle_point(1, 1.0, genome_length)[1], 1.0, places=10)
        quarter = circle_point(251, 1.0, genome_length)
        self.assertAlmostEqual(quarter[0], 1.0, places=10)
        self.assertAlmostEqual(quarter[1], 0.0, places=10)
        path = chord_path(100, 800, 0.9, genome_length)
        self.assertTrue((path.vertices[0] == circle_point(100, 0.9, genome_length)).all())
        self.assertTrue((path.vertices[-1] == circle_point(800, 0.9, genome_length)).all())

    def test_circular_chord_metadata_includes_affected_features(self):
        import tempfile

        calls = pd.DataFrame(
            [
                {
                    "_plot_support": 0.25,
                    "supporting_reads": 3,
                    "_support_rank": 1,
                    "exact_deletion_id": "mtDel_00100_00500_00399",
                    "left_breakpoint": 100,
                    "right_breakpoint": 500,
                    "deleted_size": 399,
                    "affected_feature_label": "MT-ND1+MT-CO1",
                    "replication_arc_context": "major_and_minor_arcs",
                    "major_arc_deleted_bp": 200,
                    "minor_arc_deleted_bp": 199,
                }
            ]
        )
        node_id = chord_dom_id("all", "treated", calls.iloc[0])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chords.svg"
            path.write_text(
                f'<svg xmlns="http://www.w3.org/2000/svg"><g id="{node_id}"><path d="M 0 0"/></g></svg>',
                encoding="utf-8",
            )
            add_svg_chord_metadata(path, calls, "all", "treated")
            svg = path.read_text(encoding="utf-8")
        self.assertIn('data-affected-features="MT-ND1+MT-CO1"', svg)
        self.assertIn('data-support-label="Deletion support"', svg)

    def test_circular_comparison_plot_accepts_a_zero_row_table(self):
        calls = prepare_comparison_calls(pd.DataFrame(), {}, genome_length=1000)
        self.assertTrue(calls.empty)
        self.assertTrue(
            {
                "exact_deletion_id",
                "left_group",
                "right_group",
                "left_breakpoint",
                "right_breakpoint",
                "deleted_size",
                "difference_per_million_mt_reads",
                "left_total_supporting_reads",
                "right_total_supporting_reads",
                "_absolute_difference",
                "_total_supporting_observations",
            }.issubset(calls.columns)
        )

    def test_circular_comparison_plot_excludes_group_pairs_without_evidence(self):
        comparison = pd.DataFrame(
            [
                {
                    "exact_deletion_id": "mtDel_zero",
                    "left_group": "control",
                    "right_group": "treated",
                    "left_breakpoint": 100,
                    "right_breakpoint": 300,
                    "deleted_size": 199,
                    "difference_per_million_mt_reads": 0.0,
                    "left_total_supporting_reads": 0,
                    "right_total_supporting_reads": 0,
                },
                {
                    "exact_deletion_id": "mtDel_supported",
                    "left_group": "control",
                    "right_group": "treated",
                    "left_breakpoint": 400,
                    "right_breakpoint": 700,
                    "deleted_size": 299,
                    "difference_per_million_mt_reads": 0.5,
                    "left_total_supporting_reads": 0,
                    "right_total_supporting_reads": 2,
                },
            ]
        )
        calls = prepare_comparison_calls(comparison, {}, genome_length=1000)
        self.assertEqual(calls["exact_deletion_id"].tolist(), ["mtDel_supported"])
        self.assertEqual(calls["_total_supporting_observations"].tolist(), [2])

    def test_circular_comparison_plot_replaces_existing_replication_arc_columns(self):
        comparison = pd.DataFrame(
            [
                {
                    "exact_deletion_id": "mtDel_supported",
                    "left_group": "control",
                    "right_group": "treated",
                    "left_breakpoint": 100,
                    "right_breakpoint": 300,
                    "deleted_size": 199,
                    "difference_per_million_mt_reads": 0.5,
                    "left_total_supporting_reads": 1,
                    "right_total_supporting_reads": 2,
                    "replication_arc_context": "stale",
                    "minor_arc_deleted_bp": -1,
                    "major_arc_deleted_bp": -1,
                }
            ]
        )
        calls = prepare_comparison_calls(comparison, {}, genome_length=1000)
        for column in ["replication_arc_context", "minor_arc_deleted_bp", "major_arc_deleted_bp"]:
            self.assertEqual(calls.columns.tolist().count(column), 1)
        self.assertEqual(calls.loc[0, "replication_arc_context"], "not_configured")

    def test_inline_svg_ids_and_references_are_namespaced(self):
        svg = (
            '<svg aria-labelledby="title description">'
            '<style>#glyph{fill:#ffffff} g[id^="static-points-"]{display:none;}</style>'
            '<defs><path id="glyph"/><clipPath id="clip"><path/></clipPath></defs>'
            '<title id="title">Title</title><desc id="description">Description</desc>'
            '<g id="static-points-1" clip-path="url(#clip)">'
            '<use href="#glyph"/><use xlink:href="#glyph"/>'
            '</g><g id="glyph" class="duplicate-source-id"/></svg>'
        )
        namespaced = namespace_inline_svg(svg, "comparison/control")
        prefix = "report_svg__comparison_control__"

        for node_id in ["glyph", "clip", "title", "description", "static-points-1"]:
            self.assertIn(f'id="{prefix}{node_id}"', namespaced)
            self.assertNotIn(f'id="{node_id}"', namespaced)
        self.assertIn(f'url(#{prefix}clip)', namespaced)
        self.assertEqual(namespaced.count(f'href="#{prefix}glyph"'), 2)
        self.assertIn(f'id="{prefix}glyph__duplicate_2"', namespaced)
        self.assertIn(f'#{prefix}glyph{{fill:#ffffff}}', namespaced)
        self.assertIn(f'[id^="{prefix}static-points-"]', namespaced)
        self.assertIn(f'aria-labelledby="{prefix}title {prefix}description"', namespaced)

    def test_report_circular_chord_panels_include_interactive_controls(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            plots = Path(tmp)
            location = plots / "circular_breakpoint_chords_all.pdf"
            location.write_bytes(b"pdf")
            (plots / "circular_breakpoint_chords_all__group_a.pdf").write_bytes(b"pdf")
            (plots / "circular_breakpoint_chords_all__group_a__interactive.svg").write_text(
                '<svg data-group="group A"><g class="deletion-chord" data-support="1" '
                'data-observations="2" data-baseline="1"/></svg>',
                encoding="utf-8",
            )
            location_html = circular_location_plot_panel(str(location), "Location", "Caption", "plots")
            self.assertIn("data-support-slider", location_html)
            self.assertIn("data-observation-filter", location_html)
            self.assertIn('data-size-filter', location_html)
            self.assertIn('value="1000">&ge; 1,000 bp', location_html)
            self.assertIn('value="10000">&ge; 10,000 bp', location_html)
            self.assertIn('value="100">&ge; 100', location_html)
            self.assertIn('value="200">&ge; 200', location_html)
            self.assertIn("Choosing a numeric observation cutoff moves the support slider", location_html)
            self.assertIn("group A", location_html)

            comparison = plots / "exact_deletion_comparison_chords.pdf"
            comparison.write_bytes(b"pdf")
            (plots / "exact_deletion_comparison_chords__control_vs_treated.pdf").write_bytes(b"pdf")
            (plots / "exact_deletion_comparison_chords__control_vs_treated.svg").write_text(
                '<svg data-left-group="control" data-right-group="treated">'
                '<g class="comparison-chord" data-total-observations="3" data-absolute-difference="1"/>'
                '</svg>',
                encoding="utf-8",
            )
            comparison_html = circular_comparison_plot_panel(
                str(comparison), "Comparisons", "Caption", "plots"
            )
            self.assertIn("Replicate-significant", comparison_html)
            self.assertIn("Read-depth enriched", comparison_html)
            self.assertIn("treated compared with control", comparison_html)
            self.assertNotIn("major arc", comparison_html.lower())

    def test_rainfall_interactive_sidecar_has_all_point_metadata_and_controls(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "deletion_rainfall_left_breakpoint.pdf"
            calls = pd.DataFrame(
                [
                    {
                        "_plot_group": "treated",
                        "left_breakpoint": 100,
                        "right_breakpoint": 500,
                        "deleted_size": 399,
                        "_plot_x": 100,
                        "_plot_support": 0.2,
                        "supporting_reads": 3,
                        "exact_deletion_id": "mtDel_00100_00500_00399",
                        "crosses_origin": False,
                        "affected_feature_label": "MT-CO1+MT-CO2",
                        "replication_arc_context": "major_arc_only",
                        "major_arc_deleted_bp": 399,
                        "minor_arc_deleted_bp": 0,
                        "known_deletion_label": "",
                    },
                    {
                        "_plot_group": "treated",
                        "left_breakpoint": 16000,
                        "right_breakpoint": 200,
                        "deleted_size": 768,
                        "_plot_x": 16000,
                        "_plot_support": 4.0,
                        "supporting_reads": 8,
                        "exact_deletion_id": "mtDel_16000_00200_00768",
                        "crosses_origin": True,
                        "affected_feature_label": "MT-ND1",
                        "replication_arc_context": "major_and_minor_arcs",
                        "major_arc_deleted_bp": 500,
                        "minor_arc_deleted_bp": 268,
                        "known_deletion_label": "configured target",
                    },
                ]
            )
            location_rainfall(
                calls,
                ["treated"],
                pd.DataFrame(),
                {},
                16569,
                str(path),
                "Rainfall",
                "left_breakpoint",
                "Left breakpoint",
                "Normalized support",
            )
            interactive = path.with_name("deletion_rainfall_left_breakpoint__treated__interactive.svg")
            svg = interactive.read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="rainfall-point"'), 2)
            self.assertIn('data-exact-deletion-id="mtDel_00100_00500_00399"', svg)
            self.assertIn('data-crosses-origin="yes"', svg)
            self.assertIn('data-arc-context="major_and_minor_arcs"', svg)
            self.assertIn('data-call-count="2"', svg)
            self.assertIn('data-support-label="Normalized support"', svg)
            self.assertNotIn("support rank within this group", svg)

            panel = rainfall_location_plot_panel(str(path), "Rainfall", "Caption", "plots")
            self.assertIn("data-rainfall-controls", panel)
            self.assertIn("data-rainfall-support-slider", panel)
            self.assertIn("data-observation-filter", panel)
            self.assertIn('data-size-filter', panel)
            self.assertIn('value="200">&ge; 200', panel)
            self.assertIn('value="1000">&ge; 1,000 bp', panel)
            self.assertIn('value="10000">&ge; 10,000 bp', panel)
            self.assertIn("eligible call set loaded for this view contains 2 calls", panel)

    def test_rainfall_support_scaling_uses_observed_range(self):
        support = pd.Series([0.0034, 0.01, 0.1, 0.348])
        support_min, support_max = rainfall_support_limits(support)
        self.assertAlmostEqual(support_min, 0.0034)
        self.assertAlmostEqual(support_max, 0.348)

        sizes = rainfall_point_sizes(support, support_min, support_max)
        self.assertLess(sizes[0], 3)
        self.assertGreater(sizes[-1], 500)
        self.assertTrue(all(a < b for a, b in zip(sizes, sizes[1:])))

        scale_min, scale_max = support_scale_limits(support_min, support_max)
        self.assertEqual((scale_min, scale_max), (0.002, 0.5))
        ticks = support_legend_values(scale_min, scale_max)
        self.assertEqual(ticks, [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5])

        wide_ticks = support_legend_values(*support_scale_limits(0.16, 821))
        self.assertEqual(wide_ticks, [0.1, 1.0, 10.0, 100.0, 1000.0])
        wide_size_ticks = support_size_legend_values(*support_scale_limits(0.16, 821))
        self.assertEqual(wide_size_ticks[0], 0.1)
        self.assertEqual(wide_size_ticks[-1], 1000.0)
        self.assertLess(len(wide_size_ticks), 7)

    def test_endpoint_density_interactive_sidecar_has_bin_hover_metadata(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pooled_breakpoint_support_density.pdf"
            calls = pd.DataFrame(
                [
                    {"_plot_group": "treated", "left_breakpoint": 100, "right_breakpoint": 250, "_plot_support": 2.0, "supporting_reads": 4},
                    {"_plot_group": "treated", "left_breakpoint": 150, "right_breakpoint": 650, "_plot_support": 3.0, "supporting_reads": 5},
                ]
            )
            endpoint_density_pages(
                calls,
                ["treated"],
                pd.DataFrame(),
                {},
                1000,
                str(path),
                "Pooled Breakpoint Support Density",
                "Normalized support",
                bin_size=100,
                smooth_bins=3,
            )
            interactive = path.with_name("pooled_breakpoint_support_density__treated__interactive.svg")
            svg = interactive.read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="endpoint-density-bin"'), 10)
            self.assertIn('data-plot-type="endpoint-density"', svg)
            self.assertIn('data-bin-start="101.0"', svg)
            self.assertIn('data-left-support="2.0"', svg)
            self.assertIn('data-left-raw-supporting-reads="4.0"', svg)
            self.assertIn('data-smoothed-support=', svg)

            panel = endpoint_density_plot_panel(str(path), "Density", "Caption", "plots")
            self.assertIn("Hover a density bin", panel)
            self.assertIn("10 bins", panel)

    def test_endpoint_density_separates_exact_call_count_from_raw_support(self):
        calls = pd.DataFrame(
            [
                {"_plot_group": "treated", "left_breakpoint": 50, "right_breakpoint": 250, "_plot_support": 0.4, "supporting_reads": 110},
                {"_plot_group": "treated", "left_breakpoint": 60, "right_breakpoint": 650, "_plot_support": 0.02, "supporting_reads": 2},
            ]
        )
        density = pooled_endpoint_density(calls, genome_length=1000, bin_size=100, smooth_bins=1)
        first = density.loc[density["bin_start"] == 1].iloc[0]
        self.assertEqual(first["left_endpoint_count"], 2)
        self.assertEqual(first["left_raw_supporting_reads"], 112)
        self.assertAlmostEqual(first["left_support"], 0.42)
        self.assertEqual(first["endpoint_count"], 2)
        self.assertEqual(first["raw_supporting_reads"], 112)

    def test_breakpoint_pair_map_has_point_hover_metadata(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "breakpoint_pair_support_map.pdf"
            calls = pd.DataFrame(
                [
                    {
                        "_plot_group": "treated",
                        "left_breakpoint": 900,
                        "right_breakpoint": 100,
                        "deleted_size": 200,
                        "_plot_support": 2.5,
                        "supporting_reads": 4,
                        "_support_rank": 1,
                        "exact_deletion_id": "mtDel_00900_00100_00200",
                        "affected_feature_label": "MT-ND1+MT-CO1",
                        "replication_arc_context": "major_and_minor_arcs",
                        "major_arc_deleted_bp": 120,
                        "minor_arc_deleted_bp": 80,
                    }
                ]
            )
            breakpoint_pair_support_map(
                calls,
                ["treated"],
                pd.DataFrame(),
                {},
                1000,
                str(path),
                "Normalized support",
            )
            interactive = path.with_name("breakpoint_pair_support_map__treated__interactive.svg")
            svg = interactive.read_text(encoding="utf-8")
            self.assertIn('class="breakpoint-pair-point"', svg)
            self.assertIn('g[id^="breakpoint-pair-static-points-"]{display:none;}', svg)
            self.assertIn('id="breakpoint-pair-static-points-', svg)
            self.assertIn('class="breakpoint-pair-visible-point"', svg)
            self.assertIn('class="breakpoint-pair-hit-target"', svg)
            self.assertIn('fill-opacity="0.86"', svg)
            self.assertIn('fill-opacity="0"', svg)
            self.assertIn('pointer-events="all"', svg)
            self.assertNotIn('id="breakpoint-pair-rank-1"', svg)
            self.assertIn('data-exact-deletion-id="mtDel_00900_00100_00200"', svg)
            self.assertIn('data-affected-features="MT-ND1+MT-CO1"', svg)
            self.assertIn('data-rank="1"', svg)
            self.assertIn('data-plot-type="breakpoint-pair-map"', svg)
            self.assertIn('data-support-label="Normalized support"', svg)

            panel = breakpoint_pair_plot_panel(str(path), "Breakpoint pairs", "Caption", "plots")
            self.assertIn("Hover a point to inspect its breakpoint pair", panel)
            self.assertIn("data-breakpoint-pair-controls", panel)
            self.assertIn('value="1000">&ge; 1,000 bp', panel)
            self.assertIn("Open PDF", panel)

    def test_ordination_report_panel_mentions_sample_hover(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "exact_deletion_pca.pdf"
            path.write_bytes(b"pdf")
            path.with_suffix(".svg").write_text(
                '<svg data-plot-type="ordination" data-point-count="2">'
                '<circle class="ordination-point" data-sample="S1"/></svg>',
                encoding="utf-8",
            )
            panel = plot_panel(str(path), "PCA", "Caption", "plots")
            self.assertIn("Hover a sample to inspect its coordinates", panel)

    def test_sample_point_report_panel_mentions_sample_hover(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "burden.pdf"
            path.write_bytes(b"pdf")
            path.with_suffix(".svg").write_text(
                '<svg data-plot-type="sample-points" data-point-count="1">'
                '<circle class="sample-point" data-sample="S1"/></svg>',
                encoding="utf-8",
            )
            panel = plot_panel(str(path), "Burden", "Caption", "plots")
            self.assertIn("Hover a sample point or group-mean marker", panel)

    def test_bar_report_panel_mentions_bar_hover(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "features.pdf"
            path.write_bytes(b"pdf")
            path.with_suffix(".svg").write_text(
                '<svg data-plot-type="bar-chart" data-bar-count="1">'
                '<rect class="bar-plot-bar" data-category="MT-ND1" data-value="2"/></svg>',
                encoding="utf-8",
            )
            panel = plot_panel(str(path), "Features", "Caption", "plots")
            self.assertIn("Hover a bar to inspect its category", panel)

    def test_rainfall_size_filter_has_a_change_handler(self):
        js = (Path(__file__).resolve().parents[1] / "report_assets" / "circular_chords.js").read_text(encoding="utf-8")
        self.assertIn("syncRainfallSupportToObservations();", js)
        self.assertIn("syncChordSupportToObservations();", js)
        self.assertIn("syncBreakpointSupportToObservations();", js)
        self.assertIn("supportSliderValueForThreshold", js)
        self.assertIn("target.classList.contains('bar-plot-bar')", js)

    def test_bar_plots_have_hover_metadata(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3", "S4"],
                    "group": ["control", "control", "deletion", "deletion"],
                }
            )
            reads = pd.DataFrame(
                {
                    "sample": ["S1", "S1", "S2", "S3", "S4"],
                    "deleted_size": [100, 200, 5100, 5200, 1000],
                }
            )
            size_distribution(reads, samples, "group", str(root / "size.pdf"), "Size", weighted=True)
            size_svg = (root / "size.svg").read_text(encoding="utf-8")
            self.assertIn('data-plot-type="bar-chart"', size_svg)
            self.assertIn('class="bar-plot-bar"', size_svg)
            self.assertIn('data-bin-start=', size_svg)
            self.assertIn('data-group-values=', size_svg)
            self.assertIn('supportingReads', size_svg)
            self.assertIn('control', size_svg)
            self.assertIn('deletion', size_svg)

            matrix = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3", "S4"],
                    "MT-ND1": [1, 2, 3, 4],
                    "MT-CO1": [2, 1, 4, 3],
                }
            )
            category_bar(matrix, samples, "group", str(root / "features.pdf"), "Features", "Support")
            feature_svg = (root / "features.svg").read_text(encoding="utf-8")
            self.assertIn('data-category="MT-ND1"', feature_svg)
            self.assertIn('data-group="control"', feature_svg)
            self.assertIn('data-value-label="Support"', feature_svg)

            per_gene = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3", "S4"],
                    "group": ["control", "control", "deletion", "deletion"],
                    "feature": ["MT-ND1"] * 4,
                    "support_per_million_mt_reads": [1, 2, 3, 4],
                }
            )
            per_gene_plot(per_gene, pd.DataFrame(), "group", str(root / "genes.pdf"), "Support")
            gene_svg = (root / "genes.svg").read_text(encoding="utf-8")
            self.assertIn('data-feature="MT-ND1"', gene_svg)
            self.assertIn('class="bar-plot-bar"', gene_svg)

            recurrence = pd.DataFrame(
                {
                    "exact_deletion_id": ["d1", "d2"],
                    "left_breakpoint": [100, 200],
                    "right_breakpoint": [300, 400],
                    "deleted_size": [200, 200],
                    "total_supporting_reads": [5, 4],
                    "affected_feature_label": ["MT-ND1", "MT-CO1"],
                }
            )
            exact_matrix = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3", "S4"],
                    "d1": [1, 2, 3, 4],
                    "d2": [2, 1, 4, 3],
                }
            )
            exact_recurrence(recurrence, exact_matrix, samples, "group", str(root / "recurrence.pdf"), "Support")
            recurrence_svg = (root / "recurrence.svg").read_text(encoding="utf-8")
            self.assertIn('data-deletion-id="d1"', recurrence_svg)
            self.assertIn('data-value-label="Support"', recurrence_svg)

    def test_category_bar_without_group_uses_default_seaborn_color(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples = pd.DataFrame({"sample": ["S1"]})
            matrix = pd.DataFrame(
                {
                    "sample": ["S1"],
                    "MT-ND1": [2],
                    "MT-CO1": [1],
                }
            )

            category_bar(matrix, samples, "", str(root / "features.pdf"), "Features", "Support")

            feature_svg = (root / "features.svg").read_text(encoding="utf-8")
            self.assertIn('data-category="MT-ND1"', feature_svg)
            self.assertIn('data-category="MT-CO1"', feature_svg)
            self.assertNotIn('data-group=', feature_svg)

    def test_ordination_has_sample_hover_metadata(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "exact_deletion_pca.pdf"
            matrix = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3"],
                    "mtDel_a": [1.0, 0.0, 2.0],
                    "mtDel_b": [0.0, 2.0, 1.0],
                }
            )
            samples = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3"],
                    "condition": ["control", "treated", "treated"],
                    "biological_replicate": ["1", "1", "2"],
                    "layout": ["single", "single", "single"],
                    "tissue": ["muscle", "muscle", "muscle"],
                }
            )
            ordination(matrix, samples, "condition", str(path), "Exact Deletion PCA", "pca")
            svg = path.with_suffix(".svg").read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="ordination-point"'), 3)
            self.assertIn('data-plot-type="ordination"', svg)
            self.assertIn('data-sample="S1"', svg)
            self.assertIn('data-group="control"', svg)

    def test_sample_point_plots_have_sample_hover_metadata(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            burden_path = Path(tmp) / "burden.pdf"
            burden = pd.DataFrame(
                {
                    "sample": ["S1", "S2", "S3", "S4"],
                    "group": ["control", "control", "deletion", "deletion"],
                    "deletion_support_per_million_mt_reads": [1.0, 2.0, 3.0, 4.0],
                    "biological_replicate": ["R1", "R2", "R1", "R2"],
                    "layout": ["single"] * 4,
                    "tissue": ["muscle"] * 4,
                }
            )
            burden_plot(
                burden,
                "group",
                str(burden_path),
                "deletion_support_per_million_mt_reads",
                "Burden",
                "Normalized support",
            )
            burden_svg = burden_path.with_suffix(".svg").read_text(encoding="utf-8")
            self.assertEqual(burden_svg.count('class="sample-point"'), 4)
            self.assertEqual(burden_svg.count('class="group-mean-point"'), 2)
            self.assertIn('data-sample-count="2"', burden_svg)
            self.assertIn('data-sample="S1"', burden_svg)
            self.assertIn('data-biological-replicate="R1"', burden_svg)

            factorial_path = Path(tmp) / "factorial.pdf"
            factorial = burden.assign(
                age=["young", "young", "old", "old"],
                treatment=["control", "deletion", "control", "deletion"],
            )
            factorial_interaction_plot(
                factorial,
                str(factorial_path),
                "deletion_support_per_million_mt_reads",
                "Factorial",
                "Normalized support",
            )
            factorial_svg = factorial_path.with_suffix(".svg").read_text(encoding="utf-8")
            self.assertEqual(factorial_svg.count('class="sample-point"'), 4)
            self.assertEqual(factorial_svg.count('class="group-mean-point"'), 4)
            self.assertIn('data-age="young"', factorial_svg)
            self.assertIn('data-treatment="deletion"', factorial_svg)
    def test_pooled_endpoint_density_tracks_left_and_right_support(self):
        calls = pd.DataFrame(
            [
                {"left_breakpoint": 10, "right_breakpoint": 120, "_plot_support": 2.0},
                {"left_breakpoint": 15, "right_breakpoint": 180, "_plot_support": 3.0},
                {"left_breakpoint": 220, "right_breakpoint": 20, "_plot_support": 5.0},
            ]
        )
        density = pooled_endpoint_density(calls, genome_length=300, bin_size=100, smooth_bins=1)
        first = density.loc[density["bin_index"] == 0].iloc[0]
        second = density.loc[density["bin_index"] == 1].iloc[0]
        third = density.loc[density["bin_index"] == 2].iloc[0]
        self.assertEqual(first["left_support"], 5.0)
        self.assertEqual(first["right_support"], 5.0)
        self.assertEqual(second["right_support"], 5.0)
        self.assertEqual(third["left_support"], 5.0)
        self.assertEqual(first["summed_support"], 10.0)

    def test_endpoint_density_low_support_uses_observed_y_range(self):
        calls = pd.DataFrame(
            [{"left_breakpoint": 100, "right_breakpoint": 800, "_plot_support": 0.04}]
        )
        fig = endpoint_density_figure(
            calls,
            pd.DataFrame(),
            {},
            1000,
            "density",
            "support per million usable reads",
            bin_size=50,
            smooth_bins=3,
        )
        try:
            self.assertGreater(fig.axes[0].get_ylim()[1], 0.04)
            self.assertLess(fig.axes[0].get_ylim()[1], 0.1)
        finally:
            plt.close(fig)

    def test_endpoint_density_hotspot_spacing_is_configurable(self):
        density = pd.DataFrame(
            {
                "bin_midpoint": [100, 200, 300, 400, 500],
                "smoothed_summed_support": [0, 5, 0, 4, 0],
            }
        )
        hotspots = endpoint_density_hotspots(density, genome_length=1000, min_spacing_bp=250)
        self.assertEqual([row["coord"] for row in hotspots], [200.0])

    def test_rainfall_y_axis_min_tracks_deletion_cutoff(self):
        self.assertEqual(rainfall_y_axis_min(pd.Series([104, 500, 5000])), 100)
        self.assertEqual(rainfall_y_axis_min(pd.Series([12, 500])), 10)

    def test_report_table_html_wraps_tables_for_horizontal_scroll(self):
        import pandas as pd

        html = table_html(pd.DataFrame([{"wide_column_name": "value"}]))
        self.assertIn('class="table-wrap"', html)
        self.assertIn('class="dataframe data-table"', html)
        self.assertEqual(html.count('class="table-wrap"'), 1)

    def test_report_assay_limitations_are_configuration_driven(self):
        nanopore_unknown = assay_limitations({"dataset": {"read_technology": "nanopore", "molecule_type": "unknown"}})
        self.assertIn("coordinate-focused deletion-like evidence", nanopore_unknown)
        self.assertIn("Molecule type is not specified", nanopore_unknown)

        illumina_rna = assay_limitations({"dataset": {"read_technology": "illumina", "molecule_type": "rna"}})
        self.assertIn("RNA read support does not directly measure", illumina_rna)

        nanopore_single_cell_rna = assay_limitations(
            {
                "dataset": {
                    "read_technology": "nanopore",
                    "molecule_type": "rna",
                    "assay_type": "single_cell_rna_seq",
                }
            }
        )
        self.assertIn("RNA read support does not directly measure", nanopore_single_cell_rna)
        self.assertIn("Single-cell RNA-seq support", nanopore_single_cell_rna)
        self.assertNotIn("Molecule type is not specified", nanopore_single_cell_rna)

    def test_report_alternative_explanations_are_configuration_driven(self):
        illumina_rna = potential_alternative_explanations(
            {"dataset": {"read_technology": "illumina", "molecule_type": "rna", "assay_type": "bulk_rna_seq"}}
        )
        illumina_applies = set(illumina_rna["Applies to"])
        self.assertIn("All datasets", illumina_applies)
        self.assertIn("Circular remapping", illumina_applies)
        self.assertIn("Illumina", illumina_applies)
        self.assertIn("RNA", illumina_applies)
        self.assertNotIn("Nanopore", illumina_applies)
        self.assertNotIn("DNA", illumina_applies)
        self.assertNotIn("Single-cell RNA-seq", illumina_applies)

        nanopore_single_cell_rna = potential_alternative_explanations(
            {
                "dataset": {
                    "read_technology": "nanopore",
                    "molecule_type": "rna",
                    "assay_type": "single_cell_rna_seq",
                }
            }
        )
        nanopore_applies = set(nanopore_single_cell_rna["Applies to"])
        self.assertIn("Nanopore", nanopore_applies)
        self.assertIn("RNA", nanopore_applies)
        self.assertIn("Single-cell RNA-seq", nanopore_applies)
        self.assertNotIn("Illumina", nanopore_applies)
        self.assertNotIn("DNA", nanopore_applies)

        nanopore_dna = potential_alternative_explanations(
            {"dataset": {"read_technology": "nanopore", "molecule_type": "dna", "assay_type": "genomic_dna"}}
        )
        dna_applies = set(nanopore_dna["Applies to"])
        self.assertIn("Nanopore", dna_applies)
        self.assertIn("DNA", dna_applies)
        self.assertNotIn("RNA", dna_applies)

        unknown = potential_alternative_explanations({"dataset": {}})
        unknown_applies = set(unknown["Applies to"])
        self.assertIn("Unknown read technology", unknown_applies)
        self.assertIn("Unknown molecule type", unknown_applies)
        self.assertNotIn("Illumina", unknown_applies)
        self.assertNotIn("Nanopore", unknown_applies)
        self.assertNotIn("RNA", unknown_applies)
        self.assertNotIn("DNA", unknown_applies)

    def test_report_assumptions_render_alternatives_before_arc_explanation(self):
        report_html = assumptions_section(
            {"dataset": {"read_technology": "illumina", "molecule_type": "rna", "assay_type": "bulk_rna_seq"}},
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )
        alternatives_heading = "Potential Alternative Explanations For Deletion-like Evidence"
        arc_heading = "How The Deleted Arc Is Assigned"
        self.assertIn(alternatives_heading, report_html)
        self.assertIn("Mitochondrial transcript processing", report_html)
        self.assertIn("Short or non-unique split anchors", report_html)
        self.assertNotIn("Basecalling errors and difficult sequence contexts", report_html)
        self.assertLess(report_html.index(alternatives_heading), report_html.index(arc_heading))

    def test_dataset_input_validation_rejects_layout_and_coordinate_mismatches(self):
        config = {
            "dataset": {
                "name": "example",
                "species": "human",
                "library_strategy": "single_end_short_read",
                "group_columns": ["condition"],
            },
            "references": {"human": {"mt_length": 16569}},
            "analysis": {
                "known_deletions": [
                    {
                        "name": "consistent",
                        "left_breakpoint": 8469,
                        "right_breakpoint": 13447,
                        "deleted_size": 4977,
                    }
                ]
            },
        }
        sample = {
            "sample": "s1",
            "dataset": "example",
            "species": "human",
            "layout": "single",
            "fastq_1": "reads.fastq.gz",
            "fastq_2": "",
            "condition": "case",
        }
        validate_dataset_inputs(config, [sample])

        paired_config = {**config, "dataset": {**config["dataset"], "library_strategy": "paired_end_short_read"}}
        with self.assertRaisesRegex(SystemExit, "single-end layout"):
            validate_dataset_inputs(paired_config, [sample])

        bad_target_config = {
            **config,
            "analysis": {
                "known_deletions": [
                    {
                        "name": "inconsistent",
                        "left_breakpoint": 8470,
                        "right_breakpoint": 13447,
                        "deleted_size": 4977,
                    }
                ]
            },
        }
        with self.assertRaisesRegex(SystemExit, "imply 4976 bases"):
            validate_dataset_inputs(bad_target_config, [sample])

    def test_report_describes_boolean_trimming_behavior(self):
        enabled = method_section({"dataset": {}, "qc": {"trim_reads": True}}, pd.DataFrame())
        disabled = method_section({"dataset": {}, "qc": {"trim_reads": False}}, pd.DataFrame())
        self.assertIn("fastp enabled by dataset configuration", enabled)
        self.assertIn("disabled by dataset configuration", disabled)
        self.assertNotIn("adapter rate", enabled)

    def test_known_sequence_summary_links_matching_reads_to_read_names(self):
        import pandas as pd
        import tempfile

        summary = pd.DataFrame(
            [{"sample": "s1", "deletion_id": "del1", "matching_reads": 2, "reads_examined": 10}]
        )
        hits = pd.DataFrame(
            [
                {"sample": "s1", "deletion_id": "del1", "read_id": "readA", "mate": "R1"},
                {"sample": "s1", "deletion_id": "del1", "read_id": "readB", "mate": "R1"},
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            html_cells = write_configured_sequence_read_lists(summary, hits, Path(tmp))
            table = table_html(summary, html_cells=html_cells)
            read_list = Path(tmp) / "configured_sequence__s1__del1.read_names.tsv"
            self.assertTrue(read_list.exists())
            content = read_list.read_text(encoding="utf-8")
        self.assertIn('href="read_lists/configured_sequence__s1__del1.read_names.tsv"', table)
        self.assertIn(">2</a>", table)
        self.assertIn("readA", content)
        self.assertIn("readB", content)

    def test_exact_deletion_read_lists_and_support_links(self):
        import pandas as pd
        import tempfile

        reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "exact_deletion_id": "mtDel_1", "junction_id": "mtDel_1", "left_breakpoint": 1, "right_breakpoint": 10},
                {"sample": "s1", "read_id": "readB", "exact_deletion_id": "mtDel_1", "junction_id": "mtDel_1", "left_breakpoint": 1, "right_breakpoint": 10},
                {"sample": "s2", "read_id": "readC", "exact_deletion_id": "mtDel_2", "junction_id": "mtDel_2", "left_breakpoint": 20, "right_breakpoint": 30},
            ]
        )
        clusters = pd.DataFrame(
            [
                {"exact_deletion_id": "mtDel_1", "total_supporting_reads": 2},
                {"exact_deletion_id": "mtDel_2", "total_supporting_reads": 1},
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            manifest = write_exact_deletion_read_lists(reads, Path(tmp))
            self.assertTrue((Path(tmp) / "mtDel_1.read_names.tsv").exists())
            self.assertTrue((Path(tmp) / "manifest.tsv").exists())
            content = (Path(tmp) / "mtDel_1.read_names.tsv").read_text(encoding="utf-8")
            self.assertIn("readA", content)
            links = exact_deletion_support_read_links(clusters, manifest)
        self.assertIn((0, "total_supporting_reads"), links)
        self.assertIn('href="read_lists/mtDel_1.read_names.tsv"', links[(0, "total_supporting_reads")])

    def test_sequence_remap_overlap_table_links_read_sets(self):
        import pandas as pd
        import tempfile

        config = {
            "analysis": {
                "known_sequence_searches": [
                    {"id": "mtDNA_100_500", "name": "configured deletion 100-500"}
                ]
            }
        }
        known_hits = pd.DataFrame(
            [
                {"sample": "s1", "deletion_id": "mtDNA_100_500", "read_id": "readA/1 cell=1"},
                {"sample": "s1", "deletion_id": "mtDNA_100_500", "read_id": "readB cell=2"},
            ]
        )
        junction_reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "exact_deletion_id": "mtDel_00100_00500_00399", "junction_id": "j1", "left_breakpoint": 100, "right_breakpoint": 500, "deleted_size": 399},
                {"sample": "s1", "read_id": "readC", "exact_deletion_id": "mtDel_00100_00500_00399", "junction_id": "j1", "left_breakpoint": 101, "right_breakpoint": 499, "deleted_size": 397},
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            overlap, html_cells = sequence_remap_overlap_table(config, known_hits, junction_reads, Path(tmp))
            self.assertEqual(int(overlap.loc[0, "sequence_search_reads"]), 2)
            self.assertEqual(int(overlap.loc[0, "remap_nearby_reads"]), 2)
            self.assertEqual(int(overlap.loc[0, "shared_reads"]), 1)
            self.assertEqual(int(overlap.loc[0, "sequence_only_reads"]), 1)
            self.assertEqual(int(overlap.loc[0, "remap_only_reads"]), 1)
            self.assertIn((0, "shared_reads"), html_cells)
            self.assertTrue((Path(tmp) / "configured_vs_remap__s1__mtDNA_100_500__shared_reads.tsv").exists())
            shared = (Path(tmp) / "configured_vs_remap__s1__mtDNA_100_500__shared_reads.tsv").read_text(encoding="utf-8")
            self.assertIn("readA/1", shared)
            self.assertIn("readA", shared)

    def test_analysis_deduplicates_same_read_from_rotated_references(self):
        import pandas as pd

        reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "junction_id": "j1", "rotation_name": "normal", "left_anchor_length": 20, "right_anchor_length": 20},
                {"sample": "s1", "read_id": "readA", "junction_id": "j1", "rotation_name": "half", "left_anchor_length": 25, "right_anchor_length": 25},
                {"sample": "s1", "read_id": "readB", "junction_id": "j1", "rotation_name": "normal", "left_anchor_length": 20, "right_anchor_length": 20},
            ]
        )
        deduped = deduplicate_evidence_reads(reads)
        self.assertEqual(len(deduped), 2)
        self.assertEqual(deduped.loc[deduped["read_id"] == "readA", "rotation_name"].iloc[0], "half")

    def test_expected_transcript_mask_identifies_configured_annotations(self):
        import pandas as pd

        reads = pd.DataFrame(
            [
                {"junction_id": "j1", "junction_interpretation": "expected_transcript_junction"},
                {"junction_id": "j2", "junction_interpretation": "candidate_deletion_junction"},
            ]
        )
        self.assertEqual(expected_transcript_mask(reads).tolist(), [True, False])

    def test_known_deletion_match_is_config_driven(self):
        config = {
            "analysis": {
                "known_deletions": [
                    {
                        "name": "common",
                        "left_breakpoint": 8470,
                        "right_breakpoint": 13447,
                        "deleted_size": 4977,
                        "breakpoint_tolerance_bp": 10,
                        "size_tolerance_bp": 20,
                    }
                ]
            }
        }
        label, reason = known_deletion_match(8472, 13448, 4975, config, 16569)
        self.assertEqual(label, "common")
        self.assertEqual(reason, "configured_known_deletion_match")

    def test_known_deletion_match_infers_targets_from_sequence_search_coordinates(self):
        config = {
            "analysis": {
                "known_sequence_searches": [
                    {"id": "mtDNA_100_500", "name": "Configured deletion 100-500"}
                ]
            }
        }
        label, reason = known_deletion_match(101, 499, 397, config, 1000)
        self.assertEqual(label, "Configured deletion 100-500")
        self.assertEqual(reason, "configured_sequence_search_target_match")

    def test_sequence_search_target_does_not_duplicate_explicit_known_deletion(self):
        config = {
            "analysis": {
                "known_deletions": [
                    {
                        "name": "explicit common",
                        "left_breakpoint": 100,
                        "right_breakpoint": 500,
                        "deleted_size": 399,
                        "breakpoint_tolerance_bp": 10,
                        "size_tolerance_bp": 20,
                    }
                ],
                "known_sequence_searches": [
                    {"id": "mtDNA_101_501", "name": "overlapping configured search"}
                ],
            }
        }
        targets = configured_deletion_targets(config, 1000)
        self.assertEqual([target["name"] for target in targets], ["explicit common"])

    def test_parse_star_junction_line(self):
        line = "MT_circular\t1100\t+\tMT_circular\t1500\t+\t1\t0\t0\treadA\t1101\t5S25M\t1476\t24M6S"
        row = parse_star_junction_line(line, "sample1", "rat", 16313, 1000)
        self.assertEqual(row["left_breakpoint"], 100)
        self.assertEqual(row["right_breakpoint"], 500)
        self.assertEqual(row["deleted_size"], 399)
        self.assertEqual(row["read_id"], "readA")
        self.assertEqual(row["left_anchor_length"], 25)
        self.assertEqual(row["right_anchor_length"], 24)

    def test_parse_star_junction_line_skips_header(self):
        line = "chr_donorA\tbrkpt_donorA\tstrand_donorA\tchr_acceptorB\tbrkpt_acceptorB\tstrand_acceptorB\tjunction_type\trepeat_left_lenA\trepeat_right_lenB\tread_name\tstart_alnA\tcigar_alnA\tstart_alnB\tcigar_alnB"
        self.assertIsNone(parse_star_junction_line(line, "sample1", "rat", 16313, 1000))

    def test_aligned_bases_from_cigar(self):
        self.assertEqual(aligned_bases_from_cigar("5S25M"), 25)
        self.assertEqual(aligned_bases_from_cigar("10M1I10M2D5="), 25)

    def test_rat_sra_sample_derivation(self):
        row = {
            "Run": "SRR17380091",
            "Age": "18mo",
            "treatment": "GPA",
            "Biological_Replicate": "18G_replicate_5",
        }
        self.assertEqual(derive_age(row), "18mo")
        self.assertEqual(derive_treatment(row), "GPA")
        self.assertEqual(derive_replicate(row), "5")
        self.assertEqual(make_sample_id("rat_aging_muscle", row, "rat_{age}_{treatment}_{replicate}"), "rat_18mo_GPA_5")

    def test_prepare_reads_layout_normalization(self):
        self.assertEqual(normalized_layout({"layout": "single"}), "single")
        self.assertEqual(normalized_layout({"LibraryLayout": "SINGLE"}), "single")
        self.assertEqual(normalized_layout({"layout": "paired"}), "paired")
        self.assertEqual(normalized_layout({"LibraryLayout": "PAIRED"}), "paired")

    def test_classify_mt_alignment(self):
        mt_names = {"MT", "chrM"}
        self.assertEqual(classify_alignment(FakeRead("MT", 60, nh=1), mt_names, 20), "high_confidence_mt")
        self.assertEqual(classify_alignment(FakeRead("MT", 60, nh=2), mt_names, 20), "ambiguous_mt")
        self.assertEqual(classify_alignment(FakeRead("MT", 10, nh=1), mt_names, 20), "low_quality_mt")
        self.assertEqual(classify_alignment(FakeRead("1", 60, nh=1), mt_names, 20), "non_mt_primary")
        self.assertEqual(classify_alignment(FakeRead(unmapped=True), mt_names, 20), "unmapped")

    def test_whole_genome_selector_keeps_mt_primary_but_flags_nuclear_competitor(self):
        mt_names = {"MT", "chrM"}
        category, primary = classify_group(
            [FakeWholeGenomeRead("MT", mapq=60)],
            mt_names,
            min_mt_mapq=0,
            min_mt_aligned_fraction=0.5,
            ambiguous_mapq_below=10,
            competing_nuclear_aligned_fraction=0.5,
        )
        self.assertEqual(category, "mt_primary_best")
        self.assertEqual(primary.reference_name, "MT")

        category, _ = classify_group(
            [
                FakeWholeGenomeRead("MT", mapq=60),
                FakeWholeGenomeRead("1", mapq=60, secondary=True, query_alignment_length=90),
            ],
            mt_names,
            min_mt_mapq=0,
            min_mt_aligned_fraction=0.5,
            ambiguous_mapq_below=10,
            competing_nuclear_aligned_fraction=0.5,
        )
        self.assertEqual(category, "mt_primary_ambiguous")

    def test_whole_genome_selector_rejects_nuclear_primary_with_mt_competitor(self):
        mt_names = {"MT", "chrM"}
        category, primary = classify_group(
            [
                FakeWholeGenomeRead("1", mapq=60),
                FakeWholeGenomeRead("MT", mapq=60, secondary=True),
            ],
            mt_names,
            min_mt_mapq=0,
            min_mt_aligned_fraction=0.5,
            ambiguous_mapq_below=10,
            competing_nuclear_aligned_fraction=0.5,
        )
        self.assertEqual(category, "nuclear_primary_with_mt_competitor")
        self.assertEqual(primary.reference_name, "1")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGTNacgtn"), "nacgtNACGT")


if __name__ == "__main__":
    unittest.main()
