import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from classify_mt_reads import classify_alignment, reverse_complement
from annotate_junctions import append_configured_regions, apply_feature_aliases, biological_features
from analyze_deletions import deduplicate_evidence_reads, expected_transcript_mask
from circular_deletions import affected_feature_impact, configured_deletion_targets, known_deletion_match, normalize_pos as normalize_rotated_pos
from cluster_junctions import canonical_junction
from consolidate_deletions import cluster_rows
from make_rotated_mt_reference import rotate_sequence
from parse_split_alignments import aligned_bases_from_cigar, deletion_size, normalize_pos, parse_star_junction_line
from prepare_reads import normalized_layout
from plot_deletion_results import (
    apply_cluster_coordinates,
    draw_feature_track_axis,
    location_features,
    mitochondrial_axis_bounds,
    pooled_endpoint_density,
    rainfall_point_sizes,
    support_scale_limits,
    rainfall_support_limits,
    rainfall_y_axis_min,
    prepare_location_plot_data,
    support_legend_values,
    support_size_legend_values,
    value_columns,
)
from estimate_breakpoint_reference_support import circular_window, window_covered
from make_deletion_report import (
    exact_deletion_display_table,
    exact_deletion_support_read_links,
    sequence_remap_overlap_table,
    table_html,
    write_configured_sequence_read_lists,
    write_exact_deletion_read_lists,
)
from resolve_samples import derive_age, derive_replicate, derive_treatment, make_sample_id
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


class CoreTests(unittest.TestCase):
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

    def test_canonical_junction_merges_reciprocal_breakpoint_pair(self):
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
        self.assertEqual((reciprocal["left_breakpoint"], reciprocal["right_breakpoint"], reciprocal["deleted_size"]), (8472, 13448, 4975))
        self.assertEqual(reciprocal["canonical_orientation"], "reversed_to_shorter_interval")

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

    def test_cluster_coordinate_adapter_recomputes_size_from_representative_breakpoints(self):
        reads = pd.DataFrame(
            [
                {"sample": "s1", "read_id": "readA", "junction_id": "j1", "left_breakpoint": 100, "right_breakpoint": 200, "deleted_size": 99},
                {"sample": "s1", "read_id": "readB", "junction_id": "j1", "left_breakpoint": 110, "right_breakpoint": 210, "deleted_size": 99},
            ]
        )
        clusters = pd.DataFrame({"junction_id": ["j1"], "left_breakpoint": [105], "right_breakpoint": [207], "deleted_size": [999]})
        corrected = apply_cluster_coordinates(reads, clusters, mt_length=1000)
        self.assertEqual(corrected["left_breakpoint"].tolist(), [105, 105])
        self.assertEqual(corrected["right_breakpoint"].tolist(), [207, 207])
        self.assertEqual(corrected["deleted_size"].tolist(), [101, 101])

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

    def test_rainfall_y_axis_min_tracks_deletion_cutoff(self):
        self.assertEqual(rainfall_y_axis_min(pd.Series([104, 500, 5000])), 100)
        self.assertEqual(rainfall_y_axis_min(pd.Series([12, 500])), 10)

    def test_report_table_html_wraps_tables_for_horizontal_scroll(self):
        import pandas as pd

        html = table_html(pd.DataFrame([{"wide_column_name": "value"}]))
        self.assertIn('class="table-wrap"', html)
        self.assertIn('class="dataframe data-table"', html)
        self.assertEqual(html.count('class="table-wrap"'), 1)

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
