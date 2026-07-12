import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import yaml
import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_quality_evidence import (
    cigar_query_metrics,
    clean_read_id,
    collapse_physical_observations,
    parse_star_chimeric_line,
    query_pair_metrics,
)
from analyze_deletions import count_matrix, read_tsv_safe
from filter_quality_profile import filter_profile
from finalize_quality_evidence import finalize, gene_pair_label
from estimate_breakpoint_reference_support import (
    candidate_positions,
    eligible_cluster_samples,
    original_to_rotated_pos,
    reuse_reference_support,
)
from make_deletion_report import (
    circular_validation_section,
    experimental_design_section,
    method_concordance_section,
    quality_profile_section,
    reference_support_explanation,
)
from make_deliverables import (
    create_zip_archive,
    lightweight_report_html,
    package_light_quality_results,
    package_quality_results,
)
from make_quality_report_index import profile_rows
from plot_deletion_results import gene_pair_pca_enabled, rank_label_boxes_overlap


def star_line(
    left=8102,
    right=12936,
    left_contig="MT",
    right_contig="MT",
    strand="+",
    read_id="readA",
    left_cigar="25M24S",
    right_cigar="25S24M",
):
    fields = [
        left_contig,
        str(left),
        strand,
        right_contig,
        str(right),
        strand,
        "0",
        "0",
        "0",
        read_id,
        str(left + 1),
        left_cigar,
        str(right - 23),
        right_cigar,
        "1",
        "49",
        "25",
        "45",
        "45",
        "0",
    ]
    return "\t".join(fields)


def parse(line, mt_length=16313):
    return parse_star_chimeric_line(
        line,
        "rat_sample",
        "rat",
        {"MT", "chrM"},
        mt_length,
        min_anchor=12,
        min_deletion_size=100,
        max_deletion_size=16000,
        max_query_overlap=20,
        max_query_gap=20,
        require_same_orientation=True,
    )


class QualityEvidenceTests(unittest.TestCase):
    def test_quality_deliverables_package_all_profiles_and_valid_links(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "results" / "test_dataset"
            shared = root / "quality" / "shared"
            shared.mkdir(parents=True)
            config = {
                "dataset": {"name": "test_dataset", "title": "Test dataset"},
                "quality": {
                    "primary_report_profile": "standard",
                    "report_profiles": {
                        "stringent": {"include_tiers": ["strong"]},
                        "standard": {"include_tiers": ["strong", "supported"]},
                        "exploratory": {"include_tiers": ["strong", "supported", "review"]},
                    },
                },
            }
            membership_rows = [
                "exact_deletion_id\treport_profile\tincluded\tdistinct_observation_count",
                "mtDel_1\tstringent\tyes\t2",
                "mtDel_1\tstandard\tyes\t2",
                "mtDel_1\texploratory\tyes\t2",
            ]
            (shared / "report_profile_membership.tsv").write_text("\n".join(membership_rows) + "\n", encoding="utf-8")
            (shared / "canonical_clusters.tsv").write_text("exact_deletion_id\nmtDel_1\n", encoding="utf-8")
            (shared / "source_candidates.tsv").write_text("read_id\nr1\n", encoding="utf-8")
            (shared / "evidence_build_summary.tsv").write_text("metric\tvalue\nrows\t1\n", encoding="utf-8")
            (shared / "quality_tier_summary.tsv").write_text("tier\tclusters\nstrong\t1\n", encoding="utf-8")
            (shared / "resolved_quality_config.yaml").write_text("quality:\n  enabled: true\n", encoding="utf-8")

            profiles = ["stringent", "standard", "exploratory"]
            for profile in profiles:
                source = root / "quality" / "profiles" / profile
                for relative in (".report/read_lists", "plots", "matrices", "junctions", "analysis"):
                    (source / relative).mkdir(parents=True)
                (source / ".report" / "index.html").write_text(
                    '<main><a href="plots/example.pdf">Plot</a>'
                    '<a class="read-list-link" href="read_lists/manifest.tsv">Reads</a></main>',
                    encoding="utf-8",
                )
                (source / ".report" / "read_lists" / "manifest.tsv").write_text("id\nmtDel_1\n", encoding="utf-8")
                (source / "plots" / "example.pdf").write_bytes(b"%PDF-1.4\n")
                (source / "matrices" / "exact_deletion_raw_counts.tsv").write_text("sample\tmtDel_1\ns1\t2\n", encoding="utf-8")
                (source / "junctions" / "junction_clusters.tsv").write_text("exact_deletion_id\nmtDel_1\n", encoding="utf-8")
                (source / "junctions" / "canonical_observations.tsv").write_text("exact_deletion_id\nmtDel_1\n", encoding="utf-8")
                (source / "junctions" / "junction_id_map.tsv").write_text("exact_deletion_id\nmtDel_1\n", encoding="utf-8")
                (source / "analysis" / "deletion_burden.tsv").write_text("sample\tdeletion_reads\ns1\t2\n", encoding="utf-8")

            package = root / "test_dataset_deliverables"
            package_quality_results(root, package, "test_dataset", config, profiles)

            selector = (package / "index.html").read_text(encoding="utf-8")
            for profile in profiles:
                self.assertIn(f'profiles/{profile}/index.html', selector)
                self.assertTrue((package / "profiles" / profile / "index.html").is_file())
                self.assertTrue((package / "profiles" / profile / "plots" / "example.pdf").is_file())
                self.assertTrue((package / "profiles" / profile / "tables" / "exact_deletions.tsv").is_file())
                self.assertTrue((package / "profiles" / profile / "matrices" / "exact_deletion_raw_counts.tsv").is_file())
                self.assertTrue((package / "profiles" / profile / "read_lists" / "manifest.tsv").is_file())
            self.assertEqual(
                (package / "shared" / "canonical_clusters.tsv").read_text(encoding="utf-8"),
                (shared / "canonical_clusters.tsv").read_text(encoding="utf-8"),
            )

            light = root / "test_dataset_deliverables_light"
            package_light_quality_results(root, light, "test_dataset", config, profiles)
            (light / "DELIVERABLES_COMPLETE.txt").write_text("complete\n", encoding="utf-8")
            archive = root / "test_dataset_deliverables_light.zip"
            create_zip_archive(light, archive)

            self.assertFalse((light / "shared" / "source_candidates.tsv").exists())
            self.assertFalse((light / "shared" / "canonical_clusters.tsv").exists())
            self.assertTrue((light / "shared" / "quality_tier_summary.tsv").is_file())
            for profile in profiles:
                profile_root = light / "profiles" / profile
                self.assertFalse((profile_root / "read_lists").exists())
                self.assertFalse((profile_root / "tables" / "canonical_observations.tsv").exists())
                self.assertTrue((profile_root / "tables" / "exact_deletions.tsv").is_file())
                report = (profile_root / "index.html").read_text(encoding="utf-8")
                self.assertIn("Light deliverable package", report)
                self.assertNotIn('href="read_lists/', report)
            with zipfile.ZipFile(archive, "r") as zipped:
                self.assertIsNone(zipped.testzip())
                self.assertIn(
                    "test_dataset_deliverables_light/profiles/standard/index.html",
                    zipped.namelist(),
                )

    def test_lightweight_report_rejects_unclassified_read_list_links(self):
        with self.assertRaises(ValueError):
            lightweight_report_html('<main><a href="read_lists/unclassified.tsv">Reads</a></main>')

    def test_analysis_accepts_empty_audit_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.tsv"
            path.write_text("", encoding="utf-8")
            self.assertTrue(read_tsv_safe(str(path)).empty)

    def test_dual_caller_is_explicit_for_short_read_rna_only(self):
        root = Path(__file__).resolve().parents[1]
        enabled = []
        for path in sorted((root / "config" / "datasets").glob("*.yaml")):
            config = yaml.safe_load(path.read_text(encoding="utf-8"))
            dual = bool(config.get("quality", {}).get("short_read_rna_dual_caller", {}).get("enabled", False))
            if dual:
                enabled.append(config["dataset"]["name"])
                self.assertEqual(config["dataset"]["read_technology"], "illumina")
                self.assertEqual(config["dataset"]["molecule_type"], "rna")
                self.assertIn("short_read", config["dataset"]["library_strategy"])
        self.assertEqual(
            enabled,
            ["human_bulkseq_matched_nanopore", "human_common_deletion", "rat_aging_muscle"],
        )

    def test_documentation_describes_current_workflow_without_change_history(self):
        root = Path(__file__).resolve().parents[1]
        documentation = "\n".join(
            (root / path).read_text(encoding="utf-8").lower()
            for path in [
                "README.md",
                "docs/directed_circular_deletion_workflow.md",
                "docs/workflow_methods_and_assumptions.md",
            ]
        )
        for artifact in [
            "legacy workflow",
            "previous workflow",
            "still summarizes",
            "what went wrong",
            "during development",
            "phase 2 currently",
            "worth evaluating",
        ]:
            self.assertNotIn(artifact, documentation)
        self.assertIn("3.0-quality-evidence-multi-caller", documentation)
        self.assertIn("results/<dataset>/quality/report/index.html", documentation)
        self.assertIn("results/<dataset>/<dataset>_deliverables/index.html", documentation)
        self.assertIn("results/<dataset>/<dataset>_deliverables_light.zip", documentation)
        self.assertIn("not_available_from_retained_intermediates", documentation)

    def test_cigar_query_metrics_tracks_short_read_coordinates(self):
        self.assertEqual(
            cigar_query_metrics("25M24S"),
            {
                "query_length": 49,
                "query_start": 0,
                "query_end": 25,
                "aligned_query_length": 25,
                "aligned_reference_length": 25,
            },
        )
        self.assertEqual(cigar_query_metrics("25S24M")["query_start"], 25)

    def test_query_union_coverage_does_not_count_overlap_twice(self):
        metrics = query_pair_metrics(
            {"query_start": 0, "query_end": 30, "query_length": 49},
            {"query_start": 20, "query_end": 49, "query_length": 49},
        )
        self.assertEqual(metrics["query_overlap_bp"], 10)
        self.assertAlmostEqual(metrics["query_union_coverage"], 1.0)
        self.assertEqual(metrics["query_union_aligned_length"], 49)
        self.assertEqual(metrics["total_aligned_query_length"], 59)

    def test_star_parser_accepts_49bp_single_end_mt_junction(self):
        row = parse(star_line())
        self.assertIsNotNone(row)
        self.assertEqual(row["filter_status"], "pass")
        self.assertEqual(row["read_length"], 49)
        self.assertEqual(row["left_anchor_length"], 25)
        self.assertEqual(row["right_anchor_length"], 24)
        self.assertEqual(row["query_gap_bp"], 0)
        self.assertEqual(row["physical_observation_id"], "readA")
        self.assertEqual(row["evidence_source"], "star_chimeric")

    def test_star_parser_accepts_twelve_base_short_read_anchor(self):
        row = parse(star_line(left_cigar="12M37S", right_cigar="12S37M"))
        self.assertEqual(row["filter_status"], "pass")
        self.assertEqual(row["min_anchor_length"], 12)

    def test_star_parser_rejects_non_mt_pair(self):
        self.assertIsNone(parse(star_line(right_contig="1")))

    def test_star_gene_anchor_filter_rejects_same_gene_and_missing_gene(self):
        same_gene = parse_star_chimeric_line(
            star_line(),
            "rat_sample",
            "rat",
            {"MT"},
            16313,
            12,
            100,
            16000,
            20,
            20,
            True,
            [{"name": "Mt-long", "start": 8000, "end": 13000}],
            {},
            True,
            True,
        )
        self.assertEqual(same_gene["filter_status"], "fail")
        self.assertIn("same_gene_chimeric_alignment", same_gene["filter_reason"])

        missing = parse_star_chimeric_line(
            star_line(),
            "rat_sample",
            "rat",
            {"MT"},
            16313,
            12,
            100,
            16000,
            20,
            20,
            True,
            [],
            {},
            True,
            True,
        )
        self.assertIn("missing_annotated_gene_anchor", missing["filter_reason"])

    def test_star_gene_anchor_filter_accepts_distinct_genes(self):
        row = parse_star_chimeric_line(
            star_line(),
            "rat_sample",
            "rat",
            {"MT"},
            16313,
            12,
            100,
            16000,
            20,
            20,
            True,
            [
                {"name": "Mt-left", "start": 8000, "end": 8500},
                {"name": "Mt-right", "start": 12800, "end": 13500},
            ],
            {},
            True,
            True,
        )
        self.assertEqual(row["filter_status"], "pass")
        self.assertEqual(row["star_gene_pair_label"], "Mt-left--Mt-right")

    def test_star_parser_preserves_directed_origin_spanning_arc(self):
        row = parse(star_line(left=15000, right=1000))
        self.assertEqual(row["filter_status"], "pass")
        self.assertEqual(row["wraps_origin"], "yes")
        self.assertEqual(row["deleted_size"], 2312)
        self.assertEqual(row["left_breakpoint"], 15000)
        self.assertEqual(row["right_breakpoint"], 1000)

    def test_clean_read_id_deduplicates_mates_as_one_fragment(self):
        self.assertEqual(clean_read_id("fragment/1"), "fragment")
        self.assertEqual(clean_read_id("fragment/2"), "fragment")
        self.assertEqual(clean_read_id("single_end_read"), "single_end_read")

    def test_same_physical_read_from_two_callers_is_one_observation(self):
        star = parse(star_line(left=8102, right=12936, read_id="shared"))
        minimap = dict(star)
        minimap.update(
            {
                "left_breakpoint": 8104,
                "right_breakpoint": 12934,
                "evidence_source": "minimap2_remap",
                "source": "minimap2_split_alignment",
                "rotation_name": "normal",
                "min_mapq": 10,
                "primary_chain_evidence": "yes",
                "secondary_only_evidence": "no",
            }
        )
        observations = collapse_physical_observations([star, minimap], slop=5, mt_length=16313)
        self.assertEqual(len(observations), 1)
        self.assertEqual(observations[0]["both_callers_support"], "yes")
        self.assertEqual(observations[0]["source_record_count"], 2)
        self.assertEqual(observations[0]["physical_observation_id"], "shared")
        self.assertEqual(observations[0]["star_gene_pair_label"], "")

    def test_physical_observation_preserves_layout_and_segment_metrics(self):
        base = parse(star_line())
        base.update(
            {
                "library_layout": "paired",
                "observation_unit": "sequenced_fragment",
                "paired_end_collapsed": "yes",
                "mate_context_status": "not_available_from_retained_intermediates",
                "both_mates_support_same_cluster": "not_available_from_retained_intermediates",
                "non_supporting_mate_near_breakpoint": "not_available_from_retained_intermediates",
                "mate_mapping_class": "not_available_from_retained_intermediates",
                "mate_placement_discordant": "not_available_from_retained_intermediates",
            }
        )
        observation = collapse_physical_observations([base], 10, 16313)[0]
        self.assertEqual(observation["library_layout"], "paired")
        self.assertEqual(observation["observation_unit"], "sequenced_fragment")
        self.assertEqual(observation["paired_end_collapsed"], "yes")
        self.assertEqual(observation["left_cigar"], "25M24S")
        self.assertEqual(observation["right_cigar"], "25S24M")
        self.assertEqual(observation["left_mapq"], "")
        self.assertEqual(observation["left_alignment_score"], 45.0)
        self.assertEqual(observation["strand"], "+:+")
        self.assertIn("star_chimeric", observation["source_coordinate_pairs"])

    def test_reciprocal_hypotheses_are_not_collapsed(self):
        forward = parse(star_line(left=8102, right=12936, read_id="ambiguous"))
        reverse = parse(star_line(left=12936, right=8102, read_id="ambiguous"))
        observations = collapse_physical_observations([forward, reverse], slop=5, mt_length=16313)
        self.assertEqual(len(observations), 2)
        self.assertTrue(all(row["deletion_hypotheses_from_read"] == 2 for row in observations))

    def test_gene_pair_is_derived_without_changing_exact_identity(self):
        cluster = {
            "exact_deletion_id": "mtDel_08102_12936_04833",
            "nearest_left_feature": "Mt-atp8-6:protein_coding:7800-8500;distance=0",
            "nearest_right_feature": "Mt-cyb:protein_coding:14000-15100;distance=10",
        }
        self.assertEqual(gene_pair_label(cluster), "Mt-atp8-6--Mt-cyb")
        self.assertEqual(cluster["exact_deletion_id"], "mtDel_08102_12936_04833")

    def test_star_gene_pair_is_preserved_separately_from_flanking_pair(self):
        clusters = [
            {
                "exact_deletion_id": "event",
                "junction_id": "event",
                "nearest_left_feature": "flank-left:gene:1-10",
                "nearest_right_feature": "flank-right:gene:20-30",
            }
        ]
        observations = [
            {
                "exact_deletion_id": "event",
                "junction_id": "event",
                "sample": "s1",
                "read_id": "r1",
                "star_gene_pair_label": "star-left--star-right",
                "star_support": "yes",
                "minimap2_support": "no",
                "evidence_sources": "star_chimeric",
                "direction_status": "directed",
            }
        ]
        final_clusters, final_observations, _, _ = finalize(clusters, observations, {})
        self.assertEqual(final_clusters[0]["gene_pair_label"], "star-left--star-right")
        self.assertEqual(final_clusters[0]["breakpoint_flanking_gene_pair_label"], "flank-left--flank-right")
        self.assertEqual(final_observations[0]["gene_pair_label"], "star-left--star-right")

    def test_reference_support_uses_only_minimap_eligible_clusters_and_samples(self):
        clusters = pd.DataFrame(
            [
                {"exact_deletion_id": "star", "left_breakpoint": 10, "right_breakpoint": 20},
                {"exact_deletion_id": "remap", "left_breakpoint": 30, "right_breakpoint": 40},
            ]
        )
        observations = pd.DataFrame(
            [
                {"exact_deletion_id": "star", "sample": "s1", "minimap2_support": "no"},
                {"exact_deletion_id": "remap", "sample": "s2", "minimap2_support": "yes"},
            ]
        )
        eligible = eligible_cluster_samples(clusters, observations)
        self.assertNotIn("star", eligible)
        self.assertEqual(eligible["remap"], {"s2"})
        self.assertEqual(candidate_positions(clusters, eligible), {"s2": {30, 40}})
        self.assertEqual(original_to_rotated_pos(51, 100, 51), 1)

    def test_existing_reference_support_is_reused_without_assigning_star_only_denominator(self):
        clusters = pd.DataFrame(
            [
                {"exact_deletion_id": "cached"},
                {"exact_deletion_id": "star-only"},
            ]
        )
        observations = pd.DataFrame(
            [
                {"exact_deletion_id": "cached", "sample": "s1", "minimap2_support": "yes"},
                {"exact_deletion_id": "star-only", "sample": "s1", "minimap2_support": "no"},
            ]
        )
        existing = pd.DataFrame(
            [
                {
                    "exact_deletion_id": "cached",
                    "left_reference_spanning_reads": 10,
                    "reference_support_method": "primary_alignment_depth_max_across_rotations",
                }
            ]
        )
        reused = reuse_reference_support(clusters, observations, existing, 20).set_index("exact_deletion_id")
        self.assertEqual(reused.loc["cached", "left_reference_spanning_reads"], 10)
        self.assertTrue(reused.loc["cached", "reference_support_method"].startswith("reused_existing_"))
        self.assertEqual(
            reused.loc["star-only", "reference_support_method"],
            "not_available_without_minimap2_remap_evidence",
        )

    def test_quality_tiers_use_independent_observation_support(self):
        clusters = [
            {
                "exact_deletion_id": "event",
                "junction_id": "event",
                "nearest_left_feature": "Mt-co1:gene:1-10",
                "nearest_right_feature": "Mt-cyb:gene:20-30",
                "junction_interpretation": "candidate_deletion_junction",
            }
        ]
        base = {
            "exact_deletion_id": "event",
            "junction_id": "event",
            "evidence_sources": "star_chimeric",
            "star_support": "yes",
            "minimap2_support": "no",
            "both_callers_support": "no",
            "source_record_count": 1,
            "min_anchor_length": 12,
            "query_segments_adjacent": "yes",
            "deletion_hypotheses_from_read": 1,
            "direction_status": "directed",
        }
        one = [dict(base, sample="s1", read_id="r1")]
        final_clusters, _, _, _ = finalize(clusters, one, {})
        self.assertEqual(final_clusters[0]["quality_tier"], "review")

        cross_sample = one + [dict(base, sample="s2", read_id="r2")]
        final_clusters, _, _, _ = finalize(clusters, cross_sample, {})
        self.assertEqual(final_clusters[0]["quality_tier"], "supported")

        replicated = one + [dict(base, sample="s1", read_id="r2")]
        final_clusters, _, _, _ = finalize(clusters, replicated, {})
        self.assertEqual(final_clusters[0]["quality_tier"], "strong")

    def test_cross_caller_detection_does_not_turn_one_read_into_two(self):
        clusters = [
            {
                "exact_deletion_id": "event",
                "junction_id": "event",
                "nearest_left_feature": "Mt-co1:gene:1-10",
                "nearest_right_feature": "Mt-cyb:gene:20-30",
                "junction_interpretation": "candidate_deletion_junction",
            }
        ]
        observations = [
            {
                "exact_deletion_id": "event",
                "junction_id": "event",
                "sample": "s1",
                "read_id": "one_read",
                "evidence_sources": "minimap2_remap;star_chimeric",
                "star_support": "yes",
                "minimap2_support": "yes",
                "both_callers_support": "yes",
                "source_record_count": 2,
                "direction_status": "directed",
            }
        ]
        final_clusters, _, _, _ = finalize(clusters, observations, {})
        self.assertEqual(final_clusters[0]["distinct_observation_count"], 1)
        self.assertEqual(final_clusters[0]["both_caller_supporting_observations"], 1)
        self.assertEqual(final_clusters[0]["quality_tier"], "review")

    def test_expected_transcript_cluster_is_rejected_when_configured(self):
        clusters = [
            {
                "exact_deletion_id": "transcript",
                "junction_id": "transcript",
                "junction_interpretation": "expected_transcript_junction",
            }
        ]
        observations = [
            {
                "exact_deletion_id": "transcript",
                "junction_id": "transcript",
                "sample": "s1",
                "read_id": "r1",
                "direction_status": "directed",
            },
            {
                "exact_deletion_id": "transcript",
                "junction_id": "transcript",
                "sample": "s1",
                "read_id": "r2",
                "direction_status": "directed",
            },
        ]
        final_clusters, _, _, _ = finalize(
            clusters,
            observations,
            {"junctions": {"exclude_expected_transcript_junctions": True}},
        )
        self.assertEqual(final_clusters[0]["quality_tier"], "rejected")

    def test_report_profiles_filter_shared_stable_ids(self):
        clusters = [
            {"exact_deletion_id": "strong_event", "quality_tier": "strong"},
            {"exact_deletion_id": "supported_event", "quality_tier": "supported"},
            {"exact_deletion_id": "review_event", "quality_tier": "review"},
        ]
        observations = [
            {"exact_deletion_id": row["exact_deletion_id"], "read_id": f"read_{index}"}
            for index, row in enumerate(clusters)
        ]
        id_map = [dict(row) for row in observations]
        standard = filter_profile(clusters, observations, id_map, "standard", {})
        exploratory = filter_profile(clusters, observations, id_map, "exploratory", {})
        self.assertEqual([row["exact_deletion_id"] for row in standard[0]], ["strong_event", "supported_event"])
        self.assertEqual([row["exact_deletion_id"] for row in exploratory[0]], ["strong_event", "supported_event", "review_event"])
        self.assertEqual(len({row["exact_deletion_id"] for row in exploratory[0]}), 3)

    def test_profile_matrix_contains_only_profile_observations(self):
        clusters = [
            {"exact_deletion_id": "strong_event", "quality_tier": "strong"},
            {"exact_deletion_id": "review_event", "quality_tier": "review"},
        ]
        observations = [
            {"exact_deletion_id": "strong_event", "sample": "s1", "read_id": "strong_read"},
            {"exact_deletion_id": "review_event", "sample": "s1", "read_id": "review_read"},
        ]
        _, standard_observations, _ = filter_profile(clusters, observations, observations, "standard", {})
        matrix = count_matrix(
            pd.DataFrame({"sample": ["s1"], "group": ["g"]}),
            pd.DataFrame(standard_observations),
            "exact_deletion_id",
        )
        self.assertIn("strong_event", matrix.columns)
        self.assertNotIn("review_event", matrix.columns)
        self.assertEqual(matrix.loc[0, "strong_event"], 1)

    def test_report_index_counts_stable_cluster_membership(self):
        config = {
            "quality": {
                "primary_report_profile": "standard",
                "report_profiles": {
                    "stringent": {"include_tiers": ["strong"]},
                    "standard": {"include_tiers": ["strong", "supported"]},
                },
            }
        }
        membership = [
            {"exact_deletion_id": "a", "report_profile": "stringent", "included": "yes", "distinct_observation_count": "2"},
            {"exact_deletion_id": "a", "report_profile": "standard", "included": "yes", "distinct_observation_count": "2"},
            {"exact_deletion_id": "b", "report_profile": "standard", "included": "yes", "distinct_observation_count": "3"},
        ]
        rows = profile_rows(config, membership)
        by_name = {row["profile"]: row for row in rows}
        self.assertEqual(by_name["stringent"]["deletion_clusters"], 1)
        self.assertEqual(by_name["standard"]["deletion_clusters"], 2)
        self.assertEqual(by_name["standard"]["distinct_observations"], 5)
        self.assertEqual(by_name["standard"]["role"], "Primary interpretation")

    def test_quality_profile_report_states_profile_specific_pca(self):
        config = {
            "quality": {
                "report_profiles": {
                    "standard": {"include_tiers": ["strong", "supported"]},
                }
            }
        }
        clusters = pd.DataFrame(
            {
                "exact_deletion_id": ["a"],
                "evidence_status": ["star_and_minimap2"],
            }
        )
        observations = pd.DataFrame({"read_id": ["r1", "r2"]})
        rendered = quality_profile_section(config, "standard", clusters, observations)
        self.assertIn("strong, supported", rendered)
        self.assertIn("axes are not assumed equivalent across profiles", rendered)
        self.assertIn("star_and_minimap2", rendered)

    def test_report_concordance_states_star_filter_boundary(self):
        rendered = method_concordance_section(
            {"quality": {"short_read_rna_dual_caller": {"enabled": True}}},
            "standard",
            pd.DataFrame(
                {
                    "evidence_source": ["star_chimeric", "star_chimeric"],
                    "filter_status": ["pass", "fail"],
                    "filter_reason": ["", "same_gene_chimeric_alignment"],
                }
            ),
            pd.DataFrame({"evidence_status": ["star_only"]}),
            pd.DataFrame({"evidence_sources": ["star_chimeric"]}),
        )
        self.assertIn("STAR-Fusion is not run", rendered)
        self.assertIn("same_gene_chimeric_alignment", rendered)

    def test_combined_report_layout_note_names_both_callers_without_history_language(self):
        rendered = experimental_design_section(
            pd.DataFrame({"sample": ["s1"], "layout": ["single"], "group": ["g"]}),
            pd.DataFrame(),
            "group",
            {"quality": {"short_read_rna_dual_caller": {"enabled": True}}},
            "standard",
        )
        self.assertIn("STAR chimeric or mitochondrial-remap alignments", rendered)
        self.assertIn("No mate exists, so mate evidence is neither required nor used", rendered)
        self.assertNotIn("mate-pair consistency is unavailable", rendered)
        self.assertNotIn("still", rendered.lower())

    def test_paired_report_states_unavailable_mate_context(self):
        rendered = experimental_design_section(
            pd.DataFrame({"sample": ["s1"], "layout": ["paired"], "group": ["g"]}),
            pd.DataFrame(),
            "group",
            {"quality": {"short_read_rna_dual_caller": {"enabled": True}}},
            "standard",
        )
        self.assertIn("Mate-placement context is marked unavailable", rendered)

    def test_gene_pair_pca_requires_explicit_dual_caller_configuration(self):
        self.assertTrue(gene_pair_pca_enabled({"quality": {"short_read_rna_dual_caller": {"enabled": True}}}))
        self.assertFalse(gene_pair_pca_enabled({"quality": {"short_read_rna_dual_caller": {"enabled": False}}}))
        self.assertFalse(gene_pair_pca_enabled({}))

    def test_reference_support_text_distinguishes_star_only_calls(self):
        clusters = pd.DataFrame(
            {
                "left_reference_spanning_reads": [1],
                "right_reference_spanning_reads": [2],
                "reference_spanning_reads_min": [1],
                "local_split_support_fraction": [0.5],
                "evidence_status": ["star_only"],
            }
        )
        rendered = reference_support_explanation(clusters)
        self.assertIn("minimap2 deletion-supporting split reads", rendered)
        self.assertIn("STAR-only exact deletions are marked unavailable", rendered)

        minimap_rendered = reference_support_explanation(clusters.assign(evidence_status="minimap2_only"))
        self.assertNotIn("STAR-only exact deletions are marked unavailable", minimap_rendered)

    def test_combined_circular_checks_describe_both_coordinate_sources(self):
        rendered = circular_validation_section(
            {"quality": {"short_read_rna_dual_caller": {"enabled": True}}},
            pd.DataFrame(),
        )
        self.assertIn("Minimap2 mitochondrial-remap observations", rendered)
        self.assertIn("STAR chimeric observations", rendered)

    def test_rank_label_collision_prefers_omitting_overlapping_label(self):
        self.assertTrue(rank_label_boxes_overlap((10, 10, 4), (14, 10, 4)))
        self.assertFalse(rank_label_boxes_overlap((10, 10, 4), (20, 10, 4)))


if __name__ == "__main__":
    unittest.main()
