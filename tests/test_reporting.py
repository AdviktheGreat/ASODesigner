"""Tests for the reporting module: data structures, JSON, HTML, assembly."""

import json

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.config_loader import load_config
from aso_designer.reporting import (
    ReportAssembler,
    ReportData,
    ReportSummary,
    assemble_report,
    assemble_report_from_config,
    build_report_data,
    compute_summary_statistics,
    generate_html_report,
    generate_json_report,
    generate_report,
    generate_report_from_ranking,
    generate_reports_batch,
    load_json_report,
)
from aso_designer.scoring import (
    CandidateScore,
    CompositeScore,
    ScoreResult,
    score_rank_and_filter,
)
from aso_designer.thermodynamics import ThermodynamicSummary

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_cs(score: float, seq: str = "ATGC", rank: int | None = None) -> CandidateScore:
    """Build a CandidateScore with a given final score."""
    cs = CandidateScore(
        sequence=seq,
        start=0,
        end=20,
        chemistry="gapmer",
        composite_score=CompositeScore(
            component_results=[
                ScoreResult("gc_content", score, 1.0, score),
            ],
            composite_score=score,
            chemistry_bonus=0.0,
            final_score=score,
        ),
        final_score=score,
        rank=rank,
    )
    return cs


# ---------------------------------------------------------------------------
# ReportSummary
# ---------------------------------------------------------------------------


class TestReportSummary:
    def test_defaults(self):
        s = ReportSummary()
        assert s.target_name == ""
        assert s.target_id == ""
        assert s.total_candidates == 0
        assert s.top_score == 0.0
        assert s.chemistry == "gapmer"

    def test_custom_values(self):
        s = ReportSummary(
            target_name="BRCA1",
            target_id="NM_007294",
            total_candidates=100,
            top_score=0.95,
            chemistry="mixmer",
        )
        assert s.target_name == "BRCA1"
        assert s.total_candidates == 100
        assert s.top_score == pytest.approx(0.95)
        assert s.chemistry == "mixmer"

    def test_to_dict(self):
        s = ReportSummary(target_name="TP53", total_candidates=50, top_score=0.9)
        d = s.to_dict()
        assert d["target_name"] == "TP53"
        assert d["total_candidates"] == 50
        assert d["top_score"] == pytest.approx(0.9)

    def test_round_trip(self):
        s = ReportSummary(
            target_name="EGFR",
            target_id="NM_005228",
            total_candidates=200,
            passed_filters=150,
            scored_candidates=100,
            ranked_candidates=20,
            top_score=0.92,
            median_score=0.7,
            mean_score=0.65,
            min_score=0.3,
            chemistry="2ome",
            timestamp="2026-10-03T12:00:00Z",
        )
        d = s.to_dict()
        s2 = ReportSummary.from_dict(d)
        assert s2.target_name == "EGFR"
        assert s2.target_id == "NM_005228"
        assert s2.total_candidates == 200
        assert s2.top_score == pytest.approx(0.92)
        assert s2.chemistry == "2ome"


# ---------------------------------------------------------------------------
# ReportData
# ---------------------------------------------------------------------------


class TestReportData:
    def test_defaults(self):
        r = ReportData()
        assert r.summary.target_name == ""
        assert r.candidates == []
        assert r.thermodynamic_details == []
        assert r.off_target_details == []
        assert r.variant_details == []
        assert r.metadata == {}

    def test_with_data(self):
        s = ReportSummary(target_name="BRCA1", total_candidates=10)
        r = ReportData(
            summary=s,
            candidates=[{"sequence": "ATGC", "final_score": 0.9}],
            metadata={"version": "1.0"},
        )
        assert r.summary.target_name == "BRCA1"
        assert len(r.candidates) == 1
        assert r.metadata["version"] == "1.0"

    def test_round_trip(self):
        s = ReportSummary(target_name="TP53", total_candidates=50, top_score=0.9)
        r = ReportData(
            summary=s,
            candidates=[{"sequence": "ATGC", "final_score": 0.9, "rank": 1}],
            thermodynamic_details=[{"tm": 55.0, "binding_class": "strong"}],
            metadata={"version": "1.0"},
        )
        d = r.to_dict()
        r2 = ReportData.from_dict(d)
        assert r2.summary.target_name == "TP53"
        assert len(r2.candidates) == 1
        assert len(r2.thermodynamic_details) == 1
        assert r2.metadata["version"] == "1.0"


# ---------------------------------------------------------------------------
# compute_summary_statistics
# ---------------------------------------------------------------------------


class TestComputeSummaryStatistics:
    def test_basic(self):
        scored = [
            _make_cs(0.9, rank=1),
            _make_cs(0.7, rank=2),
            _make_cs(0.5, rank=3),
        ]
        s = compute_summary_statistics(
            scored,
            total_candidates=100,
            passed_filters=80,
            target_name="BRCA1",
            target_id="NM_007294",
        )
        assert s.target_name == "BRCA1"
        assert s.total_candidates == 100
        assert s.passed_filters == 80
        assert s.scored_candidates == 3
        assert s.ranked_candidates == 3
        assert s.top_score == pytest.approx(0.9)
        assert s.min_score == pytest.approx(0.5)
        assert s.mean_score == pytest.approx(0.7)
        assert s.median_score == pytest.approx(0.7)
        assert s.timestamp != ""

    def test_empty(self):
        s = compute_summary_statistics([], total_candidates=0)
        assert s.scored_candidates == 0
        assert s.top_score == pytest.approx(0.0)
        assert s.median_score == pytest.approx(0.0)

    def test_median_even(self):
        scored = [
            _make_cs(0.9),
            _make_cs(0.7),
            _make_cs(0.5),
            _make_cs(0.3),
        ]
        s = compute_summary_statistics(scored)
        # Sorted: 0.3, 0.5, 0.7, 0.9 → median = (0.5 + 0.7) / 2 = 0.6
        assert s.median_score == pytest.approx(0.6)

    def test_median_odd(self):
        scored = [_make_cs(0.9), _make_cs(0.5), _make_cs(0.3)]
        s = compute_summary_statistics(scored)
        assert s.median_score == pytest.approx(0.5)

    def test_defaults_from_candidates(self):
        scored = [_make_cs(0.8), _make_cs(0.6)]
        s = compute_summary_statistics(scored)
        # total_candidates defaults to n
        assert s.total_candidates == 2
        assert s.passed_filters == 2

    def test_chemistry(self):
        scored = [_make_cs(0.8)]
        s = compute_summary_statistics(scored, chemistry="morpholino")
        assert s.chemistry == "morpholino"


# ---------------------------------------------------------------------------
# build_report_data
# ---------------------------------------------------------------------------


class TestBuildReportData:
    def test_basic(self):
        scored = [_make_cs(0.9, rank=1), _make_cs(0.7, rank=2)]
        report = build_report_data(
            scored_candidates=scored,
            total_candidates=100,
            passed_filters=80,
            target_name="BRCA1",
            target_id="NM_007294",
            chemistry="gapmer",
        )
        assert report.summary.target_name == "BRCA1"
        assert report.summary.total_candidates == 100
        assert len(report.candidates) == 2
        assert report.candidates[0]["final_score"] == pytest.approx(0.9)

    def test_with_thermo(self):
        scored = [_make_cs(0.9)]
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            self_structure_penalty=0.3,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
        )
        report = build_report_data(
            scored_candidates=scored,
            thermodynamic_summaries=[thermo],
        )
        assert len(report.thermodynamic_details) == 1
        assert report.thermodynamic_details[0]["tm"] == 55.0

    def test_exclude_thermo(self):
        scored = [_make_cs(0.9)]
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            self_structure_penalty=0.3,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
        )
        report = build_report_data(
            scored_candidates=scored,
            thermodynamic_summaries=[thermo],
            include_thermo_details=False,
        )
        assert len(report.thermodynamic_details) == 0

    def test_with_offtarget_tuple(self):
        scored = [_make_cs(0.9)]
        report = build_report_data(
            scored_candidates=scored,
            off_target_results=[(None, 0.85, "specific")],
        )
        assert len(report.off_target_details) == 1
        assert report.off_target_details[0]["off_target_score"] == pytest.approx(0.85)

    def test_metadata(self):
        scored = [_make_cs(0.9)]
        report = build_report_data(
            scored_candidates=scored,
            metadata={"version": "2.0", "pipeline": "ASODesigner"},
        )
        assert report.metadata["version"] == "2.0"
        assert report.metadata["pipeline"] == "ASODesigner"


# ---------------------------------------------------------------------------
# JSON report generation
# ---------------------------------------------------------------------------


class TestJsonReport:
    def test_generate_and_load(self, tmp_path):
        s = ReportSummary(
            target_name="BRCA1",
            total_candidates=50,
            scored_candidates=40,
            top_score=0.92,
        )
        r = ReportData(
            summary=s,
            candidates=[{"sequence": "ATGC", "final_score": 0.92, "rank": 1}],
            metadata={"version": "1.0"},
        )
        path = generate_json_report(r, tmp_path / "report.json")
        assert path.exists()

        loaded = load_json_report(path)
        assert loaded.summary.target_name == "BRCA1"
        assert loaded.summary.total_candidates == 50
        assert len(loaded.candidates) == 1
        assert loaded.metadata["version"] == "1.0"

    def test_creates_parent_dirs(self, tmp_path):
        r = ReportData()
        path = generate_json_report(r, tmp_path / "subdir" / "report.json")
        assert path.exists()

    def test_load_missing_raises(self):
        with pytest.raises(FileNotFoundError):
            load_json_report("/nonexistent/path/report.json")

    def test_valid_json(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="TP53", total_candidates=10),
            candidates=[{"sequence": "GCAT", "final_score": 0.8, "rank": 1}],
        )
        path = generate_json_report(r, tmp_path / "report.json")
        data = json.loads(path.read_text())
        assert "summary" in data
        assert "candidates" in data
        assert data["summary"]["target_name"] == "TP53"

    def test_round_trip_with_details(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="EGFR", total_candidates=20),
            candidates=[{"sequence": "ATGC", "final_score": 0.9, "rank": 1}],
            thermodynamic_details=[{"tm": 55.0, "binding_class": "strong"}],
            off_target_details=[{"off_target_score": 0.85, "total_hits": 3}],
            variant_details=[{"penalty": 0.2, "has_pathogenic": False}],
            metadata={"version": "1.0"},
        )
        path = generate_json_report(r, tmp_path / "report.json")
        loaded = load_json_report(path)
        assert len(loaded.thermodynamic_details) == 1
        assert len(loaded.off_target_details) == 1
        assert len(loaded.variant_details) == 1
        assert loaded.thermodynamic_details[0]["tm"] == 55.0
        assert loaded.off_target_details[0]["total_hits"] == 3


# ---------------------------------------------------------------------------
# HTML report generation
# ---------------------------------------------------------------------------


class TestHtmlReport:
    def test_basic(self, tmp_path):
        s = ReportSummary(
            target_name="BRCA1",
            target_id="NM_007294",
            total_candidates=100,
            top_score=0.95,
            chemistry="gapmer",
            timestamp="2026-10-03T12:00:00Z",
        )
        candidates = [
            {
                "sequence": "ATGCATGCATGCATGCATGC",
                "start": 100,
                "end": 120,
                "chemistry": "gapmer",
                "final_score": 0.92,
                "rank": 1,
                "component_scores": {
                    "accessibility": 0.9,
                    "thermodynamics": 0.85,
                    "gc_content": 1.0,
                    "off_target": 0.95,
                    "variant": 0.8,
                    "self_structure": 0.7,
                    "positional": 1.0,
                },
            },
        ]
        r = ReportData(summary=s, candidates=candidates)
        path = generate_html_report(r, tmp_path / "report.html")
        assert path.exists()
        content = path.read_text()
        assert "<!DOCTYPE html>" in content
        assert "BRCA1" in content
        assert "summary-grid" in content
        assert "Ranked Candidates" in content

    def test_score_distribution_chart(self, tmp_path):
        candidates = [
            {
                "sequence": f"SEQ{i}",
                "start": 0,
                "end": 20,
                "chemistry": "gapmer",
                "final_score": 0.1 * i,
                "rank": 10 - i,
                "component_scores": {},
            }
            for i in range(1, 11)
        ]
        r = ReportData(
            summary=ReportSummary(total_candidates=10, scored_candidates=10),
            candidates=candidates,
        )
        path = generate_html_report(r, tmp_path / "report.html")
        content = path.read_text()
        assert "Score Distribution" in content
        assert "<svg" in content
        assert "rect" in content

    def test_empty_candidates(self, tmp_path):
        r = ReportData()
        path = generate_html_report(r, tmp_path / "report.html")
        content = path.read_text()
        assert "No candidates to display" in content

    def test_detail_sections(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="TP53"),
            candidates=[
                {"sequence": "ATGC", "final_score": 0.9, "rank": 1, "component_scores": {}}
            ],
            thermodynamic_details=[
                {
                    "tm": 55.0,
                    "tm_score": 0.8,
                    "dg_duplex": -15.0,
                    "dg_self": -2.0,
                    "self_structure_penalty": 0.3,
                    "dg_target": -3.0,
                    "dg_net": -10.0,
                    "binding_class": "strong",
                }
            ],
            off_target_details=[
                {
                    "off_target_score": 0.95,
                    "specificity_class": "specific",
                    "total_hits": 5,
                    "off_target_hits": 2,
                    "perfect_hits": 0,
                }
            ],
            variant_details=[
                {
                    "penalty": 0.2,
                    "overlapping_variants": 1,
                    "seed_variants": 0,
                    "non_seed_variants": 1,
                    "has_pathogenic": False,
                }
            ],
            metadata={"version": "1.0"},
        )
        path = generate_html_report(r, tmp_path / "report.html")
        content = path.read_text()
        assert "Thermodynamic Details" in content
        assert "Off-Target Details" in content
        assert "Variant Details" in content
        assert "Run Metadata" in content
        assert "badge-strong" in content

    def test_creates_parent_dirs(self, tmp_path):
        r = ReportData()
        path = generate_html_report(r, tmp_path / "subdir" / "report.html")
        assert path.exists()

    def test_custom_title(self, tmp_path):
        r = ReportData()
        path = generate_html_report(r, tmp_path / "report.html", title="Custom Title")
        content = path.read_text()
        assert "Custom Title" in content


# ---------------------------------------------------------------------------
# generate_report (combined JSON + HTML)
# ---------------------------------------------------------------------------


class TestGenerateReport:
    def test_both_formats(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="BRCA1", total_candidates=10),
            candidates=[
                {"sequence": "ATGC", "final_score": 0.9, "rank": 1, "component_scores": {}}
            ],
        )
        paths = generate_report(r, tmp_path, basename="brca1")
        assert "json" in paths
        assert "html" in paths
        assert paths["json"].exists()
        assert paths["html"].exists()

    def test_json_only(self, tmp_path):
        r = ReportData()
        paths = generate_report(r, tmp_path, generate_html=False)
        assert "json" in paths
        assert "html" not in paths

    def test_html_only(self, tmp_path):
        r = ReportData()
        paths = generate_report(r, tmp_path, generate_json=False)
        assert "html" in paths
        assert "json" not in paths

    def test_csv_export(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="BRCA1"),
            candidates=[
                {
                    "rank": 1,
                    "sequence": "ATGCATGCATGCATGCATGC",
                    "start": 10,
                    "end": 30,
                    "chemistry": "gapmer",
                    "composite_score": 0.75,
                    "chemistry_bonus": 0.0,
                    "final_score": 0.75,
                    "component_scores": {"gc_content": 0.8, "accessibility": 0.6},
                    "features": {
                        "cpg_count": 2,
                        "g_quadruplex_motifs": 0,
                        "motif_risk": "moderate",
                    },
                },
                {
                    "rank": 2,
                    "sequence": "GCGCGCGCGCGCGCGCGCGC",
                    "start": 40,
                    "end": 60,
                    "chemistry": "gapmer",
                    "composite_score": 0.5,
                    "chemistry_bonus": 0.1,
                    "final_score": 0.6,
                    "component_scores": {"gc_content": 0.9},
                },
            ],
        )
        paths = generate_report(
            r, tmp_path, basename="csvout", generate_html=False, generate_csv=True
        )
        assert "csv" in paths
        assert paths["csv"].exists()

        lines = paths["csv"].read_text().strip().splitlines()
        header = lines[0].split(",")
        # Fixed core columns, then the union of component columns (sorted),
        # then the fixed motif-screen columns.
        assert header[:9] == [
            "rank",
            "sequence",
            "backbone",
            "start",
            "end",
            "chemistry",
            "composite_score",
            "chemistry_bonus",
            "final_score",
        ]
        assert header[9:11] == ["accessibility", "gc_content"]
        assert header[11:] == ["cpg_count", "g_quadruplex_motifs", "motif_risk"]
        row1 = lines[1].split(",")
        assert row1[0] == "1"
        assert row1[1] == "ATGCATGCATGCATGCATGC"
        assert row1[9] == "0.6"  # accessibility present in row 1
        assert row1[11:] == ["2", "0", "moderate"]  # motif features present
        row2 = lines[2].split(",")
        assert row2[9] == ""  # accessibility missing in row 2
        assert row2[10] == "0.9"
        assert row2[11:] == ["", "", ""]  # motif features missing

    def test_csv_empty_candidates(self, tmp_path):
        r = ReportData(summary=ReportSummary(target_name="X"), candidates=[])
        paths = generate_report(
            r, tmp_path, basename="empty", generate_html=False, generate_csv=True
        )
        lines = paths["csv"].read_text().strip().splitlines()
        assert len(lines) == 1  # header only

    def test_fasta_export(self, tmp_path):
        r = ReportData(
            summary=ReportSummary(target_name="BRCA1"),
            candidates=[
                {
                    "rank": 1,
                    "sequence": "ATGCATGCATGCATGCATGC",
                    "start": 10,
                    "end": 30,
                    "chemistry": "gapmer",
                    "final_score": 0.75,
                },
                {
                    "sequence": "GCGCGCGCGCGCGCGCGCGC",
                    "start": 40,
                    "end": 60,
                    "chemistry": "2ome",
                },
            ],
        )
        paths = generate_report(
            r, tmp_path, basename="fa", generate_html=False, generate_fasta=True
        )
        assert "fasta" in paths
        text = paths["fasta"].read_text()
        blocks = text.strip().split(">")
        assert blocks[0] == ""
        first = blocks[1].splitlines()
        assert first[0].startswith("aso_1 rank=1 start=10 end=30 chemistry=gapmer score=0.7500")
        assert first[1] == "ATGCATGCATGCATGCATGC"
        second = blocks[2].splitlines()
        # Unranked candidate: unique index id, no rank/score fields.
        assert second[0].startswith("aso_2 start=40 end=60 chemistry=2ome")
        assert "score=" not in second[0]

    def test_fasta_empty_candidates(self, tmp_path):
        r = ReportData(summary=ReportSummary(target_name="X"), candidates=[])
        paths = generate_report(
            r, tmp_path, basename="nofa", generate_html=False, generate_fasta=True
        )
        assert paths["fasta"].read_text() == ""


# ---------------------------------------------------------------------------
# generate_report_from_ranking
# ---------------------------------------------------------------------------


class TestGenerateReportFromRanking:
    def test_basic(self, tmp_path):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="3UTR"),
        ]
        features = [
            {"accessibility_score": 0.8, "tm_score": 0.9, "binding_class": "strong"},
            {"accessibility_score": 0.3, "tm_score": 0.5, "binding_class": "moderate"},
        ]
        ranking = score_rank_and_filter(cands, features, top_n=2)
        paths = generate_report_from_ranking(
            ranking,
            tmp_path,
            target_name="BRCA1",
            target_id="NM_007294",
            chemistry="gapmer",
            total_candidates=100,
            basename="brca1",
        )
        assert "json" in paths
        assert "html" in paths
        loaded = load_json_report(paths["json"])
        assert loaded.summary.target_name == "BRCA1"
        assert loaded.summary.total_candidates == 100


# ---------------------------------------------------------------------------
# ReportAssembler
# ---------------------------------------------------------------------------


class TestReportAssembler:
    def test_fluent_interface(self):
        scored = [_make_cs(0.9, rank=1), _make_cs(0.7, rank=2)]
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            self_structure_penalty=0.3,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
        )
        assembler = ReportAssembler(
            target_name="BRCA1",
            target_id="NM_007294",
            chemistry="gapmer",
        )
        assembler.set_candidate_counts(total=100, passed=80)
        assembler.set_scored_candidates(scored)
        assembler.set_thermodynamic_summaries([thermo])
        assembler.add_metadata("version", "1.0")

        report = assembler.assemble()
        assert report.summary.target_name == "BRCA1"
        assert report.summary.total_candidates == 100
        assert len(report.candidates) == 2
        assert len(report.thermodynamic_details) == 1
        assert report.metadata["version"] == "1.0"

    def test_exclude_sections(self):
        scored = [_make_cs(0.9)]
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            self_structure_penalty=0.3,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
        )
        assembler = ReportAssembler(target_name="TP53")
        assembler.set_scored_candidates(scored)
        assembler.set_thermodynamic_summaries([thermo])

        report = assembler.assemble(include_thermo=False)
        assert len(report.thermodynamic_details) == 0

    def test_set_metadata(self):
        assembler = ReportAssembler(target_name="EGFR")
        assembler.set_metadata({"key1": "val1", "key2": "val2"})
        report = assembler.assemble()
        assert report.metadata["key1"] == "val1"
        assert report.metadata["key2"] == "val2"


# ---------------------------------------------------------------------------
# assemble_report
# ---------------------------------------------------------------------------


class TestAssembleReport:
    def test_basic(self, tmp_path):
        scored = [_make_cs(0.9, rank=1), _make_cs(0.7, rank=2)]
        paths = assemble_report(
            scored_candidates=scored,
            output_dir=tmp_path,
            target_name="BRCA1",
            target_id="NM_007294",
            chemistry="gapmer",
            total_candidates=100,
            basename="brca1",
        )
        assert "json" in paths
        assert "html" in paths
        loaded = load_json_report(paths["json"])
        assert loaded.summary.target_name == "BRCA1"
        assert loaded.summary.total_candidates == 100

    def test_exclude_details(self, tmp_path):
        scored = [_make_cs(0.9)]
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            self_structure_penalty=0.3,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
        )
        paths = assemble_report(
            scored_candidates=scored,
            output_dir=tmp_path,
            thermodynamic_summaries=[thermo],
            include_thermo=False,
            basename="test",
        )
        loaded = load_json_report(paths["json"])
        assert len(loaded.thermodynamic_details) == 0

    def test_json_only(self, tmp_path):
        scored = [_make_cs(0.9)]
        paths = assemble_report(
            scored_candidates=scored,
            output_dir=tmp_path,
            generate_html=False,
            basename="test",
        )
        assert "json" in paths
        assert "html" not in paths


# ---------------------------------------------------------------------------
# assemble_report_from_config
# ---------------------------------------------------------------------------


class TestAssembleReportFromConfig:
    def test_basic(self, tmp_path):
        config = load_config(chemistry="gapmer")
        config.reporting.output_dir = str(tmp_path)
        scored = [_make_cs(0.9, rank=1), _make_cs(0.7, rank=2)]
        paths = assemble_report_from_config(
            scored_candidates=scored,
            config=config,
            target_name="BRCA1",
            target_id="NM_007294",
            chemistry="gapmer",
            total_candidates=100,
            basename="brca1_config",
        )
        assert "json" in paths
        assert paths["json"].exists()

    def test_html_disabled(self, tmp_path):
        config = load_config(chemistry="gapmer")
        config.reporting.output_dir = str(tmp_path)
        config.reporting.generate_html = False
        scored = [_make_cs(0.9)]
        paths = assemble_report_from_config(
            scored_candidates=scored,
            config=config,
            target_name="TP53",
            basename="tp53",
        )
        assert "json" in paths
        assert "html" not in paths


# ---------------------------------------------------------------------------
# generate_reports_batch
# ---------------------------------------------------------------------------


class TestGenerateReportsBatch:
    def test_basic(self, tmp_path):
        targets = [
            {
                "target_name": "BRCA1",
                "target_id": "NM_007294",
                "scored_candidates": [_make_cs(0.9, rank=1), _make_cs(0.7, rank=2)],
                "total_candidates": 100,
            },
            {
                "target_name": "TP53",
                "target_id": "NM_000546",
                "scored_candidates": [_make_cs(0.85, rank=1)],
                "total_candidates": 80,
            },
        ]
        all_paths = generate_reports_batch(targets, tmp_path)
        assert len(all_paths) == 2
        for paths in all_paths:
            assert "json" in paths
            assert "html" in paths
            assert paths["json"].exists()
            assert paths["html"].exists()

    def test_skip_empty(self, tmp_path):
        targets = [
            {"target_name": "BRCA1", "scored_candidates": [_make_cs(0.9)]},
            {"target_name": "Empty", "scored_candidates": []},
        ]
        all_paths = generate_reports_batch(targets, tmp_path)
        assert len(all_paths) == 1

    def test_per_target_subdir(self, tmp_path):
        targets = [
            {"target_name": "GeneA", "scored_candidates": [_make_cs(0.9)]},
            {"target_name": "GeneB", "scored_candidates": [_make_cs(0.8)]},
        ]
        all_paths = generate_reports_batch(targets, tmp_path)
        assert all_paths[0]["json"].parent.name == "GeneA"
        assert all_paths[1]["json"].parent.name == "GeneB"

    def test_json_only(self, tmp_path):
        targets = [
            {"target_name": "GeneA", "scored_candidates": [_make_cs(0.9)]},
        ]
        all_paths = generate_reports_batch(targets, tmp_path, generate_html=False)
        assert "json" in all_paths[0]
        assert "html" not in all_paths[0]
