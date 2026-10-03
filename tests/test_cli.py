"""Tests for the CLI module: argument parsing, commands, and pipeline."""

import json
import os
import subprocess

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.candidate_serializer import save_candidates
from aso_designer.cli import (
    build_parser,
    format_score_table,
    format_summary,
    is_accession,
    main,
    validate_input_file,
    validate_output_dir,
    validate_target,
)
from aso_designer.pipeline import PipelineConfig, PipelineResult, run_pipeline

# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


class TestBuildParser:
    def test_parser_has_subcommands(self):
        parser = build_parser()
        # Check subcommands exist
        actions = {a.dest: a for a in parser._actions}
        assert "command" in actions

    def test_no_command_returns_0(self):
        assert main([]) == 0

    def test_version_flag(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_design_requires_target(self):
        with pytest.raises(SystemExit):
            main(["design"])

    def test_score_requires_candidates(self):
        with pytest.raises(SystemExit):
            main(["score"])

    def test_report_requires_results(self):
        with pytest.raises(SystemExit):
            main(["report"])


# ---------------------------------------------------------------------------
# Argument validation helpers
# ---------------------------------------------------------------------------


class TestValidateTarget:
    def test_gene_symbol(self):
        assert validate_target("brca1") == "BRCA1"

    def test_accession(self):
        assert validate_target("nm_007294") == "NM_007294"

    def test_strips_whitespace(self):
        assert validate_target("  BRCA1  ") == "BRCA1"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_target("")

    def test_invalid_chars_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            validate_target("BRCA1; rm -rf /")


class TestIsAccession:
    def test_nm_accession(self):
        assert is_accession("NM_007294") is True

    def test_nr_accession(self):
        assert is_accession("NR_123456") is True

    def test_xm_with_version(self):
        assert is_accession("XM_456789.1") is True

    def test_gene_symbol(self):
        assert is_accession("BRCA1") is False

    def test_empty(self):
        assert is_accession("") is False


class TestValidateInputFile:
    def test_valid_file(self, tmp_path):
        f = tmp_path / "test.json"
        f.write_text("{}")
        result = validate_input_file(str(f))
        assert result == str(f.resolve())

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            validate_input_file("/nonexistent/file.json")

    def test_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            validate_input_file(str(tmp_path), "test dir")


class TestValidateOutputDir:
    def test_creates_dir(self, tmp_path):
        path = validate_output_dir(str(tmp_path / "newdir"))
        import os

        assert os.path.isdir(path)

    def test_existing_dir(self, tmp_path):
        path = validate_output_dir(str(tmp_path))
        assert os.path.isdir(path)


# ---------------------------------------------------------------------------
# CLI utility functions
# ---------------------------------------------------------------------------


class TestFormatScoreTable:
    def test_with_candidates(self):
        from aso_designer.scoring import CandidateScore, CompositeScore

        cs = CandidateScore(
            "ATGCATGCATGCATGCATGC",
            0,
            20,
            "gapmer",
            CompositeScore([], 0.9, 0, 0.9),
            0.9,
            rank=1,
        )
        table = format_score_table([cs])
        assert "Rank" in table
        assert "ATGC" in table

    def test_empty(self):
        assert "No candidates" in format_score_table([])


class TestFormatSummary:
    def test_basic(self):
        stats = {"target_name": "BRCA1", "total_candidates": 100, "top_score": 0.95}
        text = format_summary(stats)
        assert "BRCA1" in text
        assert "100" in text

    def test_floats(self):
        stats = {"score": 0.95}
        text = format_summary(stats)
        assert "0.9500" in text

    def test_empty(self):
        assert format_summary({}) == ""


# ---------------------------------------------------------------------------
# Pipeline tests
# ---------------------------------------------------------------------------


def _write_test_fasta(path: str, seq_id: str = "NM_TEST_100", length: int = 200) -> str:
    """Write a test FASTA file and return its path."""
    seq = "ATGCATGCATGCATGCATGC" * (length // 20 + 1)
    seq = seq[:length]
    with open(path, "w") as f:
        f.write(f">{seq_id} Test gene\n{seq}\n")
    return path


class TestPipelineConfig:
    def test_defaults(self):
        cfg = PipelineConfig()
        assert cfg.target == ""
        assert cfg.chemistry == "gapmer"
        assert cfg.run_accessibility is True

    def test_from_args(self):
        import argparse

        args = argparse.Namespace(
            target="BRCA1",
            chemistry="mixmer",
            config=None,
            local_fasta=None,
            output_dir="out",
            top_n=5,
            min_score=0.3,
            no_offtarget=True,
            no_variant=False,
            no_html=False,
        )
        cfg = PipelineConfig.from_args(args)
        assert cfg.target == "BRCA1"
        assert cfg.chemistry == "mixmer"
        assert cfg.top_n == 5
        assert cfg.run_off_target is False
        assert cfg.run_variant is True


class TestPipelineResult:
    def test_defaults(self):
        r = PipelineResult()
        assert r.success is False
        assert r.ranked_candidates == 0

    def test_success(self):
        r = PipelineResult(ranked_candidates=5, top_score=0.9)
        assert r.success is True

    def test_to_dict(self):
        r = PipelineResult(target="BRCA1", total_candidates=100, top_score=0.95)
        d = r.to_dict()
        assert d["target"] == "BRCA1"
        assert d["total_candidates"] == 100


class TestRunPipeline:
    def test_with_local_fasta(self, tmp_path):
        fasta_path = _write_test_fasta(str(tmp_path / "test.fasta"))
        output_dir = str(tmp_path / "output")
        result = run_pipeline(
            target="NM_TEST_100",
            local_fasta=fasta_path,
            output_dir=output_dir,
            run_accessibility=False,
            run_thermodynamics=False,
            run_off_target=False,
            run_variant=False,
            top_n=3,
        )
        assert result.success
        assert result.transcript_length > 0
        assert result.total_candidates > 0
        assert result.ranked_candidates == 3
        assert "json" in result.report_paths
        assert os.path.exists(result.report_paths["json"])
        assert os.path.exists(result.report_paths["html"])

    def test_no_transcript_fails(self):
        result = run_pipeline(
            target="NONEXISTENT",
            run_accessibility=False,
            run_thermodynamics=False,
            run_off_target=False,
            run_variant=False,
        )
        assert not result.success
        assert len(result.errors) > 0


# ---------------------------------------------------------------------------
# CLI command tests (subprocess)
# ---------------------------------------------------------------------------


class TestCliDesign:
    def test_design_success(self, tmp_path):
        fasta_path = _write_test_fasta(str(tmp_path / "test.fasta"))
        output_dir = str(tmp_path / "output")
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "design",
                "--target",
                "NM_TEST_100",
                "--local-fasta",
                fasta_path,
                "--output-dir",
                output_dir,
                "--top-n",
                "3",
                "--no-offtarget",
                "--no-variant",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert os.path.exists(os.path.join(output_dir, "NM_TEST_100_report.json"))

    def test_design_invalid_target(self):
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "design",
                "--target",
                "INVALID; rm -rf /",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1

    def test_design_missing_fasta(self):
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "design",
                "--target",
                "BRCA1",
                "--local-fasta",
                "/nonexistent/file.fasta",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1


class TestCliScore:
    def test_score_to_stdout(self, tmp_path):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="3UTR"),
        ]
        cand_path = str(tmp_path / "cands.json")
        save_candidates(cands, cand_path)

        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "score",
                "--candidates",
                cand_path,
                "--top-n",
                "2",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "ranked_candidates" in result.stdout

    def test_score_to_file(self, tmp_path):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
        ]
        cand_path = str(tmp_path / "cands.json")
        save_candidates(cands, cand_path)

        output_path = str(tmp_path / "output" / "scored.json")
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "score",
                "--candidates",
                cand_path,
                "--output",
                output_path,
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert os.path.exists(output_path)

    def test_score_missing_file(self):
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "score",
                "--candidates",
                "/nonexistent/file.json",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1


class TestCliReport:
    def test_report_from_scored(self, tmp_path):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
        ]
        cand_path = str(tmp_path / "cands.json")
        save_candidates(cands, cand_path)

        scored_path = str(tmp_path / "output" / "scored.json")
        subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "score",
                "--candidates",
                cand_path,
                "--output",
                scored_path,
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        report_dir = str(tmp_path / "reports")
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "report",
                "--results",
                scored_path,
                "--output-dir",
                report_dir,
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert os.path.exists(os.path.join(report_dir, "aso_report.json"))
        assert os.path.exists(os.path.join(report_dir, "aso_report.html"))

    def test_report_no_html(self, tmp_path):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
        ]
        cand_path = str(tmp_path / "cands.json")
        save_candidates(cands, cand_path)

        scored_path = str(tmp_path / "output" / "scored.json")
        subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "score",
                "--candidates",
                cand_path,
                "--output",
                scored_path,
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        report_dir = str(tmp_path / "reports")
        result = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "report",
                "--results",
                scored_path,
                "--output-dir",
                report_dir,
                "--no-html",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert os.path.exists(os.path.join(report_dir, "aso_report.json"))
        assert not os.path.exists(os.path.join(report_dir, "aso_report.html"))


class TestCliInfo:
    def test_info_basic(self):
        result = subprocess.run(
            ["python3", "-m", "aso_designer", "info", "--quiet"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Version" in result.stdout
        assert "Candidate" in result.stdout

    def test_info_weights(self):
        result = subprocess.run(
            ["python3", "-m", "aso_designer", "info", "--show-weights", "--quiet"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Scoring Weights" in result.stdout

    def test_info_presets(self):
        result = subprocess.run(
            ["python3", "-m", "aso_designer", "info", "--show-presets", "--quiet"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Chemistry Presets" in result.stdout


class TestCliHelp:
    def test_no_command(self):
        result = subprocess.run(
            ["python3", "-m", "aso_designer"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "design" in result.stdout
        assert "score" in result.stdout

    def test_version(self):
        result = subprocess.run(
            ["python3", "-m", "aso_designer", "--version"],
            capture_output=True,
            text=True,
        )
        assert "ASODesigner" in result.stdout


# ---------------------------------------------------------------------------
# End-to-end integration test
# ---------------------------------------------------------------------------


class TestEndToEnd:
    def test_full_pipeline_design_to_report(self, tmp_path):
        """Full pipeline: design → score → report → verify."""
        # 1. Design
        fasta_path = _write_test_fasta(
            str(tmp_path / "test.fasta"), seq_id="NM_E2E_001", length=300
        )
        design_dir = str(tmp_path / "design_output")
        r1 = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "design",
                "--target",
                "NM_E2E_001",
                "--local-fasta",
                fasta_path,
                "--output-dir",
                design_dir,
                "--top-n",
                "5",
                "--no-offtarget",
                "--no-variant",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert r1.returncode == 0
        design_json = os.path.join(design_dir, "NM_E2E_001_report.json")
        assert os.path.exists(design_json)

        # 2. Load and verify the report
        with open(design_json) as f:
            report_data = json.loads(f.read())
        assert report_data["summary"]["target_name"] == "NM_E2E_001"
        assert report_data["summary"]["total_candidates"] > 0
        assert len(report_data["candidates"]) == 5

        # 3. Report from the design output
        report_dir = str(tmp_path / "final_reports")
        r3 = subprocess.run(
            [
                "python3",
                "-m",
                "aso_designer",
                "report",
                "--results",
                design_json,
                "--output-dir",
                report_dir,
                "--basename",
                "final_e2e",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )
        assert r3.returncode == 0
        assert os.path.exists(os.path.join(report_dir, "final_e2e.json"))
        assert os.path.exists(os.path.join(report_dir, "final_e2e.html"))
