"""Tests for pipeline feature stages (accessibility, thermodynamics, off-target, variant).

These tests verify that the four feature-computation stages write real
per-candidate features that flow into the scoring stage. ViennaRNA and
network access are mocked; the thermodynamics stage runs the fast
nearest-neighbor path which requires no external tools.
"""

from unittest.mock import MagicMock, patch

import pytest

from aso_designer.pipeline import ASOPipeline, PipelineConfig
from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.viennarna import RNAplfoldResult, ViennaRNAError

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

# 126-nt synthetic transcript with a clean ORF structure:
#   positions 0-2     ATG (start codon)
#   positions 3-98    coding sequence, 32 codons (no in-frame stops)
#   positions 99-101  TAA (stop codon)
#   positions 102-125 3' UTR (24 nt, long enough for candidate windows)
# The CDS is 102 nt = 34 codons, above the 30-codon minimum used by
# find_orfs/detect_cds.
TRANSCRIPT_SEQ = "ATG" + "CATGCATGC" * 10 + "CATGCA" + "TAA" + "GCGCGCGCGCAAGCGCGCGCGCAA"

# All 20-mers in the ATGC-repeat region have 50% GC content and no
# homopolymers, so most candidates pass the default filters.


def make_pipeline(sequence: str = TRANSCRIPT_SEQ) -> ASOPipeline:
    """Build an ASOPipeline with the transcript pre-loaded (no fetching)."""
    config = PipelineConfig(target="TEST")
    pipeline = ASOPipeline(config)
    pipeline._transcript = TranscriptRecord(
        accession="NM_TEST.1",
        sequence=sequence,
        gene_name="TESTGENE",
        organism="Homo sapiens",
    )
    return pipeline


def run_through_filter(pipeline: ASOPipeline) -> None:
    """Run generate + filter stages so candidates carry region/flags."""
    pipeline._stage_generate()
    pipeline._stage_filter()


def passing_candidates(pipeline: ASOPipeline) -> list:
    return [c for c in pipeline._candidates if c.passes_filters]


@pytest.fixture
def mock_rnaplfold_result() -> RNAplfoldResult:
    """Accessibility profile: unpaired (0.9) in the first half of the
    transcript, paired (0.1) in the second half."""
    n = len(TRANSCRIPT_SEQ)
    probs = [0.9] * (n // 2) + [0.1] * (n - n // 2)
    return RNAplfoldResult(
        positions=list(range(1, n + 1)),
        unpaired_probs=probs,
    )


# ---------------------------------------------------------------------------
# Generation stage
# ---------------------------------------------------------------------------


class TestStageGenerate:
    def test_generate_matches_direct_generation(self):
        """Batched pipeline generation equals direct generate_candidates."""
        from aso_designer.candidate_generator import generate_candidates
        from aso_designer.quality import extract_metadata

        pipeline = make_pipeline()
        direct = generate_candidates(
            sequence=pipeline._transcript.sequence,
            config=pipeline.aso_config.candidate,
            metadata=extract_metadata(pipeline._transcript),
        )

        pipeline._stage_generate()
        assert pipeline.result.total_candidates == len(direct)
        for got, want in zip(pipeline._candidates, direct, strict=True):
            assert got.sequence == want.sequence
            assert got.start == want.start
            assert got.end == want.end
            assert got.region == want.region
            assert got.features["transcript_length"] == len(pipeline._transcript.sequence)

    def test_batched_generation_equals_direct(self):
        """Flattened batches equal the one-shot generator, any batch size."""
        from aso_designer.candidate_generator import (
            CandidateConfig,
            generate_candidates,
            generate_candidates_batched,
        )

        seq = TRANSCRIPT_SEQ
        cfg = CandidateConfig(length=16, step=3)
        direct = generate_candidates(seq, cfg)
        for batch_size in (1, 5, 10000):
            flattened = [
                c
                for batch in generate_candidates_batched(seq, cfg, batch_size=batch_size)
                for c in batch
            ]
            assert [(c.sequence, c.start, c.end, c.region) for c in flattened] == [
                (c.sequence, c.start, c.end, c.region) for c in direct
            ]


# ---------------------------------------------------------------------------
# Accessibility stage
# ---------------------------------------------------------------------------


class TestAccessibilityStage:
    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_accessibility_scores_written(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        passed = passing_candidates(pipeline)
        assert len(passed) > 0
        for c in passed:
            assert "accessibility_score" in c.features
            score = c.features["accessibility_score"]
            assert 0.0 <= score <= 1.0

    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_accessible_region_scores_higher(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        """A candidate in the unpaired first half should score higher
        than one in the paired second half."""
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        first_half = [
            c.features["accessibility_score"]
            for c in passing_candidates(pipeline)
            if c.start < len(TRANSCRIPT_SEQ) // 2
        ]
        second_half = [
            c.features["accessibility_score"]
            for c in passing_candidates(pipeline)
            if c.start >= len(TRANSCRIPT_SEQ) // 2
        ]
        assert first_half and second_half
        assert min(first_half) > max(second_half)

    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_filtered_candidates_skipped(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        """Candidates that failed filters should not get accessibility scores."""
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        failed = [c for c in pipeline._candidates if not c.passes_filters]
        for c in failed:
            assert "accessibility_score" not in c.features

    @patch("aso_designer.viennarna.verify_viennarna")
    def test_missing_viennarna_raises(self, mock_verify):
        mock_verify.side_effect = ViennaRNAError("RNAfold not found")

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        with pytest.raises(ViennaRNAError):
            pipeline._stage_accessibility()

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_accessibility()


# ---------------------------------------------------------------------------
# Thermodynamics stage
# ---------------------------------------------------------------------------


class TestThermodynamicsStage:
    def test_thermodynamic_features_written(self):
        """NN method requires no ViennaRNA; all features should be populated."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        passed = passing_candidates(pipeline)
        assert len(passed) > 0
        feature_keys = (
            "tm",
            "tm_score",
            "dg_duplex",
            "dg_self",
            "dg_target",
            "dg_net",
            "binding_class",
            "self_structure_penalty",
        )
        for c in passed:
            for key in feature_keys:
                assert key in c.features, f"missing {key} for candidate at {c.start}"
            assert c.features["tm"] > 0.0
            assert 0.0 <= c.features["tm_score"] <= 1.0
            assert c.features["binding_class"] in ("strong", "moderate", "weak")
            assert 0.0 <= c.features["self_structure_penalty"]

    def test_summaries_collected_for_reporting(self):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        passed = passing_candidates(pipeline)
        assert len(pipeline._thermo_summaries) == len(passed)

    def test_filtered_candidates_skipped(self):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        for c in pipeline._candidates:
            if not c.passes_filters:
                assert "tm" not in c.features

    def test_config_parameters_applied(self, tmp_path):
        """Custom Tm bounds from a user config should flow into the summaries."""
        config_file = tmp_path / "custom.yaml"
        config_file.write_text(
            "thermodynamics:\n  tm_min: 30.0\n  tm_max: 60.0\n  tm_optimal: 45.0\n"
        )
        config = PipelineConfig(target="TEST", config_path=str(config_file))
        pipeline = ASOPipeline(config)
        pipeline._transcript = TranscriptRecord(
            accession="NM_TEST.1",
            sequence=TRANSCRIPT_SEQ,
            gene_name="TESTGENE",
            organism="Homo sapiens",
        )
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        assert len(pipeline._thermo_summaries) == len(passing_candidates(pipeline))
        for summary in pipeline._thermo_summaries:
            assert 0.0 <= summary.tm_score <= 1.0

    def test_invalid_method_rejected(self, tmp_path):
        config_file = tmp_path / "bad.yaml"
        config_file.write_text("thermodynamics:\n  method: bogus\n")
        pipeline = ASOPipeline(PipelineConfig(target="TEST", config_path=str(config_file)))
        with pytest.raises(Exception, match="Config validation failed"):
            _ = pipeline.aso_config

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_thermodynamics()


# ---------------------------------------------------------------------------
# Off-target stage
# ---------------------------------------------------------------------------


class TestOffTargetStage:
    @staticmethod
    def _write_fasta(path, entries):
        with open(path, "w") as f:
            for header, seq in entries:
                f.write(f">{header}\n{seq}\n")

    def test_skips_without_transcriptome(self):
        """No transcriptome configured -> warn and return, no crash."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_off_target()

        for c in passing_candidates(pipeline):
            assert "off_target_score" not in c.features

    def test_scores_with_transcriptome(self, tmp_path):
        from aso_designer.utils.sequence import reverse_complement

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        passing = passing_candidates(pipeline)
        assert passing
        probe = passing[0].sequence

        fasta = tmp_path / "transcriptome.fna"
        self._write_fasta(
            fasta,
            [
                ("NM_TEST.1 TESTGENE", TRANSCRIPT_SEQ),
                ("NM_DECOY.1 DECOYGENE", reverse_complement(probe) * 3),
            ],
        )
        pipeline.config.transcriptome_path = str(fasta)
        pipeline._stage_off_target()

        scored = [c for c in passing if "off_target_score" in c.features]
        assert len(scored) == len(passing)
        for c in scored:
            assert 0.0 <= c.features["off_target_score"] <= 1.0
            assert c.features["off_target_class"] in (
                "specific",
                "moderate",
                "non-specific",
            )

    def test_off_target_hit_reduces_score(self, tmp_path):
        """A candidate with an exact off-target match should score lower
        than candidates without one."""
        from aso_designer.utils.sequence import reverse_complement

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        passing = passing_candidates(pipeline)
        probe = passing[0].sequence

        fasta = tmp_path / "transcriptome.fna"
        self._write_fasta(
            fasta,
            [
                ("NM_TEST.1 TESTGENE", TRANSCRIPT_SEQ),
                ("NM_DECOY.1 DECOYGENE", reverse_complement(probe) * 3),
            ],
        )
        pipeline.config.transcriptome_path = str(fasta)
        pipeline._stage_off_target()

        probe_score = next(c.features["off_target_score"] for c in passing if c.sequence == probe)
        other_scores = [c.features["off_target_score"] for c in passing if c.sequence != probe]
        assert other_scores
        assert probe_score < max(other_scores)

    def test_results_collected_for_reporting(self, tmp_path):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        fasta = tmp_path / "transcriptome.fna"
        self._write_fasta(fasta, [("NM_TEST.1 TESTGENE", TRANSCRIPT_SEQ)])
        pipeline.config.transcriptome_path = str(fasta)
        pipeline._stage_off_target()

        assert len(pipeline._off_target_results) == len(passing_candidates(pipeline))

    def test_empty_fasta_skips(self, tmp_path):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        fasta = tmp_path / "empty.fna"
        fasta.write_text("")
        pipeline.config.transcriptome_path = str(fasta)
        pipeline._stage_off_target()

        for c in passing_candidates(pipeline):
            assert "off_target_score" not in c.features

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_off_target()


# ---------------------------------------------------------------------------
# Variant stage
# ---------------------------------------------------------------------------

MINIMAL_VCF = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chrT\t26\trs1\tC\tT\t.\t.\tAF=0.05;CLNSIG=Benign\n"
    "chrT\t52\trs2\tA\tG\t.\t.\tAF=0.30;CLNSIG=Pathogenic\n"
)


class TestVariantStage:
    @staticmethod
    def _write_vcf(path, content=MINIMAL_VCF):
        path.write_text(content)

    def test_skips_without_variant_source(self):
        """No VCF and no gene name -> warn and return, no crash."""
        pipeline = make_pipeline()
        pipeline._transcript.gene_name = ""
        run_through_filter(pipeline)
        pipeline._stage_variant()

        for c in passing_candidates(pipeline):
            assert "variant_penalty" not in c.features

    def test_annotates_with_vcf(self, tmp_path):
        vcf = tmp_path / "variants.vcf"
        self._write_vcf(vcf)

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline.config.vcf_path = str(vcf)
        pipeline._stage_variant()

        passed = passing_candidates(pipeline)
        annotated = [c for c in passed if "variant_penalty" in c.features]
        assert len(annotated) == len(passed)
        for c in annotated:
            assert 0.0 <= c.features["variant_penalty"]
            assert c.features["variant_risk"] in ("low_risk", "moderate_risk", "high_risk")

    def test_variant_overlap_increases_penalty(self, tmp_path):
        """Candidates overlapping a variant should carry a higher penalty
        than candidates in clean regions."""
        vcf = tmp_path / "variants.vcf"
        self._write_vcf(vcf)

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline.config.vcf_path = str(vcf)
        pipeline._stage_variant()

        # Variants at transcript positions 26 and 52 (1-based); a candidate
        # [start, end) covers 1-based [start+1, end].
        def overlaps(c, pos):
            return c.start + 1 <= pos <= c.end

        overlapping = [
            c.features["variant_penalty"]
            for c in passing_candidates(pipeline)
            if overlaps(c, 26) or overlaps(c, 52)
        ]
        clean = [
            c.features["variant_penalty"]
            for c in passing_candidates(pipeline)
            if not overlaps(c, 26) and not overlaps(c, 52)
        ]
        assert overlapping and clean
        assert max(overlapping) > 0.0
        assert all(p == 0.0 for p in clean)

    def test_annotations_collected_for_reporting(self, tmp_path):
        vcf = tmp_path / "variants.vcf"
        self._write_vcf(vcf)

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline.config.vcf_path = str(vcf)
        pipeline._stage_variant()

        assert len(pipeline._variant_annotations) == len(passing_candidates(pipeline))

    def test_maf_below_threshold_ignored(self, tmp_path):
        """A variant with MAF below the threshold should not be annotated."""
        vcf_content = MINIMAL_VCF.replace("AF=0.05", "AF=0.001").replace("AF=0.30", "AF=0.001")
        vcf = tmp_path / "variants.vcf"
        self._write_vcf(vcf, vcf_content)

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline.config.vcf_path = str(vcf)
        pipeline._stage_variant()

        # Default maf_threshold=0.01 filters both variants out
        for c in passing_candidates(pipeline):
            assert c.features.get("variant_penalty", 0.0) == 0.0
            assert c.features.get("variant_risk", "") == "low_risk"

    def test_invalid_vcf_skips_gracefully(self, tmp_path):
        vcf = tmp_path / "bad.vcf"
        vcf.write_text("not a vcf at all\n")

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline.config.vcf_path = str(vcf)
        # Should not raise — bad VCF data is a skip condition
        pipeline._stage_variant()

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_variant()


# ---------------------------------------------------------------------------
# Score stage: features must flow into scores
# ---------------------------------------------------------------------------


class TestScoreStageIntegration:
    def test_features_flow_into_scores(self):
        """Scores computed with real thermodynamic features must differ
        from scores with neutral defaults."""
        from aso_designer.scoring import create_score_engine, score_candidates_batch

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        passed = passing_candidates(pipeline)
        assert passed

        w = pipeline.aso_config.scoring.weights
        engine = create_score_engine(
            weights={
                "accessibility": w.accessibility,
                "thermodynamics": w.thermodynamics,
                "gc_content": w.gc_content,
                "off_target": w.off_target,
                "variant": w.variant,
                "self_structure": w.self_structure,
                "positional": w.positional,
            },
            register_all=True,
        )

        # Score with the computed features
        scored_real = score_candidates_batch(passed, [c.features for c in passed], engine=engine)
        # Score with empty features (neutral defaults)
        scored_neutral = score_candidates_batch(passed, [{} for _ in passed], engine=engine)

        real_thermo = [s.component_scores["thermodynamics"] for s in scored_real]
        neutral_thermo = [s.component_scores["thermodynamics"] for s in scored_neutral]
        # Real Tm scores vary across candidates; neutral is always 0.5
        assert any(t != 0.5 for t in real_thermo)
        assert all(t == 0.5 for t in neutral_thermo)
        assert real_thermo != neutral_thermo

    def test_regions_annotated_from_orf_detection(self):
        """Candidates should be tagged 5UTR/CDS/3UTR, not all 'unknown'."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)

        regions = {c.region for c in pipeline._candidates}
        assert "unknown" not in regions
        assert "CDS" in regions
        assert "3UTR" in regions

    def test_transcript_length_in_features(self):
        pipeline = make_pipeline()
        run_through_filter(pipeline)

        for c in pipeline._candidates:
            assert c.features["transcript_length"] == len(TRANSCRIPT_SEQ)

    def test_full_stage_sequence_ranks_consistently(self):
        """thermodynamics → score → verify ranking is descending."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()
        pipeline._stage_score()

        assert pipeline._scored
        scores = [s.final_score for s in pipeline._scored]
        assert scores == sorted(scores, reverse=True)

    def test_report_metadata_carries_package_version(self, tmp_path):
        import json

        import aso_designer

        pipeline = make_pipeline()
        pipeline.config.output_dir = str(tmp_path / "results")
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()
        pipeline._stage_score()
        pipeline._stage_report()

        report_json = tmp_path / "results" / "TEST_report.json"
        assert report_json.exists()
        data = json.loads(report_json.read_text())
        assert data.get("metadata", {}).get("version") == aso_designer.__version__


# ---------------------------------------------------------------------------
# Pipeline error handling contract
# ---------------------------------------------------------------------------


class TestPipelineErrorContract:
    @staticmethod
    def _offline(pipeline: ASOPipeline) -> None:
        """Bypass fetching (transcript pre-loaded) and network stages."""
        pipeline._stage_fetch = lambda: None  # type: ignore[method-assign]
        pipeline.config.run_accessibility = False
        pipeline.config.run_off_target = False
        pipeline.config.run_variant = False

    def test_thermodynamics_failure_is_non_fatal(self):
        """A thermodynamics failure should not abort the pipeline; scoring
        and reporting continue with neutral features."""
        from aso_designer.utils.exceptions import ThermodynamicsError

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        self._offline(pipeline)

        def boom():
            raise ThermodynamicsError("NN table unavailable")

        pipeline._stage_thermodynamics = boom  # type: ignore[method-assign]
        result = pipeline.run()

        assert any("thermodynamics" in e for e in result.errors)
        assert "score" in result.stages_run
        assert "report" in result.stages_run
        assert result.ranked_candidates > 0

    def test_accessibility_failure_is_non_fatal(self):
        from aso_designer.utils.exceptions import ViennaRNAError

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        self._offline(pipeline)
        # Re-enable accessibility: the patched stage replaces the real
        # ViennaRNA call entirely.
        pipeline.config.run_accessibility = True

        def boom():
            raise ViennaRNAError("RNAfold not found")

        pipeline._stage_accessibility = boom  # type: ignore[method-assign]
        result = pipeline.run()

        assert any("accessibility" in e for e in result.errors)
        assert result.ranked_candidates > 0

    def test_programming_error_propagates(self):
        """Bugs (TypeError etc.) must not be swallowed by stage handling."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        self._offline(pipeline)

        def buggy():
            raise TypeError("programming error, not a data problem")

        pipeline._stage_thermodynamics = buggy  # type: ignore[method-assign]
        with pytest.raises(TypeError, match="programming error"):
            pipeline.run()

    def test_missing_transcript_aborts(self):
        """A RuntimeError in a fatal stage stops the pipeline (recorded,
        not raised) and scoring is skipped."""
        pipeline = make_pipeline()
        pipeline._transcript = None
        self._offline(pipeline)
        result = pipeline.run()

        assert any("generate" in e for e in result.errors)
        assert "score" not in result.stages_run
