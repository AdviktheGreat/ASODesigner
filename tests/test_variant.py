"""Tests for variant/SNP awareness (data structures, VCF, annotation, scoring)."""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from aso_designer.utils.exceptions import VariantError
from aso_designer.variant import (
    ClinVarClient,
    DbSNPClient,
    Variant,
    VariantAnnotation,
    VariantAnnotator,
    VariantFetchManager,
    VariantSet,
    _classify_variant_type,
    _classify_vcf_variant_type,
    _deduplicate_variant_set,
    _extract_clinical_from_info,
    _extract_maf_from_info,
    _normalize_clinical_significance,
    _parse_clinvar_alleles,
    _parse_clinvar_position,
    _parse_info_field,
    _parse_maf_from_summary,
    _parse_maf_value,
    _parse_spdi,
    classify_variant_risk,
    filter_candidates_by_variants,
    load_vcf,
    load_vcf_filtered,
    parse_vcf,
    score_variant_penalty,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_vcf_content() -> str:
    """Standard VCF content for testing."""
    return (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical Significance">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\trs123\tA\tG\t100\tPASS\tAF=0.05\n"
        "1\t200\trs456\tAT\tA\t50\tPASS\tAF=0.001;CLNSIG=pathogenic\n"
        "1\t300\trs789\tC\tT\t80\tPASS\tAF=0.02;CLNSIG=likely pathogenic\n"
        "1\t400\trs012\tG\tC\t90\tPASS\tCLNSIG=benign\n"
        "2\t500\trs345\tA\tT\t70\tPASS\tAF=0.1\n"
    )


@pytest.fixture
def sample_vcf_file(sample_vcf_content, tmp_path):
    """Write a temporary VCF file."""
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(sample_vcf_content)
    return str(vcf_path)


@pytest.fixture
def populated_variant_set() -> VariantSet:
    """A VariantSet with variants at known positions."""
    vs = VariantSet()
    vs.add_variant(Variant("1", 105, "A", "G", maf=0.05, variant_type="snv"))
    vs.add_variant(Variant("1", 115, "C", "T", maf=0.02, variant_type="snv"))
    vs.add_variant(Variant("1", 108, "G", "A", maf=0.01, clinical_significance="pathogenic"))
    vs.add_variant(Variant("1", 500, "T", "C", maf=0.03))
    return vs


# ---------------------------------------------------------------------------
# Variant dataclass
# ---------------------------------------------------------------------------


class TestVariant:
    def test_basic(self):
        v = Variant("1", 100, "A", "G", rs_id="rs123", maf=0.05)
        assert v.chrom == "1"
        assert v.pos == 100
        assert v.ref == "A"
        assert v.alt == "G"
        assert v.rs_id == "rs123"
        assert v.maf == 0.05

    def test_is_pathogenic(self):
        v = Variant("1", 100, "A", "G", clinical_significance="pathogenic")
        assert v.is_pathogenic
        v2 = Variant("1", 100, "A", "G", clinical_significance="likely_pathogenic")
        assert v2.is_pathogenic
        v3 = Variant("1", 100, "A", "G", clinical_significance="benign")
        assert not v3.is_pathogenic

    def test_is_snv(self):
        v = Variant("1", 100, "A", "G", variant_type="snv")
        assert v.is_snv
        v2 = Variant("1", 100, "A", "G", variant_type="indel")
        assert not v2.is_snv

    def test_has_maf(self):
        assert Variant("1", 100, "A", "G", maf=0.05).has_maf
        assert not Variant("1", 100, "A", "G", maf=-1.0).has_maf

    def test_to_dict(self):
        v = Variant("1", 100, "A", "G", maf=0.05, clinical_significance="pathogenic")
        d = v.to_dict()
        assert d["chrom"] == "1"
        assert d["pos"] == 100
        assert d["maf"] == 0.05
        assert d["is_pathogenic"] is True


# ---------------------------------------------------------------------------
# VariantSet
# ---------------------------------------------------------------------------


class TestVariantSet:
    def test_add_and_count(self):
        vs = VariantSet()
        vs.add_variant(Variant("1", 100, "A", "G"))
        vs.add_variant(Variant("1", 200, "C", "T"))
        assert vs.num_variants == 2

    def test_get_at_position(self):
        vs = VariantSet()
        vs.add_variant(Variant("1", 100, "A", "G"))
        result = vs.get_at_position("1", 100)
        assert len(result) == 1
        assert result[0].ref == "A"

    def test_get_at_position_empty(self):
        vs = VariantSet()
        assert vs.get_at_position("1", 999) == []

    def test_get_in_range(self, populated_variant_set):
        result = populated_variant_set.get_in_range("1", 100, 120)
        assert len(result) == 3  # positions 105, 108, 115

    def test_get_in_range_empty_chrom(self, populated_variant_set):
        assert populated_variant_set.get_in_range("99", 1, 1000) == []

    def test_filter_by_maf(self, populated_variant_set):
        filtered = populated_variant_set.filter_by_maf(0.01)
        assert len(filtered) == 4  # all have MAF >= 0.01
        filtered = populated_variant_set.filter_by_maf(0.03)
        assert len(filtered) == 2  # only 0.05 and 0.03

    def test_filter_pathogenic(self, populated_variant_set):
        result = populated_variant_set.filter_pathogenic()
        assert len(result) == 1

    def test_summary(self, populated_variant_set):
        s = populated_variant_set.summary()
        assert s["total_variants"] == 4
        assert s["snvs"] == 2  # only first two have variant_type="snv"
        assert s["pathogenic"] == 1
        assert s["with_maf"] == 4
        assert s["chromosomes"] == 1


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


class TestHelperFunctions:
    def test_parse_spdi(self):
        assert _parse_spdi("NC_000001.1:100:A:G") == ("A", "G")
        assert _parse_spdi("") == ("", "")

    def test_classify_variant_type(self):
        assert _classify_variant_type("snv") == "snv"
        assert _classify_variant_type("snp") == "snv"
        assert _classify_variant_type("deletion") == "indel"
        assert _classify_variant_type("mnp") == "mnv"
        assert _classify_variant_type("") == "unknown"

    def test_parse_maf_from_summary(self):
        assert _parse_maf_from_summary({"global_maf": "A/G:0.05"}) == 0.05
        assert _parse_maf_from_summary({}) == -1.0

    def test_parse_maf_from_summary_multi_value(self):
        """Multi-allelic frequency lists resolve to the maximum value."""
        assert _parse_maf_from_summary({"global_maf": "A/G:0.05,0.03"}) == 0.05
        assert _parse_maf_from_summary({"global_maf": "A/G:0.02,0.11"}) == 0.11

    def test_parse_maf_from_summary_extra_segments(self):
        """Extra colon segments fall back to the last parseable value."""
        assert _parse_maf_from_summary({"global_maf": "A/G:0.05:0.07"}) == 0.07
        assert _parse_maf_from_summary({"global_maf": "A/G"}) == -1.0
        assert _parse_maf_from_summary({"global_maf": "A/G:abc"}) == -1.0

    def test_normalize_clinical_significance(self):
        assert _normalize_clinical_significance("Pathogenic") == "pathogenic"
        assert _normalize_clinical_significance("Likely pathogenic") == "likely_pathogenic"
        assert _normalize_clinical_significance("Benign") == "benign"
        assert _normalize_clinical_significance("Uncertain significance") == "uncertain"

    def test_parse_clinvar_position(self):
        assert _parse_clinvar_position("12345") == 12345
        assert _parse_clinvar_position("12345-12350") == 12345
        assert _parse_clinvar_position("") == 0

    def test_parse_clinvar_alleles(self):
        assert _parse_clinvar_alleles("A/G") == ("A", "G")
        assert _parse_clinvar_alleles("") == ("", "")

    def test_parse_info_field(self):
        result = _parse_info_field("AF=0.05;CLNSIG=pathogenic")
        assert result == {"AF": "0.05", "CLNSIG": "pathogenic"}
        assert _parse_info_field(".") == {}

    def test_extract_maf_from_info(self):
        assert _extract_maf_from_info({"AF": "0.05"}) == 0.05
        assert _extract_maf_from_info({"MAF": "0.1"}) == 0.1
        assert _extract_maf_from_info({}) == -1.0

    def test_extract_maf_from_info_multi_allelic(self):
        """Comma-separated AF lists (multi-allelic sites) resolve to max."""
        assert _extract_maf_from_info({"AF": "0.1,0.2"}) == 0.2
        assert _extract_maf_from_info({"AF": "0.02,0.03,0.4"}) == 0.4

    def test_extract_maf_from_info_unparseable(self):
        assert _extract_maf_from_info({"AF": "not_a_number"}) == -1.0
        assert _extract_maf_from_info({"AF": "0.05,,"}) == 0.05
        assert _extract_maf_from_info({"MAF": ""}) == -1.0

    def test_parse_maf_value(self):
        """The shared value parser handles common formats."""
        assert _parse_maf_value("0.05") == 0.05
        assert _parse_maf_value(" 0.07 ") == 0.07
        assert _parse_maf_value("0.1,0.02") == 0.1
        assert _parse_maf_value("0.1|0.02") == 0.1
        assert _parse_maf_value("") == -1.0
        assert _parse_maf_value("abc") == -1.0
        assert _parse_maf_value("0.05,") == 0.05

    def test_extract_clinical_from_info(self):
        assert _extract_clinical_from_info({"CLNSIG": "pathogenic"}) == "pathogenic"
        assert _extract_clinical_from_info({}) == ""

    def test_classify_vcf_variant_type(self):
        assert _classify_vcf_variant_type("A", "G") == "snv"
        assert _classify_vcf_variant_type("AT", "A") == "indel"
        assert _classify_vcf_variant_type("A", "AT") == "indel"
        assert _classify_vcf_variant_type("AT", "GC") == "mnv"


# ---------------------------------------------------------------------------
# VCF parser
# ---------------------------------------------------------------------------


class TestParseVcf:
    def test_basic_parse(self, sample_vcf_file):
        variants = parse_vcf(sample_vcf_file)
        assert len(variants) == 5
        assert variants[0].chrom == "1"
        assert variants[0].pos == 100
        assert variants[0].ref == "A"
        assert variants[0].alt == "G"
        assert variants[0].maf == 0.05

    def test_clinical_significance(self, sample_vcf_file):
        variants = parse_vcf(sample_vcf_file)
        assert variants[1].clinical_significance == "pathogenic"
        assert variants[2].clinical_significance == "likely_pathogenic"
        assert variants[3].clinical_significance == "benign"

    def test_variant_type(self, sample_vcf_file):
        variants = parse_vcf(sample_vcf_file)
        assert variants[0].variant_type == "snv"
        assert variants[1].variant_type == "indel"

    def test_chromosome_normalization(self, sample_vcf_file):
        """VCF with chr prefix should be normalized."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\trs1\tA\tG\t100\tPASS\tAF=0.05\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(content)
            path = f.name
        variants = parse_vcf(path)
        os.unlink(path)
        assert variants[0].chrom == "1"  # chr prefix stripped

    def test_file_not_found(self):
        with pytest.raises(VariantError, match="VCF file not found"):
            parse_vcf("/nonexistent/path.vcf")

    def test_skip_malformed_lines(self):
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\trs1\tA\tG\t100\tPASS\tAF=0.05\n"
            "bad_line\n"
            "1\t200\trs2\tC\tT\t50\tPASS\tAF=0.01\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(content)
            path = f.name
        variants = parse_vcf(path)
        os.unlink(path)
        assert len(variants) == 2  # malformed line skipped


class TestLoadVcf:
    def test_load_vcf(self, sample_vcf_file):
        vs = load_vcf(sample_vcf_file)
        assert vs.num_variants == 5
        assert isinstance(vs, VariantSet)

    def test_load_vcf_filtered_maf(self, sample_vcf_file):
        vs = load_vcf_filtered(sample_vcf_file, min_maf=0.02)
        assert vs.num_variants == 3  # 0.05, 0.02, 0.1

    def test_load_vcf_filtered_snv(self, sample_vcf_file):
        vs = load_vcf_filtered(sample_vcf_file, variant_types={"snv"})
        assert vs.num_variants == 4  # 4 SNVs, 1 indel

    def test_load_vcf_filtered_pathogenic(self, sample_vcf_file):
        vs = load_vcf_filtered(sample_vcf_file, pathogenic_only=True)
        assert vs.num_variants == 2  # pathogenic + likely_pathogenic


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


class TestDeduplication:
    def test_dedup(self):
        vs = VariantSet()
        vs.add_variant(Variant("1", 100, "A", "G", maf=0.05))
        vs.add_variant(Variant("1", 100, "A", "G", clinical_significance="pathogenic"))
        vs.add_variant(Variant("1", 200, "C", "T"))
        result = _deduplicate_variant_set(vs)
        assert result.num_variants == 2

    def test_dedup_merges_data(self):
        vs = VariantSet()
        vs.add_variant(Variant("1", 100, "A", "G", maf=0.05))
        vs.add_variant(Variant("1", 100, "A", "G", clinical_significance="pathogenic"))
        result = _deduplicate_variant_set(vs)
        v = result.get_at_position("1", 100)[0]
        assert v.maf == 0.05
        assert v.clinical_significance == "pathogenic"


# ---------------------------------------------------------------------------
# API clients (mocked)
# ---------------------------------------------------------------------------


class TestDbSNPClient:
    @patch("urllib.request.urlopen")
    def test_fetch_success(self, mock_urlopen):
        mock_resp = MagicMock()
        mock_resp.read.return_value = b'{"esearchresult": {"idlist": ["123", "456"]}}'
        mock_urlopen.return_value.__enter__ = MagicMock(return_value=mock_resp)
        mock_urlopen.return_value.__exit__ = MagicMock(return_value=False)

        client = DbSNPClient(email="test@test.com")
        # Mock _fetch_summaries to avoid second API call
        with patch.object(
            client,
            "_fetch_summaries",
            return_value=[
                Variant("1", 100, "A", "G", rs_id="rs123", maf=0.05),
                Variant("1", 200, "C", "T", rs_id="rs456", maf=0.01),
            ],
        ):
            variants = client.fetch_variants_by_gene("DMD")
        assert len(variants) == 2

    @patch("urllib.request.urlopen")
    def test_fetch_no_results(self, mock_urlopen):
        mock_resp = MagicMock()
        mock_resp.read.return_value = b'{"esearchresult": {"idlist": []}}'
        mock_urlopen.return_value.__enter__ = MagicMock(return_value=mock_resp)
        mock_urlopen.return_value.__exit__ = MagicMock(return_value=False)

        client = DbSNPClient(email="test@test.com")
        variants = client.fetch_variants_by_gene("UNKNOWNGENE")
        assert len(variants) == 0


class TestClinVarClient:
    @patch("urllib.request.urlopen")
    def test_fetch_success(self, mock_urlopen):
        mock_resp = MagicMock()
        mock_resp.read.return_value = b'{"esearchresult": {"idlist": ["1", "2"]}}'
        mock_urlopen.return_value.__enter__ = MagicMock(return_value=mock_resp)
        mock_urlopen.return_value.__exit__ = MagicMock(return_value=False)

        client = ClinVarClient(email="test@test.com")
        with patch.object(
            client,
            "_fetch_clinvar_summaries",
            return_value=[
                Variant("1", 100, "A", "G", clinical_significance="pathogenic"),
            ],
        ):
            variants = client.fetch_pathogenic_variants("DMD")
        assert len(variants) == 1
        assert variants[0].is_pathogenic


class TestVariantFetchManager:
    def test_init(self):
        mgr = VariantFetchManager(email="test@test.com", vcf_path="/tmp/test.vcf")
        assert mgr.dbsnp is not None
        assert mgr.clinvar is not None
        assert mgr.vcf_path == "/tmp/test.vcf"

    def test_vcf_fallback(self, sample_vcf_file):
        """When API returns nothing, fall back to VCF."""
        mgr = VariantFetchManager(email="test@test.com", vcf_path=sample_vcf_file)
        # Mock both APIs to return empty
        with (
            patch.object(mgr.dbsnp, "fetch_variants_by_gene", return_value=[]),
            patch.object(mgr.clinvar, "fetch_pathogenic_variants", return_value=[]),
        ):
            vs = mgr.fetch_variants("TESTGENE")
        assert vs.num_variants == 5  # From VCF


# ---------------------------------------------------------------------------
# Variant annotation
# ---------------------------------------------------------------------------


class TestVariantAnnotator:
    def test_basic_annotation(self, populated_variant_set):
        annotator = VariantAnnotator(populated_variant_set)
        ann = annotator.annotate(100, 120, chrom="1")
        assert ann.total_variants == 3  # positions 105, 108, 115
        assert ann.seed_variants == 2  # positions 105, 108 (in seed [102,108))
        assert ann.non_seed_variants == 1  # position 115
        assert ann.has_pathogenic is True
        assert ann.penalty > 0

    def test_no_variants(self):
        vs = VariantSet()
        annotator = VariantAnnotator(vs)
        ann = annotator.annotate(100, 120, chrom="1")
        assert ann.total_variants == 0
        assert ann.penalty == 0.0
        assert ann.has_pathogenic is False

    def test_transcript_coordinates(self, populated_variant_set):
        """Without chrom, search all variants."""
        annotator = VariantAnnotator(populated_variant_set)
        ann = annotator.annotate(100, 120)
        assert ann.total_variants == 3

    def test_batch_annotation(self, populated_variant_set):
        annotator = VariantAnnotator(populated_variant_set)
        results = annotator.annotate_batch([(100, 120), (490, 510)], chrom="1")
        assert len(results) == 2
        assert results[0].total_variants == 3
        assert results[1].total_variants == 1

    def test_maf_threshold_filtering(self):
        """Variants below MAF threshold should be filtered."""
        vs = VariantSet()
        vs.add_variant(Variant("1", 105, "A", "G", maf=0.001))
        vs.add_variant(Variant("1", 110, "C", "T", maf=0.05))
        annotator = VariantAnnotator(vs, maf_threshold=0.01)
        ann = annotator.annotate(100, 120, chrom="1")
        assert ann.total_variants == 1  # Only the 0.05 MAF variant

    def test_seed_region_multiplier(self):
        """Seed region multiplier should increase penalty."""
        vs = VariantSet()
        vs.add_variant(Variant("1", 105, "A", "G", maf=0.05))
        annotator_low = VariantAnnotator(vs, seed_region_multiplier=1.0)
        annotator_high = VariantAnnotator(vs, seed_region_multiplier=5.0)
        ann_low = annotator_low.annotate(100, 120, chrom="1")
        ann_high = annotator_high.annotate(100, 120, chrom="1")
        assert ann_high.penalty > ann_low.penalty


# ---------------------------------------------------------------------------
# Variant scoring
# ---------------------------------------------------------------------------


class TestVariantScoring:
    def test_score_no_variants(self):
        ann = VariantAnnotation(
            aso_start=0,
            aso_end=20,
            total_variants=0,
            seed_variants=0,
            non_seed_variants=0,
            pathogenic_variants=0,
            has_pathogenic=False,
            penalty=0.0,
            variants=[],
        )
        assert score_variant_penalty(ann) == 1.0

    def test_score_with_variants(self, populated_variant_set):
        annotator = VariantAnnotator(populated_variant_set)
        ann = annotator.annotate(100, 120, chrom="1")
        score = score_variant_penalty(ann)
        assert 0.0 <= score < 1.0

    def test_classify_low_risk(self):
        ann = VariantAnnotation(
            aso_start=0,
            aso_end=20,
            total_variants=0,
            seed_variants=0,
            non_seed_variants=0,
            pathogenic_variants=0,
            has_pathogenic=False,
            penalty=0.0,
            variants=[],
        )
        assert classify_variant_risk(ann) == "low_risk"

    def test_classify_high_risk_pathogenic(self):
        ann = VariantAnnotation(
            aso_start=0,
            aso_end=20,
            total_variants=1,
            seed_variants=1,
            non_seed_variants=0,
            pathogenic_variants=1,
            has_pathogenic=True,
            penalty=0.3,
            variants=[],
        )
        assert classify_variant_risk(ann) == "high_risk"

    def test_classify_moderate_risk(self):
        ann = VariantAnnotation(
            aso_start=0,
            aso_end=20,
            total_variants=2,
            seed_variants=0,
            non_seed_variants=2,
            pathogenic_variants=0,
            has_pathogenic=False,
            penalty=0.3,
            variants=[],
        )
        assert classify_variant_risk(ann) == "moderate_risk"

    def test_filter_candidates(self, populated_variant_set):
        annotator = VariantAnnotator(populated_variant_set)
        candidates = [(100, 120), (200, 220), (300, 320)]
        filtered = filter_candidates_by_variants(candidates, annotator, max_penalty=2.0, chrom="1")
        assert len(filtered) == 3  # All pass with max_penalty=2.0

    def test_filter_exclude_pathogenic(self, populated_variant_set):
        annotator = VariantAnnotator(populated_variant_set)
        candidates = [(100, 120), (200, 220), (300, 320)]
        filtered = filter_candidates_by_variants(
            candidates,
            annotator,
            max_penalty=2.0,
            exclude_pathogenic=True,
            chrom="1",
        )
        # (100, 120) has pathogenic variant, should be excluded
        assert len(filtered) == 2
