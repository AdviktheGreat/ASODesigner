"""Tests for the quality module: ORF detection, validation, metadata, and batch reports."""

from aso_designer.quality import (
    TranscriptSet,
    classify_position,
    detect_cds,
    extract_metadata,
    find_orfs,
    generate_quality_report,
    validate_transcript,
)
from aso_designer.transcript_fetcher import TranscriptRecord

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_record(
    seq: str,
    accession: str = "NM_TEST",
    gene_name: str = "TESTGENE",
) -> TranscriptRecord:
    """Build a TranscriptRecord for testing."""
    return TranscriptRecord(
        accession=accession,
        gene_name=gene_name,
        sequence=seq,
        organism="Homo sapiens",
    )


# A sequence with a clear ORF: ATG...TAA
_CDS_SEQ = (
    "AAA" * 20  # 60 nt 5' UTR
    + "ATG"
    + "GCT" * 30
    + "TAA"  # 96 nt CDS
    + "TTT" * 20  # 60 nt 3' UTR
)


# ---------------------------------------------------------------------------
# find_orfs
# ---------------------------------------------------------------------------


class TestFindOrfs:
    def test_finds_orf(self):
        orfs = find_orfs("ATGGCTGCTGCTTAA", min_length=5)
        assert len(orfs) >= 1
        assert orfs[0][0] == 0
        assert orfs[0][1] == 15

    def test_no_orf_short(self):
        orfs = find_orfs("ATGTAA", min_length=10)
        assert orfs == []

    def test_no_start_codon(self):
        orfs = find_orfs("GCTGCTGCTGCTGCTGCT", min_length=5)
        assert orfs == []

    def test_multiple_orfs(self):
        seq = "ATGGCTGCTGCTGCTTAATGCGCATGCGCTAA"
        orfs = find_orfs(seq, min_length=5)
        assert len(orfs) >= 1

    def test_all_three_frames(self):
        seq = "ATG" + "GCT" * 10 + "TAA" + "ATG" + "GCT" * 10 + "TAA"
        orfs = find_orfs(seq, min_length=10)
        assert len(orfs) >= 1


# ---------------------------------------------------------------------------
# detect_cds
# ---------------------------------------------------------------------------


class TestDetectCds:
    def test_detects_longest_orf(self):
        cds = detect_cds(_CDS_SEQ)
        assert cds is not None
        start, end = cds
        assert end - start >= 90

    def test_no_cds(self):
        cds = detect_cds("GCTGCTGCTGCTGCT")
        assert cds is None

    def test_short_orf_filtered(self):
        # detect_cds uses min_length=30 internally
        cds = detect_cds("ATGTAA")
        assert cds is None


# ---------------------------------------------------------------------------
# validate_transcript
# ---------------------------------------------------------------------------


class TestValidateTranscript:
    def test_valid_transcript(self):
        record = _make_record(_CDS_SEQ)
        report = validate_transcript(record, min_length=100)
        assert report.is_valid
        assert report.length == len(_CDS_SEQ)
        assert report.has_cds
        assert report.orf_count >= 1

    def test_too_short(self):
        record = _make_record("ATGCATGC")
        report = validate_transcript(record, min_length=200)
        assert not report.is_valid
        assert any("too short" in e.lower() for e in report.errors)

    def test_ambiguous_bases(self):
        record = _make_record("ATGC" * 50 + "NNNN" * 20)
        report = validate_transcript(record, min_length=100, max_ambiguous_fraction=0.05)
        assert not report.is_valid
        assert any("ambiguous" in e.lower() for e in report.errors)

    def test_require_cds(self):
        record = _make_record("GCT" * 100)
        report = validate_transcript(record, min_length=100, require_cds=True)
        assert not report.is_valid
        assert any("CDS" in e for e in report.errors)

    def test_no_cds_warning(self):
        record = _make_record("GCT" * 100)
        report = validate_transcript(record, min_length=100, require_cds=False)
        assert report.is_valid
        assert any("No CDS" in w for w in report.warnings)

    def test_quality_report_to_dict(self):
        record = _make_record(_CDS_SEQ)
        report = validate_transcript(record, min_length=100)
        d = report.to_dict()
        assert d["accession"] == "NM_TEST"
        assert d["is_valid"] is True
        assert d["length"] == len(_CDS_SEQ)
        assert "warnings" in d
        assert "errors" in d


# ---------------------------------------------------------------------------
# extract_metadata
# ---------------------------------------------------------------------------


class TestExtractMetadata:
    def test_with_cds(self):
        record = _make_record(_CDS_SEQ)
        meta = extract_metadata(record)
        assert meta.accession == "NM_TEST"
        assert meta.cds_start is not None
        assert meta.cds_end is not None
        assert meta.cds_length > 0
        assert meta.five_utr_length > 0
        assert meta.three_utr_length > 0
        assert len(meta.region_map) >= 2

    def test_without_cds(self):
        record = _make_record("GCT" * 100)
        meta = extract_metadata(record)
        assert meta.cds_start is None
        assert meta.cds_length == 0
        assert len(meta.region_map) == 1
        assert meta.region_map[0][2] == "noncoding"

    def test_metadata_to_dict(self):
        record = _make_record(_CDS_SEQ)
        meta = extract_metadata(record)
        d = meta.to_dict()
        assert d["accession"] == "NM_TEST"
        assert d["gene_name"] == "TESTGENE"
        assert "cds_start" in d
        assert "region_map" in d


# ---------------------------------------------------------------------------
# classify_position
# ---------------------------------------------------------------------------


class TestClassifyPosition:
    def test_classify_cds(self):
        record = _make_record(_CDS_SEQ)
        meta = extract_metadata(record)
        cds_pos = meta.cds_start + 10
        assert classify_position(meta, cds_pos) == "CDS"

    def test_classify_5utr(self):
        record = _make_record(_CDS_SEQ)
        meta = extract_metadata(record)
        if meta.five_utr_length > 0:
            assert classify_position(meta, 5) == "5UTR"

    def test_classify_3utr(self):
        record = _make_record(_CDS_SEQ)
        meta = extract_metadata(record)
        if meta.three_utr_length > 0:
            assert classify_position(meta, meta.cds_end + 10) == "3UTR"

    def test_classify_noncoding(self):
        record = _make_record("GCT" * 100)
        meta = extract_metadata(record)
        assert classify_position(meta, 50) == "noncoding"

    def test_classify_unknown(self):
        record = _make_record("GCT" * 100)
        meta = extract_metadata(record)
        assert classify_position(meta, 999) == "unknown"


# ---------------------------------------------------------------------------
# TranscriptSet
# ---------------------------------------------------------------------------


class TestTranscriptSet:
    def test_add_transcript(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("ATGC" * 100, "NM_001"))
        assert len(ts.transcripts) == 1

    def test_longest_transcript(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("ATGC" * 50, "NM_SHORT"))
        ts.add_transcript(_make_record("ATGC" * 100, "NM_LONG"))
        longest = ts.longest_transcript()
        assert longest.accession == "NM_LONG"

    def test_longest_transcript_empty(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        assert ts.longest_transcript() is None

    def test_longest_cds_transcript(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("GCT" * 100, "NM_NOCDS"))
        ts.add_transcript(_make_record(_CDS_SEQ, "NM_CDS"))
        best = ts.longest_cds_transcript()
        assert best.accession == "NM_CDS"

    def test_get_metadata_cached(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record(_CDS_SEQ, "NM_001"))
        m1 = ts.get_metadata()
        m2 = ts.get_metadata()
        assert m1 is m2  # Cached

    def test_summary(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("ATGC" * 100, "NM_001"))
        ts.add_transcript(_make_record("ATGC" * 80, "NM_002"))
        s = ts.summary()
        assert s["gene_name"] == "BRCA1"
        assert s["isoform_count"] == 2
        assert s["longest_length"] == 400

    def test_shared_regions_single(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("ATGC" * 100, "NM_001"))
        shared = ts.shared_regions()
        assert len(shared) == 1
        assert shared[0] == (0, 400)

    def test_shared_regions_multiple(self):
        ts = TranscriptSet(gene_name="BRCA1", organism="Homo sapiens")
        ts.add_transcript(_make_record("ATGCATGC" * 50, "NM_001"))
        ts.add_transcript(_make_record("ATGCATGC" * 40, "NM_002"))
        shared = ts.shared_regions(min_overlap=20)
        assert len(shared) >= 1


# ---------------------------------------------------------------------------
# generate_quality_report (batch)
# ---------------------------------------------------------------------------


class TestGenerateQualityReport:
    def test_basic(self):
        records = [
            _make_record(_CDS_SEQ, "NM_001"),
            _make_record("GCT" * 100, "NM_002"),
        ]
        batch = generate_quality_report(records, min_length=100)
        assert batch.total_count == 2
        assert batch.valid_count == 2
        assert batch.transcripts_with_cds >= 1
        assert batch.transcripts_without_cds >= 1
        assert len(batch.per_transcript) == 2

    def test_with_invalid(self):
        records = [
            _make_record(_CDS_SEQ, "NM_001"),
            _make_record("ATGC", "NM_SHORT"),
        ]
        batch = generate_quality_report(records, min_length=100)
        assert batch.valid_count == 1
        assert batch.invalid_count == 1
        assert len(batch.flagged) >= 1

    def test_require_cds(self):
        records = [
            _make_record("GCT" * 100, "NM_NOCDS"),
        ]
        batch = generate_quality_report(records, min_length=100, require_cds=True)
        assert batch.valid_count == 0

    def test_batch_to_dict(self):
        records = [_make_record(_CDS_SEQ, "NM_001")]
        batch = generate_quality_report(records, min_length=100)
        d = batch.to_dict()
        assert d["total_count"] == 1
        assert d["valid_count"] == 1
        assert "per_transcript" in d
        assert "flagged" in d

    def test_empty(self):
        batch = generate_quality_report([])
        assert batch.total_count == 0
        assert batch.valid_count == 0
        assert batch.min_length == 0
