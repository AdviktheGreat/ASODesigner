"""Tests for off-target analysis (reference, search, scoring)."""

import pytest

from aso_designer.off_target import (
    OffTargetHit,
    OffTargetResult,
    OffTargetSearchEngine,
    TranscriptEntry,
    TranscriptomeReference,
    _hamming_distance,
    calculate_off_target_score,
    classify_hit,
    classify_off_target,
    filter_hits,
    score_off_target,
    score_off_target_batch,
    summarize_hits,
)
from aso_designer.utils.exceptions import OffTargetError

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_reference() -> TranscriptomeReference:
    """A small 2-transcript reference for testing."""
    ref = TranscriptomeReference(kmer_size=6)
    ref.add_transcript(
        TranscriptEntry(
            "NM_001.1",
            "GENEA",
            "ATGCATGCATGCATGCATGCGGCCATGCATGCATGC",
            36,
        )
    )
    ref.add_transcript(
        TranscriptEntry(
            "NM_002.1",
            "GENEB",
            "TTTTAAAATTTTAAAATTTTGGCCATGCATGCATGC",
            36,
        )
    )
    ref.build_index()
    return ref


@pytest.fixture
def clean_reference() -> TranscriptomeReference:
    """A reference with no off-target sites for the test ASO."""
    ref = TranscriptomeReference(kmer_size=6)
    ref.add_transcript(
        TranscriptEntry(
            "NM_CLEAN.1",
            "CLEANGENE",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            36,
        )
    )
    ref.build_index()
    return ref


# ---------------------------------------------------------------------------
# Hamming distance
# ---------------------------------------------------------------------------


class TestHammingDistance:
    def test_identical(self):
        assert _hamming_distance("ATGC", "ATGC") == 0

    def test_one_mismatch(self):
        assert _hamming_distance("ATGC", "ATGT") == 1

    def test_all_mismatches(self):
        assert _hamming_distance("ATGC", "GCTA") == 4

    def test_empty(self):
        assert _hamming_distance("", "") == 0


# ---------------------------------------------------------------------------
# Transcriptome reference
# ---------------------------------------------------------------------------


class TestTranscriptomeReference:
    def test_add_transcript(self):
        ref = TranscriptomeReference()
        ref.add_transcript(TranscriptEntry("NM_001.1", "GENEA", "ATGCATGCATGC", 12))
        assert ref.num_transcripts == 1
        assert ref.total_length == 12

    def test_add_multiple(self):
        ref = TranscriptomeReference()
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGC" * 5, 20))
        ref.add_transcript(TranscriptEntry("NM_002.1", "B", "TTTT" * 5, 20))
        assert ref.num_transcripts == 2
        assert ref.total_length == 40

    def test_invalid_sequence(self):
        ref = TranscriptomeReference()
        with pytest.raises(OffTargetError, match="Invalid sequence"):
            ref.add_transcript(TranscriptEntry("NM_001.1", "X", "XYZP", 4))

    def test_build_index(self):
        ref = TranscriptomeReference(kmer_size=4)
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGCATGC", 8))
        ref.build_index()
        assert ref._indexed
        assert len(ref._index) > 0

    def test_ensure_index(self):
        ref = TranscriptomeReference(kmer_size=4)
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGCATGC", 8))
        assert not ref._indexed
        ref.ensure_index()
        assert ref._indexed

    def test_get_transcript(self):
        ref = TranscriptomeReference()
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGC", 4))
        t = ref.get_transcript(0)
        assert t.transcript_id == "NM_001.1"

    def test_get_transcript_out_of_range(self):
        ref = TranscriptomeReference()
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGC", 4))
        with pytest.raises(OffTargetError, match="out of range"):
            ref.get_transcript(5)

    def test_get_sequence(self):
        ref = TranscriptomeReference()
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGCATGC", 8))
        assert ref.get_sequence(0, 2, 6) == "GCAT"

    def test_summary(self):
        ref = TranscriptomeReference(kmer_size=6)
        ref.add_transcript(TranscriptEntry("NM_001.1", "A", "ATGCATGCATGC", 12))
        ref.build_index()
        s = ref.summary()
        assert s["num_transcripts"] == 1
        assert s["total_length"] == 12
        assert s["indexed"] is True
        assert s["unique_kmers"] > 0

    def test_add_from_records(self):
        """Test adding from TranscriptRecord-like objects."""
        from aso_designer.transcript_fetcher import TranscriptRecord

        record = TranscriptRecord(
            accession="NM_TEST.1",
            sequence="ATGCATGCATGC",
            gene_name="TESTGENE",
            organism="Homo sapiens",
            description="test",
        )
        ref = TranscriptomeReference()
        ref.add_from_records([record])
        assert ref.num_transcripts == 1
        assert ref.get_transcript(0).gene_name == "TESTGENE"


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------


class TestSearchEngine:
    def test_basic_search(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        assert len(hits) > 0
        assert all(h.mismatches <= 3 for h in hits)

    def test_hits_sorted_by_mismatches(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        for i in range(len(hits) - 1):
            assert hits[i].mismatches <= hits[i + 1].mismatches

    def test_on_target_flag(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        hits = engine.search(
            "ATGCATGCATGCATGCATGC",
            on_target_transcript="NM_001.1",
            on_target_position=2,
        )
        on_target = [h for h in hits if h.is_on_target]
        assert len(on_target) >= 1
        assert on_target[0].transcript_id == "NM_001.1"
        assert on_target[0].position == 2

    def test_no_hits(self, clean_reference):
        engine = OffTargetSearchEngine(clean_reference, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        assert len(hits) == 0

    def test_zero_mismatch_hit(self, small_reference):
        """An ASO whose rev_comp exactly matches a transcript region."""
        # rev_comp("ATGCATGCATGCATGCATGC") = "GCATGCATGCATGCATGCAT"
        # Transcript 1 has "GCATGCATGCATGCATGCAT" at position 2 (with 2 mismatches)
        # Let's create a transcript with an exact match
        ref = TranscriptomeReference(kmer_size=6)
        ref.add_transcript(TranscriptEntry("NM_EXACT.1", "EXACT", "GCATGCATGCATGCATGCATAAAA", 24))
        ref.build_index()
        engine = OffTargetSearchEngine(ref, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        perfect = [h for h in hits if h.mismatches == 0]
        assert len(perfect) >= 1
        assert perfect[0].position == 0

    def test_seed_mismatch_count(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3, seed_start=2, seed_end=8)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        for h in hits:
            assert h.seed_mismatches >= 0
            assert h.seed_mismatches <= h.mismatches or h.mismatches == 0

    def test_max_mismatches_filter(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=0)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        for h in hits:
            assert h.mismatches == 0

    def test_short_aso_error(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        with pytest.raises(OffTargetError, match="shorter than k-mer"):
            engine.search("ATGC")  # 4 nt < kmer_size 6

    def test_batch_search(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        results = engine.search_batch(["ATGCATGCATGCATGCATGC", "TTTTAAAATTTTAAAATTTT"])
        assert 0 in results
        assert 1 in results
        assert isinstance(results[0], list)

    def test_aligned_sequence_correct(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        for h in hits:
            assert len(h.aligned_sequence) == 20


# ---------------------------------------------------------------------------
# Result processing
# ---------------------------------------------------------------------------


class TestClassifyHit:
    def test_on_target(self):
        hit = OffTargetHit("NM_001.1", "GENEA", 0, 0, 0, "ATGC", is_on_target=True)
        assert classify_hit(hit) == "on_target"

    def test_perfect(self):
        hit = OffTargetHit("NM_002.1", "GENEB", 0, 0, 0, "ATGC")
        assert classify_hit(hit) == "perfect"

    def test_seed_match(self):
        hit = OffTargetHit("NM_002.1", "GENEB", 0, 2, 0, "ATGC")
        assert classify_hit(hit) == "seed_match"

    def test_seed_mismatch(self):
        hit = OffTargetHit("NM_002.1", "GENEB", 0, 2, 1, "ATGC")
        assert classify_hit(hit) == "seed_mismatch"


class TestFilterHits:
    def test_exclude_on_target(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 0, 0, "ATGC", is_on_target=True),
            OffTargetHit("NM_002.1", "B", 0, 1, 0, "ATGC"),
        ]
        filtered = filter_hits(hits, exclude_on_target=True)
        assert len(filtered) == 1
        assert filtered[0].transcript_id == "NM_002.1"

    def test_keep_on_target(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 0, 0, "ATGC", is_on_target=True),
            OffTargetHit("NM_002.1", "B", 0, 1, 0, "ATGC"),
        ]
        filtered = filter_hits(hits, exclude_on_target=False)
        assert len(filtered) == 2

    def test_max_mismatches(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 1, 0, "ATGC"),
            OffTargetHit("NM_002.1", "B", 0, 3, 0, "ATGC"),
            OffTargetHit("NM_003.1", "C", 0, 5, 0, "ATGC"),
        ]
        filtered = filter_hits(hits, max_mismatches=3)
        assert len(filtered) == 2

    def test_max_seed_mismatches(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 2, 0, "ATGC"),
            OffTargetHit("NM_002.1", "B", 0, 2, 1, "ATGC"),
            OffTargetHit("NM_003.1", "C", 0, 2, 2, "ATGC"),
        ]
        filtered = filter_hits(hits, max_seed_mismatches=1)
        assert len(filtered) == 2


class TestSummarizeHits:
    def test_basic_summary(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 0, 0, "ATGC", is_on_target=True),
            OffTargetHit("NM_002.1", "B", 0, 0, 0, "ATGC"),
            OffTargetHit("NM_003.1", "C", 0, 2, 0, "ATGC"),
            OffTargetHit("NM_004.1", "D", 0, 3, 1, "ATGC"),
        ]
        result = summarize_hits("ATGC", hits)
        assert result.total_hits == 4
        assert result.off_target_hits == 3
        assert result.perfect_hits == 1
        assert result.seed_match_hits == 1
        assert result.seed_mismatch_hits == 1
        assert result.min_mismatches == 0
        assert result.has_on_target is True

    def test_no_hits(self):
        result = summarize_hits("ATGC", [])
        assert result.total_hits == 0
        assert result.off_target_hits == 0
        assert result.min_mismatches == -1
        assert result.has_on_target is False

    def test_all_on_target(self):
        hits = [
            OffTargetHit("NM_001.1", "A", 0, 0, 0, "ATGC", is_on_target=True),
        ]
        result = summarize_hits("ATGC", hits)
        assert result.off_target_hits == 0
        assert result.has_on_target is True


# ---------------------------------------------------------------------------
# Off-target scoring
# ---------------------------------------------------------------------------


class TestCalculateOffTargetScore:
    def test_no_off_targets(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=1,
            off_target_hits=0,
            perfect_hits=0,
            seed_match_hits=0,
            seed_mismatch_hits=0,
            min_mismatches=-1,
            hits=[],
            has_on_target=True,
        )
        assert calculate_off_target_score(result) == 1.0

    def test_one_perfect_hit(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=2,
            off_target_hits=1,
            perfect_hits=1,
            seed_match_hits=0,
            seed_mismatch_hits=0,
            min_mismatches=0,
            hits=[],
            has_on_target=True,
        )
        score = calculate_off_target_score(result, max_hits=50)
        assert score == pytest.approx(1.0 - 1.0 / 50)

    def test_one_seed_match(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=2,
            off_target_hits=1,
            perfect_hits=0,
            seed_match_hits=1,
            seed_mismatch_hits=0,
            min_mismatches=2,
            hits=[],
            has_on_target=True,
        )
        score = calculate_off_target_score(result, seed_mismatch_weight=2.0, max_hits=50)
        assert score == pytest.approx(1.0 - (0.5 * 2.0) / 50)

    def test_one_seed_mismatch(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=2,
            off_target_hits=1,
            perfect_hits=0,
            seed_match_hits=0,
            seed_mismatch_hits=1,
            min_mismatches=2,
            hits=[],
            has_on_target=True,
        )
        score = calculate_off_target_score(result, max_hits=50)
        assert score == pytest.approx(1.0 - 0.25 / 50)

    def test_many_hits_saturates_to_zero(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=100,
            off_target_hits=100,
            perfect_hits=100,
            seed_match_hits=0,
            seed_mismatch_hits=0,
            min_mismatches=0,
            hits=[],
            has_on_target=False,
        )
        score = calculate_off_target_score(result, max_hits=50)
        assert score == 0.0

    def test_score_in_range(self):
        result = OffTargetResult(
            aso_sequence="ATGC",
            total_hits=10,
            off_target_hits=5,
            perfect_hits=2,
            seed_match_hits=2,
            seed_mismatch_hits=1,
            min_mismatches=0,
            hits=[],
            has_on_target=True,
        )
        score = calculate_off_target_score(result, max_hits=50)
        assert 0.0 <= score <= 1.0


class TestClassifyOffTarget:
    def test_specific(self):
        assert classify_off_target(0.9) == "specific"

    def test_at_specific_threshold(self):
        assert classify_off_target(0.8) == "specific"

    def test_moderate(self):
        assert classify_off_target(0.5) == "moderate"

    def test_at_moderate_threshold(self):
        assert classify_off_target(0.3) == "moderate"

    def test_non_specific(self):
        assert classify_off_target(0.1) == "non-specific"

    def test_custom_thresholds(self):
        assert (
            classify_off_target(0.6, specific_threshold=0.7, moderate_threshold=0.2) == "moderate"
        )


class TestScoreOffTarget:
    def test_basic(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        hits = engine.search(
            "ATGCATGCATGCATGCATGC",
            on_target_transcript="NM_001.1",
            on_target_position=2,
        )
        result, score, classification = score_off_target("ATGCATGCATGCATGCATGC", hits)
        assert isinstance(result, OffTargetResult)
        assert 0.0 <= score <= 1.0
        assert classification in ("specific", "moderate", "non-specific")

    def test_clean_candidate(self, clean_reference):
        engine = OffTargetSearchEngine(clean_reference, max_mismatches=3)
        hits = engine.search("ATGCATGCATGCATGCATGC")
        result, score, classification = score_off_target("ATGCATGCATGCATGCATGC", hits)
        assert score == 1.0
        assert classification == "specific"

    def test_batch(self, small_reference):
        engine = OffTargetSearchEngine(small_reference, max_mismatches=3)
        all_hits = engine.search_batch(["ATGCATGCATGCATGCATGC", "TTTTAAAATTTTAAAATTTT"])
        results = score_off_target_batch(["ATGCATGCATGCATGCATGC", "TTTTAAAATTTTAAAATTTT"], all_hits)
        assert 0 in results
        assert 1 in results
        for _, score, cls in results.values():
            assert 0.0 <= score <= 1.0
            assert cls in ("specific", "moderate", "non-specific")
