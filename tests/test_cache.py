"""Tests for cache management and FetchManager fallback chain."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from aso_designer.cache_manager import TranscriptCacheManager
from aso_designer.fetch_manager import FetchManager
from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.utils.exceptions import TranscriptFetchError

# ---------------------------------------------------------------------------
# TranscriptCacheManager tests
# ---------------------------------------------------------------------------


def test_cache_store_and_retrieve(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """Storing and retrieving a record round-trips correctly."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    cache.store("test_key", [sample_transcript_record], organism="Homo sapiens")

    retrieved = cache.retrieve("test_key", organism="Homo sapiens")
    assert retrieved is not None
    assert len(retrieved) == 1
    assert retrieved[0].accession == sample_transcript_record.accession
    assert retrieved[0].sequence == sample_transcript_record.sequence


def test_cache_miss_returns_none(tmp_cache_dir: Path) -> None:
    """Retrieving a non-existent key returns None."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    assert cache.retrieve("nonexistent_key") is None


def test_cache_invalidation(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """Invalidating a cache entry removes it."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    cache.store("invalidate_me", [sample_transcript_record])

    assert cache.invalidate("invalidate_me") is True
    assert cache.retrieve("invalidate_me") is None
    # Invalidating again returns False (already gone)
    assert cache.invalidate("invalidate_me") is False


def test_cache_ttl_expiry(tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord) -> None:
    """Cached entries past TTL return None."""
    from datetime import timedelta

    cache = TranscriptCacheManager(cache_base=tmp_cache_dir, default_ttl=168)
    cache.store("ttl_test", [sample_transcript_record])

    # Within TTL — should return data
    assert cache.retrieve("ttl_test") is not None

    # With zero TTL — should return None (expired)
    assert cache.retrieve("ttl_test", ttl=timedelta(seconds=0)) is None


def test_cache_organism_directory(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """Cache files are organized into organism-specific directories."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    cache.store("org_test", [sample_transcript_record], organism="Homo sapiens")

    org_dir = tmp_cache_dir / "Homo_sapiens"
    assert org_dir.exists()
    assert org_dir.is_dir()
    assert list(org_dir.glob("*.json"))


def test_cache_inspect(tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord) -> None:
    """Inspect returns a summary with per-organism details."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    cache.store("inspect1", [sample_transcript_record], organism="Homo sapiens")
    cache.store("inspect2", [sample_transcript_record], organism="Mus musculus")

    report = cache.inspect()
    assert "Homo sapiens" in report
    assert "Mus musculus" in report
    assert report["Homo sapiens"]["files"] >= 1


def test_cache_clear_organism(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """Clearing an organism removes only that organism's entries."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    cache.store("clear1", [sample_transcript_record], organism="Homo sapiens")
    cache.store("clear2", [sample_transcript_record], organism="Mus musculus")

    deleted = cache.clear_organism("Homo sapiens")
    assert deleted >= 1
    assert cache.retrieve("clear1", organism="Homo sapiens") is None
    assert cache.retrieve("clear2", organism="Mus musculus") is not None


def test_cache_total_entries(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """total_entries returns the count of all cached files."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    assert cache.total_entries() == 0
    cache.store("e1", [sample_transcript_record], organism="Homo sapiens")
    cache.store("e2", [sample_transcript_record], organism="Mus musculus")
    assert cache.total_entries() == 2


# ---------------------------------------------------------------------------
# FetchManager fallback chain tests
# ---------------------------------------------------------------------------


def test_fetch_manager_cache_hit(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """FetchManager returns cached data without hitting NCBI."""
    fm = FetchManager(cache_dir=tmp_cache_dir)
    # Pre-populate cache
    fm.cache.store("gene:DMD:Homo sapiens", [sample_transcript_record], "Homo sapiens")

    records = fm.fetch_by_gene("DMD", "Homo sapiens")
    assert len(records) == 1
    assert records[0].accession == sample_transcript_record.accession


def test_fetch_manager_local_fasta_fallback(tmp_cache_dir: Path, sample_fasta_file: Path) -> None:
    """FetchManager falls back to local FASTA when NCBI is unavailable."""
    fm = FetchManager(cache_dir=tmp_cache_dir, local_fasta=sample_fasta_file)

    # Mock the EntrezFetcher.fetch_gene to raise TranscriptFetchError
    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_gene.side_effect = TranscriptFetchError("Network down")
        records = fm.fetch_by_gene("ANYGENE", "Homo sapiens")

    assert len(records) == 2
    assert records[0].accession == "NM_FAKE01.1"


def test_fetch_manager_all_sources_fail(tmp_cache_dir: Path) -> None:
    """FetchManager raises TranscriptFetchError when all sources fail."""
    fm = FetchManager(cache_dir=tmp_cache_dir)

    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_gene.side_effect = TranscriptFetchError("Network down")
        with pytest.raises(TranscriptFetchError, match="All sources failed"):
            fm.fetch_by_gene("NOGENE", "Homo sapiens")


def test_fetch_manager_accession_cache_hit(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """FetchManager returns cached accession data without hitting NCBI."""
    fm = FetchManager(cache_dir=tmp_cache_dir)
    fm.cache.store("acc:NM_TEST.1", [sample_transcript_record], "Homo sapiens")

    record = fm.fetch_by_accession("NM_TEST.1", "Homo sapiens")
    assert record.accession == sample_transcript_record.accession


def test_fetch_from_local(sample_fasta_file: Path) -> None:
    """fetch_from_local loads a FASTA file directly."""
    fm = FetchManager()
    records = fm.fetch_from_local(sample_fasta_file)
    assert len(records) == 2


# ---------------------------------------------------------------------------
# Corruption and cache-wide utilities
# ---------------------------------------------------------------------------


def test_cache_corrupt_file_returns_none(tmp_cache_dir: Path) -> None:
    """A corrupted cache file is treated as a miss, not a crash."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    path = cache._cache_key_path("bad_key", "Homo sapiens")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not valid json", encoding="utf-8")
    assert cache.retrieve("bad_key", organism="Homo sapiens") is None


def test_cache_clear_all_and_total_size(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """clear_all removes every entry; size accounting tracks stored bytes."""
    cache = TranscriptCacheManager(cache_base=tmp_cache_dir)
    assert cache.total_size_bytes() == 0
    cache.store("a", [sample_transcript_record], organism="Homo sapiens")
    cache.store("b", [sample_transcript_record], organism="Mus musculus")
    assert cache.total_size_bytes() > 0

    deleted = cache.clear_all()
    assert deleted == 2
    assert cache.total_entries() == 0
    assert cache.total_size_bytes() == 0


# ---------------------------------------------------------------------------
# FetchManager NCBI-success and accession fallback paths
# ---------------------------------------------------------------------------


def test_fetch_manager_ncbi_success_caches(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """A successful NCBI gene fetch is stored in the cache."""
    fm = FetchManager(cache_dir=tmp_cache_dir)
    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_gene.return_value = [sample_transcript_record]
        records = fm.fetch_by_gene("DMD", "Homo sapiens")
    assert len(records) == 1

    fresh = FetchManager(cache_dir=tmp_cache_dir)
    assert fresh.cache.retrieve("gene:DMD:Homo sapiens", organism="Homo sapiens")


def test_fetch_manager_accession_success_caches(
    tmp_cache_dir: Path, sample_transcript_record: TranscriptRecord
) -> None:
    """A successful NCBI accession fetch is stored in the cache."""
    fm = FetchManager(cache_dir=tmp_cache_dir)
    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_by_accession.return_value = sample_transcript_record
        record = fm.fetch_by_accession("NM_TEST.1", "Homo sapiens")
    assert record.accession == "NM_TEST.1"

    fresh = FetchManager(cache_dir=tmp_cache_dir)
    assert fresh.cache.retrieve("acc:NM_TEST.1", organism="Homo sapiens")


def test_fetch_manager_accession_fasta_match(tmp_cache_dir: Path, sample_fasta_file: Path) -> None:
    """When NCBI fails, the accession fallback scans the local FASTA."""
    fm = FetchManager(cache_dir=tmp_cache_dir, local_fasta=sample_fasta_file)
    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_by_accession.side_effect = TranscriptFetchError("down")
        record = fm.fetch_by_accession("NM_FAKE01.1", "Homo sapiens")
    assert record.accession == "NM_FAKE01.1"


def test_fetch_manager_accession_fasta_no_match(
    tmp_cache_dir: Path, sample_fasta_file: Path
) -> None:
    """No FASTA record matching the accession → all sources failed."""
    fm = FetchManager(cache_dir=tmp_cache_dir, local_fasta=sample_fasta_file)
    with patch.object(FetchManager, "entrez", new_callable=MagicMock) as mock_entrez:
        mock_entrez.fetch_by_accession.side_effect = TranscriptFetchError("down")
        with pytest.raises(TranscriptFetchError, match="All sources failed"):
            fm.fetch_by_accession("NM_MISSING.9", "Homo sapiens")
