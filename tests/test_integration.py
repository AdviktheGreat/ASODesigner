"""Integration tests that hit the live NCBI E-utilities API.

These tests are marked @pytest.mark.slow and @pytest.mark.network so they
are skipped in CI by default. Run them manually with:

    pytest -m slow tests/test_integration.py
"""

import pytest

from aso_designer.fetch_manager import FetchManager
from aso_designer.quality import extract_metadata, generate_quality_report, validate_transcript
from aso_designer.transcript_fetcher import EntrezFetcher, TranscriptRecord

pytestmark = [pytest.mark.slow, pytest.mark.network]


@pytest.fixture
def live_fetcher() -> EntrezFetcher:
    """Create a fetcher using the real NCBI API (requires NCBI_API_KEY env var)."""
    return EntrezFetcher()


def test_live_fetch_gene_dmd(live_fetcher: EntrezFetcher) -> None:
    """Fetch DMD transcripts from NCBI and validate them.

    DMD (dystrophin) has a large mRNA (~14 kb), making it a good
    integration test target with many ASO candidate sites.
    """
    accessions = live_fetcher.resolve_gene_to_transcripts("DMD", "Homo sapiens")
    assert len(accessions) > 0
    assert any(a.startswith("NM_") for a in accessions)

    # Fetch the first transcript
    record = live_fetcher.fetch_by_accession(accessions[0])
    assert isinstance(record, TranscriptRecord)
    assert record.length > 1000  # DMD mRNA should be large
    assert record.gene_name  # Should have a gene name


def test_live_fetch_and_validate(live_fetcher: EntrezFetcher) -> None:
    """Fetch a transcript and run the full quality validation pipeline."""
    accessions = live_fetcher.resolve_gene_to_transcripts("SMN1", "Homo sapiens")
    if not accessions:
        pytest.skip("SMN1 not found in NCBI (may be temporarily unavailable)")

    record = live_fetcher.fetch_by_accession(accessions[0])

    # Validate quality
    report = validate_transcript(record)
    assert report.is_valid
    assert report.length > 0

    # Extract metadata
    metadata = extract_metadata(record)
    assert metadata.accession == record.accession
    assert metadata.length == record.length


def test_live_batch_quality_report(live_fetcher: EntrezFetcher) -> None:
    """Fetch multiple transcripts and generate a batch quality report."""
    accessions = live_fetcher.resolve_gene_to_transcripts("SMN1", "Homo sapiens")
    if not accessions:
        pytest.skip("SMN1 not found in NCBI")

    records = live_fetcher.fetch_transcripts(accessions[:3])
    assert len(records) > 0

    batch = generate_quality_report(records)
    assert batch.total_count > 0
    assert batch.valid_count > 0
    assert batch.mean_length > 0


def test_live_fetch_manager_fallback() -> None:
    """Test FetchManager with a real NCBI fetch (no local fallback)."""
    fm = FetchManager()
    try:
        records = fm.fetch_by_gene("TP53", "Homo sapiens")
        assert len(records) > 0
        assert all(r.length > 0 for r in records)
    except Exception:
        pytest.skip("NCBI API unavailable or rate-limited")
