"""Tests for the NCBI Entrez transcript fetcher with mocked HTTP responses."""

import json

import pytest
import responses

from aso_designer.transcript_fetcher import EntrezFetcher, TranscriptRecord
from aso_designer.utils.exceptions import TranscriptFetchError

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@responses.activate
def test_resolve_gene_to_transcripts_success() -> None:
    """Gene symbol resolution returns RefSeq NM_ accessions via ESearch/ESummary."""
    # ESearch response
    esearch_body = json.dumps({"esearchresult": {"idlist": ["12345", "67890"], "count": "2"}})
    responses.add(
        responses.GET,
        f"{_EUTILS}/esearch.fcgi",
        body=esearch_body,
        status=200,
        content_type="application/json",
    )

    # ESummary response
    esummary_body = json.dumps(
        {
            "result": {
                "12345": {"accessionversion": "NM_000001.3", "title": "test gene (GENEA) mRNA"},
                "67890": {"accessionversion": "NR_000002.1", "title": "test RNA"},
            }
        }
    )
    responses.add(
        responses.GET,
        f"{_EUTILS}/esummary.fcgi",
        body=esummary_body,
        status=200,
        content_type="application/json",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_1")
    accessions = fetcher.resolve_gene_to_transcripts("GENEA", "Homo sapiens")

    assert "NM_000001.3" in accessions
    # NR_ accessions (non-coding) should be excluded
    assert all(a.startswith("NM_") for a in accessions)


@responses.activate
def test_resolve_gene_no_results() -> None:
    """Gene resolution raises TranscriptFetchError when no transcripts found."""
    esearch_body = json.dumps({"esearchresult": {"idlist": [], "count": "0"}})
    responses.add(
        responses.GET,
        f"{_EUTILS}/esearch.fcgi",
        body=esearch_body,
        status=200,
        content_type="application/json",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_2")
    with pytest.raises(TranscriptFetchError, match="No RefSeq transcripts"):
        fetcher.resolve_gene_to_transcripts("FAKEGENE", "Homo sapiens")


@responses.activate
def test_fetch_by_accession_success() -> None:
    """Accession-based fetch returns a valid TranscriptRecord."""
    fasta_text = ">NM_TEST01.1 Homo sapiens test gene (TESTG) mRNA\nATGAAAGGGTTTCCCTAA\n"
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        body=fasta_text,
        status=200,
        content_type="text/plain",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_3")
    record = fetcher.fetch_by_accession("NM_TEST01.1")

    assert isinstance(record, TranscriptRecord)
    assert record.accession == "NM_TEST01.1"
    assert record.sequence == "ATGAAAGGGTTTCCCTAA"
    assert record.length == 18


@responses.activate
def test_fetch_by_accession_retry_on_failure() -> None:
    """Fetcher retries on transient HTTP errors before raising."""
    fasta_text = ">NM_TEST02.1 Homo sapiens test gene (TESTG2) mRNA\nATGAAAGGGTTTCCCTAA\n"
    # First two calls fail with 500, third succeeds
    responses.add(responses.GET, f"{_EUTILS}/efetch.fcgi", status=500)
    responses.add(responses.GET, f"{_EUTILS}/efetch.fcgi", status=500)
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        body=fasta_text,
        status=200,
        content_type="text/plain",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_4")
    record = fetcher.fetch_by_accession("NM_TEST02.1")
    assert record.accession == "NM_TEST02.1"


@responses.activate
def test_fetch_by_accession_all_retries_fail() -> None:
    """Fetcher raises TranscriptFetchError after all retries are exhausted."""
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        status=503,
    )
    # responses library reuses the same response for repeated calls

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_5")
    with pytest.raises(TranscriptFetchError, match="NCBI request failed"):
        fetcher.fetch_by_accession("NM_FAIL01.1")


@responses.activate
def test_fetch_by_accession_empty_response() -> None:
    """Empty EFetch response raises TranscriptFetchError."""
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        body="",
        status=200,
        content_type="text/plain",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir="/tmp/test_aso_cache_6")
    with pytest.raises(TranscriptFetchError, match="No data returned"):
        fetcher.fetch_by_accession("NM_EMPTY01.1")


def test_transcript_record_to_dict_roundtrip(sample_transcript_record: TranscriptRecord) -> None:
    """TranscriptRecord serializes and deserializes correctly."""
    d = sample_transcript_record.to_dict()
    restored = TranscriptRecord.from_dict(d)
    assert restored.accession == sample_transcript_record.accession
    assert restored.sequence == sample_transcript_record.sequence
    assert restored.gene_name == sample_transcript_record.gene_name
    assert restored.length == sample_transcript_record.length
