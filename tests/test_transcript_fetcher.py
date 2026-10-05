"""Tests for the NCBI Entrez transcript fetcher with mocked HTTP responses."""

import json
from unittest.mock import patch

import pytest
import responses

from aso_designer.transcript_fetcher import EntrezFetcher, TranscriptRecord
from aso_designer.utils.exceptions import TranscriptFetchError

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@responses.activate
def test_resolve_gene_to_transcripts_success(tmp_path) -> None:
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

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
    accessions = fetcher.resolve_gene_to_transcripts("GENEA", "Homo sapiens")

    assert "NM_000001.3" in accessions
    # NR_ accessions (non-coding) should be excluded
    assert all(a.startswith("NM_") for a in accessions)


@responses.activate
def test_resolve_gene_no_results(tmp_path) -> None:
    """Gene resolution raises TranscriptFetchError when no transcripts found."""
    esearch_body = json.dumps({"esearchresult": {"idlist": [], "count": "0"}})
    responses.add(
        responses.GET,
        f"{_EUTILS}/esearch.fcgi",
        body=esearch_body,
        status=200,
        content_type="application/json",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
    with pytest.raises(TranscriptFetchError, match="No RefSeq transcripts"):
        fetcher.resolve_gene_to_transcripts("FAKEGENE", "Homo sapiens")


@responses.activate
def test_fetch_by_accession_success(tmp_path) -> None:
    """Accession-based fetch returns a valid TranscriptRecord."""
    fasta_text = ">NM_TEST01.1 Homo sapiens test gene (TESTG) mRNA\nATGAAAGGGTTTCCCTAA\n"
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        body=fasta_text,
        status=200,
        content_type="text/plain",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
    record = fetcher.fetch_by_accession("NM_TEST01.1")

    assert isinstance(record, TranscriptRecord)
    assert record.accession == "NM_TEST01.1"
    assert record.sequence == "ATGAAAGGGTTTCCCTAA"
    assert record.length == 18


@responses.activate
def test_fetch_by_accession_retry_on_failure(tmp_path) -> None:
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

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
    record = fetcher.fetch_by_accession("NM_TEST02.1")
    assert record.accession == "NM_TEST02.1"


@responses.activate
def test_fetch_by_accession_all_retries_fail(tmp_path) -> None:
    """Fetcher raises TranscriptFetchError after all retries are exhausted."""
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        status=503,
    )
    # responses library reuses the same response for repeated calls

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
    with pytest.raises(TranscriptFetchError, match="NCBI request failed"):
        fetcher.fetch_by_accession("NM_FAIL01.1")


@responses.activate
def test_fetch_by_accession_empty_response(tmp_path) -> None:
    """Empty EFetch response raises TranscriptFetchError."""
    responses.add(
        responses.GET,
        f"{_EUTILS}/efetch.fcgi",
        body="",
        status=200,
        content_type="text/plain",
    )

    fetcher = EntrezFetcher(api_key="fake_key", cache_dir=str(tmp_path))
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


# ---------------------------------------------------------------------------
# Cache behavior
# ---------------------------------------------------------------------------


def _write_fasta(path, seq_id="NM_TEST01.1", seq="ATGAAAGGGTTTCCCTAA"):
    with open(path, "w") as f:
        f.write(f">{seq_id} Homo sapiens test gene (TESTG) mRNA\n{seq}\n")
    return path


class TestTranscriptCache:
    def test_fetch_then_cache_hit(self, tmp_path):
        """A second fetch of the same accession is served from the cache."""
        fasta_text = ">NM_TEST01.1 Homo sapiens test gene (TESTG) mRNA\nATGAAAGGGTTTCCCTAA\n"
        with responses.RequestsMock() as rsps:
            rsps.add(
                responses.GET,
                f"{_EUTILS}/efetch.fcgi",
                body=fasta_text,
                status=200,
                content_type="text/plain",
            )
            fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path / "cache"))
            record = fetcher.fetch_by_accession("NM_TEST01.1")
        assert record.gene_name == "TESTG"

        # No registered HTTP responses: a second call must come from cache.
        fetcher2 = EntrezFetcher(api_key="k", cache_dir=str(tmp_path / "cache"))
        cached = fetcher2.fetch_by_accession("NM_TEST01.1")
        assert cached.sequence == record.sequence

    def test_corrupted_cache_returns_none(self, tmp_path, caplog):
        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        cache_file = fetcher._cache_key("acc:NM_X.1")
        cache_file.write_text("{corrupt json", encoding="utf-8")
        assert fetcher._load_cache("acc:NM_X.1") is None

    def test_save_and_load_roundtrip(self, tmp_path):
        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        rec = TranscriptRecord(
            accession="NM_R.1",
            sequence="ATGCATGCATGC",
            gene_name="RGENE",
            organism="Homo sapiens",
        )
        fetcher._save_cache("acc:NM_R.1", [rec])
        loaded = fetcher._load_cache("acc:NM_R.1")
        assert loaded is not None
        assert loaded[0].gene_name == "RGENE"

    def test_expired_cache_invalid(self, tmp_path):
        from datetime import datetime, timedelta

        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        rec = TranscriptRecord(
            accession="NM_OLD.1",
            sequence="ATGCATGCATGC",
            gene_name="OLD",
            organism="Homo sapiens",
        )
        fetcher._save_cache("acc:NM_OLD.1", [rec])
        cache_file = fetcher._cache_key("acc:NM_OLD.1")
        old = datetime.now() - fetcher.cache_ttl - timedelta(hours=1)
        import os

        os.utime(cache_file, (old.timestamp(), old.timestamp()))
        assert fetcher._is_cache_valid(cache_file) is False


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


class TestDescriptionParsing:
    def test_gene_from_parentheses(self):
        assert (
            EntrezFetcher._parse_gene_from_description(
                "NM_000109.4 Homo sapiens dystrophin (DMD), mRNA", "NM_000109.4"
            )
            == "DMD"
        )

    def test_gene_fallback_to_accession_prefix(self):
        assert EntrezFetcher._parse_gene_from_description("no parens here", "NM_12345.1") == "NM"

    def test_gene_fallback_no_underscore(self):
        assert EntrezFetcher._parse_gene_from_description("nothing", "plain") == "unknown"

    @pytest.mark.parametrize(
        ("description", "expected"),
        [
            ("Homo sapiens dystrophin mRNA", "Homo sapiens"),
            ("Mus musculus gene mRNA", "Mus musculus"),
            ("Rattus norvegicus transcript", "Rattus norvegicus"),
            ("Danio rerio mRNA", "Danio rerio"),
            ("Drosophila melanogaster mRNA", "Drosophila melanogaster"),
            ("an unremarkable header", "unknown"),
        ],
    )
    def test_organism_from_description(self, description, expected):
        assert EntrezFetcher._parse_organism_from_description(description) == expected


# ---------------------------------------------------------------------------
# Batch and ortholog fetching
# ---------------------------------------------------------------------------


class TestBatchAndOrthologs:
    @responses.activate
    def test_fetch_transcripts_skips_failures(self, tmp_path):
        """fetch_transcripts keeps successes and warns on failures."""
        fasta_text = ">NM_GOOD.1 Homo sapiens test gene (GOODG) mRNA\nATGAAAGGGTTTCCCTAA\n"
        rsps = responses.RequestsMock(assert_all_requests_are_fired=False)
        with rsps:
            rsps.add(
                responses.GET,
                f"{_EUTILS}/efetch.fcgi",
                body=fasta_text,
                status=200,
                content_type="text/plain",
            )
            rsps.add(responses.GET, f"{_EUTILS}/efetch.fcgi", status=500)
            rsps.add(responses.GET, f"{_EUTILS}/efetch.fcgi", status=500)
            rsps.add(responses.GET, f"{_EUTILS}/efetch.fcgi", status=500)

            fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path / "c"))
            records = fetcher.fetch_transcripts(["NM_GOOD.1", "NM_BAD.1"])
        assert len(records) == 1
        assert records[0].accession == "NM_GOOD.1"

    @responses.activate
    def test_unparseable_fasta_raises(self, tmp_path):
        # Header-only FASTA yields an empty sequence, which fails validation.
        responses.add(
            responses.GET,
            f"{_EUTILS}/efetch.fcgi",
            body=">NM_BADFASTA.1 Homo sapiens\n",
            status=200,
            content_type="text/plain",
        )
        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        with pytest.raises(TranscriptFetchError):
            fetcher.fetch_by_accession("NM_BADFASTA.1")

    def test_fetch_orthologs_mixed_success(self, tmp_path):
        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        record = TranscriptRecord(
            accession="NM_MOUSE.1",
            sequence="ATGC" * 30,
            gene_name="GENE",
            organism="Mus musculus",
        )
        with patch.object(
            fetcher,
            "resolve_gene_to_transcripts",
            side_effect=lambda gene, org: (
                ["NM_MOUSE.1"]
                if org == "Mus musculus"
                else (_ for _ in ()).throw(TranscriptFetchError("none found"))
            ),
        ):
            with patch.object(fetcher, "fetch_transcripts", return_value=[record]):
                result = fetcher.fetch_orthologs("GENE", ["Mus musculus", "Rattus norvegicus"])
        assert [r.accession for r in result["Mus musculus"]] == ["NM_MOUSE.1"]
        assert result["Rattus norvegicus"] == []

    def test_fetch_gene_combines_resolution_and_fetch(self, tmp_path):
        fetcher = EntrezFetcher(api_key="k", cache_dir=str(tmp_path))
        record = TranscriptRecord(
            accession="NM_G.1",
            sequence="ATGC" * 30,
            gene_name="GENE",
            organism="Homo sapiens",
        )
        with patch.object(fetcher, "resolve_gene_to_transcripts", return_value=["NM_G.1"]):
            with patch.object(fetcher, "fetch_transcripts", return_value=[record]) as ft:
                records = fetcher.fetch_gene("GENE")
        assert [r.accession for r in records] == ["NM_G.1"]
        ft.assert_called_once_with(["NM_G.1"])
