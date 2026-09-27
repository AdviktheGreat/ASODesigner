"""Shared test fixtures for ASODesigner tests."""

import tempfile
from pathlib import Path

import pytest

from aso_designer.transcript_fetcher import TranscriptRecord


@pytest.fixture
def tmp_cache_dir() -> Path:
    """Provide a temporary cache directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def sample_transcript_record() -> TranscriptRecord:
    """A small synthetic transcript for testing."""
    # Contains an ORF: ATG...TAA starting at position 3
    seq = (
        "AAA"  # 5' UTR (3 nt)
        "ATG"  # Start codon
        "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"  # 42 nt coding
        "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
        "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
        "TAA"  # Stop codon
        "TTTTTTTTTT"  # 3' UTR (10 nt)
    )
    return TranscriptRecord(
        accession="NM_TEST.1",
        sequence=seq,
        gene_name="TESTGENE",
        organism="Homo sapiens",
        description="NM_TEST.1 Homo sapiens test gene (TESTGENE), mRNA",
    )


@pytest.fixture
def sample_fasta_content() -> str:
    """FASTA content for testing local file loading."""
    return (
        ">NM_FAKE01.1 Homo sapiens test gene A (GENEA), mRNA\n"
        "ATGAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCC\n"
        "AAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCTAA\n"
        ">NM_FAKE02.1 Homo sapiens test gene B (GENEB), mRNA\n"
        "ATGCCCgggTTTAAAcccgggTTTAAAcccgggTTTAAAcccgggTTTAAAcccgggTTTAAA\n"
        "cccgggTTTAAAcccgggTTTAAAcccgggTTTAAAcccgggTTTAAAcccgggTTTAAATAA\n"
    )


@pytest.fixture
def sample_fasta_file(sample_fasta_content, tmp_path) -> Path:
    """Write a temporary FASTA file for testing."""
    fasta_path = tmp_path / "test_transcripts.fasta"
    fasta_path.write_text(sample_fasta_content)
    return fasta_path
