"""Tests for local FASTA file loading and validation."""

from pathlib import Path

import pytest

from aso_designer.local_loader import load_fasta, load_single_fasta
from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.utils.exceptions import TranscriptFetchError


def test_load_fasta_valid_file(sample_fasta_file: Path) -> None:
    """Loading a valid multi-record FASTA returns all records."""
    records = load_fasta(sample_fasta_file)

    assert len(records) == 2
    assert all(isinstance(r, TranscriptRecord) for r in records)
    assert records[0].accession == "NM_FAKE01.1"
    assert records[1].accession == "NM_FAKE02.1"


def test_load_fasta_parses_gene_name(sample_fasta_file: Path) -> None:
    """Gene name is extracted from FASTA header via parenthetical pattern."""
    records = load_fasta(sample_fasta_file)
    assert records[0].gene_name == "GENEA"
    assert records[1].gene_name == "GENEB"


def test_load_fasta_parses_organism(sample_fasta_file: Path) -> None:
    """Organism name is extracted from FASTA header."""
    records = load_fasta(sample_fasta_file)
    assert records[0].organism == "Homo sapiens"


def test_load_fasta_uppercase_sequences(sample_fasta_file: Path) -> None:
    """Sequences are returned in uppercase."""
    records = load_fasta(sample_fasta_file)
    for r in records:
        assert r.sequence == r.sequence.upper()


def test_load_single_fasta(sample_fasta_file: Path) -> None:
    """load_single_fasta returns the first record."""
    record = load_single_fasta(sample_fasta_file)
    assert record.accession == "NM_FAKE01.1"


def test_load_fasta_min_length_filter(tmp_path: Path) -> None:
    """Sequences below min_length are skipped."""
    fasta = tmp_path / "short.fasta"
    fasta.write_text(">SHORT01.1 short seq\nATGC\n>LONG01.1 long seq\n" + "ATGC" * 50 + "\n")
    records = load_fasta(fasta, min_length=100)
    assert len(records) == 1
    assert records[0].accession == "LONG01.1"


def test_load_fasta_file_not_found() -> None:
    """Missing file raises TranscriptFetchError."""
    with pytest.raises(TranscriptFetchError, match="FASTA file not found"):
        load_fasta("/nonexistent/path.fasta")


def test_load_fasta_empty_file(tmp_path: Path) -> None:
    """Empty FASTA file raises TranscriptFetchError."""
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")
    with pytest.raises(TranscriptFetchError, match="No valid sequences"):
        load_fasta(fasta)


def test_load_fasta_invalid_chars_skipped(tmp_path: Path) -> None:
    """Records with invalid characters are skipped with a warning."""
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(">BAD01.1 bad seq\nATGCXYZATGC\n>GOOD01.1 good seq\n" + "ATGC" * 30 + "\n")
    records = load_fasta(fasta, min_length=50, allow_iupac=False)
    # The bad record should be skipped (XYZ not valid DNA)
    assert all("X" not in r.sequence for r in records)


def test_load_fasta_iupac_allowed(tmp_path: Path) -> None:
    """IUPAC ambiguity codes are accepted when allow_iupac=True."""
    fasta = tmp_path / "iupac.fasta"
    fasta.write_text(">IUPAC01.1 iupac seq\n" + "ATGCNRYSW" * 20 + "\n")
    records = load_fasta(fasta, min_length=50, allow_iupac=True)
    assert len(records) == 1
    assert "N" in records[0].sequence
