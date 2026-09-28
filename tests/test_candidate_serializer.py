"""Tests for candidate serialization, batch I/O, and deduplication."""

from pathlib import Path

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.candidate_serializer import (
    candidates_equal,
    deduplicate_candidates,
    load_candidates,
    load_candidates_batch,
    save_candidates,
    save_candidates_batch,
    validate_roundtrip,
)


def _make_candidates() -> list[ASOCandidate]:
    """Create a small set of test candidates."""
    return [
        ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, chemistry="gapmer"),
        ASOCandidate("GCTAGCTAGCTAGCTAGCTA", 20, 40, 20, gc_content=0.5, chemistry="gapmer"),
        ASOCandidate("TTTTAAAACCCCGGGGAAAA", 40, 60, 20, gc_content=0.5, chemistry="mixmer"),
    ]


def test_save_and_load_roundtrip(tmp_path: Path) -> None:
    """Save and load preserves all candidate data."""
    cands = _make_candidates()
    path = save_candidates(cands, tmp_path / "test.json", metadata={"gene": "DMD"})
    loaded, meta = load_candidates(path)
    assert len(loaded) == 3
    assert meta["gene"] == "DMD"
    for orig, load in zip(cands, loaded, strict=False):
        assert orig.sequence == load.sequence
        assert orig.start == load.start
        assert orig.chemistry == load.chemistry


def test_save_creates_parent_dirs(tmp_path: Path) -> None:
    """save_candidates creates parent directories."""
    cands = _make_candidates()
    path = tmp_path / "subdir" / "deeper" / "test.json"
    save_candidates(cands, path)
    assert path.exists()


def test_load_file_not_found() -> None:
    """Loading a missing file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_candidates("/nonexistent/path.json")


def test_load_invalid_format(tmp_path: Path) -> None:
    """Loading a file without 'candidates' key raises ValueError."""
    path = tmp_path / "bad.json"
    path.write_text('{"foo": "bar"}')
    with pytest.raises(ValueError, match="No 'candidates'"):
        load_candidates(path)


def test_batch_save_and_load(tmp_path: Path) -> None:
    """Batch save/load preserves all candidates across multiple files."""
    cands = _make_candidates()
    paths = save_candidates_batch(cands, tmp_path / "batch", batch_size=2, metadata={"gene": "T"})
    assert len(paths) == 2
    loaded, meta = load_candidates_batch(tmp_path / "batch")
    assert len(loaded) == 3
    assert meta["gene"] == "T"


def test_load_batch_empty_dir(tmp_path: Path) -> None:
    """Loading from an empty directory raises FileNotFoundError."""
    d = tmp_path / "empty"
    d.mkdir()
    with pytest.raises(FileNotFoundError, match="No JSON files"):
        load_candidates_batch(d)


def test_deduplicate_first() -> None:
    """Deduplication keeps first occurrence."""
    cands = [
        ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20),
        ASOCandidate("ATGCATGCATGCATGCATGC", 20, 40, 20),  # dup sequence
        ASOCandidate("GCTAGCTAGCTAGCTAGCTA", 40, 60, 20),
    ]
    deduped = deduplicate_candidates(cands, keep="first")
    assert len(deduped) == 2
    assert deduped[0].start == 0  # first occurrence kept


def test_deduplicate_last() -> None:
    """Deduplication with keep='last' keeps last occurrence."""
    cands = [
        ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20),
        ASOCandidate("ATGCATGCATGCATGCATGC", 20, 40, 20),  # dup
        ASOCandidate("GCTAGCTAGCTAGCTAGCTA", 40, 60, 20),
    ]
    deduped = deduplicate_candidates(cands, keep="last")
    assert len(deduped) == 2
    assert deduped[0].start == 20  # last occurrence kept


def test_deduplicate_no_duplicates() -> None:
    """Deduplication with no duplicates returns same count."""
    cands = _make_candidates()
    deduped = deduplicate_candidates(cands)
    assert len(deduped) == 3


def test_deduplicate_invalid_keep() -> None:
    """Invalid keep parameter raises ValueError."""
    cands = _make_candidates()
    with pytest.raises(ValueError, match="keep must be"):
        deduplicate_candidates(cands, keep="middle")


def test_candidates_equal() -> None:
    """candidates_equal compares sequence and position."""
    a = ASOCandidate("ATGC", 0, 4, 4)
    b = ASOCandidate("ATGC", 0, 4, 4)
    c = ASOCandidate("ATGC", 4, 8, 4)
    d = ASOCandidate("GCTA", 0, 4, 4)
    assert candidates_equal(a, b)
    assert not candidates_equal(a, c)  # different position
    assert not candidates_equal(a, d)  # different sequence


def test_validate_roundtrip(tmp_path: Path) -> None:
    """validate_roundtrip returns True for valid candidates."""
    cands = _make_candidates()
    assert validate_roundtrip(cands, tmp_path / "rt.json")
