"""Tests for ASO candidate generation, exon mapping, and memory utilities."""

from aso_designer.candidate import ASOCandidate
from aso_designer.candidate_generator import (
    ExonMap,
    annotate_with_exons,
    count_candidates,
    estimate_memory_mb,
    generate_candidates,
    generate_candidates_batched,
    should_use_batching,
)
from aso_designer.config_schema import CandidateConfig
from aso_designer.quality import TranscriptMetadata


def test_generate_candidates_basic() -> None:
    """Basic sliding window generates correct number of candidates."""
    seq = "ATGCATGCATGCATGCATGC"  # 20 nt
    cfg = CandidateConfig(length=10, step=1)
    cands = generate_candidates(seq, cfg)
    assert len(cands) == 11  # 20 - 10 + 1 = 11
    assert all(c.length == 10 for c in cands)


def test_generate_candidates_step() -> None:
    """Step size > 1 reduces candidate count."""
    seq = "ATGCATGCATGCATGCATGC"  # 20 nt
    cfg = CandidateConfig(length=10, step=5)
    cands = generate_candidates(seq, cfg)
    assert len(cands) == 3  # positions 0, 5, 10
    assert cands[0].start == 0
    assert cands[1].start == 5
    assert cands[2].start == 10


def test_generate_candidate_positions() -> None:
    """Candidate start/end positions are correct."""
    seq = "ATGCATGCATGCATGCATGCATGC"  # 24 nt
    cfg = CandidateConfig(length=10, step=1)
    cands = generate_candidates(seq, cfg)
    for i, c in enumerate(cands):
        assert c.start == i
        assert c.end == i + 10
        assert c.sequence == seq[i : i + 10]


def test_generate_candidates_too_short() -> None:
    """Transcript shorter than ASO length yields no candidates."""
    seq = "ATGC"  # 4 nt
    cfg = CandidateConfig(length=20, step=1)
    cands = generate_candidates(seq, cfg)
    assert len(cands) == 0


def test_generate_candidates_region_targeting() -> None:
    """Candidates are generated only within specified regions."""
    seq = "ATGCATGCATGCATGCATGCATGCATGCATGC"  # 32 nt
    cfg = CandidateConfig(length=10, step=1)
    cands = generate_candidates(seq, cfg, regions=[(5, 20)])
    assert len(cands) == 6  # 20 - 5 - 10 + 1 = 6
    assert all(5 <= c.start < 20 for c in cands)
    assert all(c.end <= 20 for c in cands)


def test_generate_candidates_invalid_region() -> None:
    """Invalid regions are skipped with a warning."""
    seq = "ATGCATGCATGCATGC"
    cfg = CandidateConfig(length=10, step=1)
    cands = generate_candidates(seq, cfg, regions=[(-1, 5), (0, 10)])
    assert len(cands) == 1  # only the valid region


def test_generate_candidates_with_metadata() -> None:
    """Candidates get region labels when metadata is provided."""
    seq = "ATG" + "GCT" * 20 + "TAA" + "TTT" * 10  # 3 + 60 + 3 + 30 = 96 nt
    cfg = CandidateConfig(length=10, step=10)
    meta = TranscriptMetadata(
        accession="test",
        gene_name="TEST",
        organism="Homo sapiens",
        length=len(seq),
        region_map=[(0, 3, "5UTR"), (3, 63, "CDS"), (63, 96, "3UTR")],
    )
    cands = generate_candidates(seq, cfg, metadata=meta)
    assert cands[0].region == "5UTR"
    assert cands[1].region == "CDS"
    assert any(c.region == "3UTR" for c in cands)


def test_generate_candidates_skips_invalid_chars() -> None:
    """Candidates with truly invalid characters (X, Z) are skipped."""
    seq = "ATGCXXXXXXATGCATGCATGC"  # contains X (not a valid nucleotide)
    cfg = CandidateConfig(length=10, step=1)
    cands = generate_candidates(seq, cfg)
    for c in cands:
        assert "X" not in c.sequence


def test_generate_candidates_batched() -> None:
    """Batched generation yields same total as non-batched."""
    seq = "ATGCATGCATGCATGC" * 100  # 1600 nt
    cfg = CandidateConfig(length=20, step=1)
    all_cands = generate_candidates(seq, cfg)
    batched = list(generate_candidates_batched(seq, cfg, batch_size=500))
    total_batched = sum(len(b) for b in batched)
    assert total_batched == len(all_cands)
    assert len(batched) > 1  # multiple batches


def test_count_candidates() -> None:
    """count_candidates matches actual generation count."""
    seq = "ATGCATGCATGCATGCATGCATGCATGCATGC"  # 32 nt
    cfg = CandidateConfig(length=15, step=2)
    expected = count_candidates(seq, cfg)
    actual = generate_candidates(seq, cfg)
    assert expected == len(actual)


def test_estimate_memory_mb() -> None:
    """Memory estimation returns positive value."""
    seq = "ATGC" * 1000  # 4000 nt
    cfg = CandidateConfig(length=20, step=1)
    mem = estimate_memory_mb(seq, cfg)
    assert mem > 0


def test_should_use_batching() -> None:
    """Batching is recommended for large transcripts."""
    seq = "ATGC" * 5000  # 20000 nt
    cfg = CandidateConfig(length=20, step=1)
    assert should_use_batching(seq, cfg, memory_threshold_mb=0.001)


# ---------------------------------------------------------------------------
# Exon mapping tests
# ---------------------------------------------------------------------------


def test_exon_map_basic() -> None:
    """ExonMap correctly identifies exon indices."""
    em = ExonMap(exons=[(0, 50), (50, 100), (100, 150)])
    assert em.exon_index(10) == 0
    assert em.exon_index(60) == 1
    assert em.exon_index(120) == 2
    assert em.exon_index(200) is None  # outside all exons


def test_exon_map_spans_junction() -> None:
    """spans_junction detects candidates crossing exon boundaries."""
    em = ExonMap(exons=[(0, 50), (50, 100)])
    c_within = ASOCandidate("ATGC" * 10, 10, 50, 40)  # within exon 0
    c_spanning = ASOCandidate("ATGC" * 10, 40, 80, 40)  # spans junction at 50
    assert not em.spans_junction(c_within)
    assert em.spans_junction(c_spanning)


def test_exon_map_junction_positions() -> None:
    """junction_positions returns all junctions within a candidate."""
    em = ExonMap(exons=[(0, 30), (30, 60), (60, 90)])
    c = ASOCandidate("ATGC" * 20, 25, 105, 80)  # spans junctions at 30 and 60
    junctions = em.junction_positions(c)
    assert 30 in junctions
    assert 60 in junctions


def test_annotate_with_exons() -> None:
    """annotate_with_exons populates features on all candidates."""
    seq = "ATGCATGCATGCATGCATGCATGCATGCATGC"  # 32 nt
    cfg = CandidateConfig(length=10, step=5)
    cands = generate_candidates(seq, cfg)
    em = ExonMap(exons=[(0, 16), (16, 32)])
    annotate_with_exons(cands, em)
    for c in cands:
        assert "exon_index" in c.features
        assert "spans_junction" in c.features
        assert "junction_positions" in c.features
    junction_cands = [c for c in cands if c.features["spans_junction"]]
    assert len(junction_cands) > 0
