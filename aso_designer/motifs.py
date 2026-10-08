"""Immunostimulatory-motif and G-quadruplex screening for ASO candidates.

Flags sequence motifs that are associated with innate-immune stimulation
or non-canonical structure formation:

- CpG dinucleotide content: unmethylated CpG motifs can stimulate innate
  immunity through Toll-like receptor 9 (TLR9). The canonical
  immunostimulatory hexamers are TCGTCG (human-type) and GTCGTT
  (mouse-type) (Krieg et al. 1995, Nature 374:546).
- G-quadruplex-prone motifs: four runs of at least three guanines
  separated by short loops can fold into G-quadruplex structures, which
  may affect hybridization behavior.

Scientific note: these counts are sequence-based heuristics used to
surface risk signals alongside candidate scores — they are not
immunogenicity or structure predictions. CpG methylation state and
flanking context, which strongly modulate TLR9 response, are not
modeled.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Canonical TLR9-stimulatory CpG hexamers (heuristic screen).
TLR9_MOTIFS: tuple[str, ...] = ("TCGTCG", "GTCGTT")

# G-quadruplex-prone motif: four G-tracts of >=3 G separated by 1-7 nt.
_G_QUADRUPLEX_RE = re.compile(r"(?:G{3,}[ACGT]{1,7}){3}G{3,}")

# Thresholds for the coarse motif-risk label.
_ELEVATED_CPG_FRACTION = 0.10


@dataclass
class MotifReport:
    """Sequence-motif risk summary for one ASO candidate.

    Attributes
    ----------
    sequence : str
        The analyzed sequence (uppercase).
    cpg_count : int
        Number of CpG (CG) dinucleotides.
    cpg_fraction : float
        CpG count divided by the number of dinucleotide positions
        (0.0 for sequences shorter than 2 nt).
    g_quadruplex_motifs : int
        Number of G-quadruplex-prone motifs found.
    longest_g_run : int
        Length of the longest consecutive guanine run.
    tlr9_stimulatory : bool
        True if any canonical TLR9-stimulatory hexamer is present.
    """

    sequence: str
    cpg_count: int
    cpg_fraction: float
    g_quadruplex_motifs: int
    longest_g_run: int
    tlr9_stimulatory: bool

    def risk_label(self) -> str:
        """Coarse motif-risk label: "elevated", "moderate", or "low".

        - elevated: contains a canonical TLR9 motif or a
          G-quadruplex-prone motif.
        - moderate: CpG-rich (fraction above 0.10) but no flagged motif.
        - low: everything else.
        """
        if self.tlr9_stimulatory or self.g_quadruplex_motifs > 0:
            return "elevated"
        if self.cpg_fraction > _ELEVATED_CPG_FRACTION:
            return "moderate"
        return "low"


def cpg_count(sequence: str) -> int:
    """Count CpG (CG) dinucleotides in the sequence."""
    seq = sequence.upper()
    return sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG")


def cpg_fraction(sequence: str) -> float:
    """CpG count divided by the number of dinucleotide positions.

    Returns 0.0 for sequences shorter than 2 nt.
    """
    n = cpg_count(sequence)
    positions = max(len(sequence) - 1, 0)
    return n / positions if positions else 0.0


def find_g_quadruplex_motifs(sequence: str) -> list[str]:
    """Return the G-quadruplex-prone motifs in the sequence.

    A motif is four runs of at least three guanines separated by 1-7
    nucleotides (the canonical GGG N1-7 GGG N1-7 GGG N1-7 GGG pattern).
    """
    seq = sequence.upper()
    return _G_QUADRUPLEX_RE.findall(seq)


def longest_g_run(sequence: str) -> int:
    """Return the length of the longest consecutive guanine run."""
    longest = 0
    current = 0
    for base in sequence.upper():
        if base == "G":
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def assess_motifs(sequence: str) -> MotifReport:
    """Assess immunostimulatory-motif and G-quadruplex risk for a sequence.

    Parameters
    ----------
    sequence : str
        ASO sequence (DNA alphabet, case-insensitive).

    Returns
    -------
    MotifReport
        Counts and a coarse risk label (see ``MotifReport.risk_label``).
    """
    seq = sequence.upper()
    count = cpg_count(seq)
    report = MotifReport(
        sequence=seq,
        cpg_count=count,
        cpg_fraction=count / max(len(seq) - 1, 1) if len(seq) > 1 else 0.0,
        g_quadruplex_motifs=len(find_g_quadruplex_motifs(seq)),
        longest_g_run=longest_g_run(seq),
        tlr9_stimulatory=any(motif in seq for motif in TLR9_MOTIFS),
    )
    return report
