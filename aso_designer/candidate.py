"""ASO candidate data structures for ASODesigner.

Defines the core ASOCandidate dataclass that carries a candidate
oligonucleotide through the entire pipeline: from generation, through
filtering, feature computation, scoring, and ranking. Designed from
scratch for ASO-specific needs — GuideX uses a simple gRNACandidate
with sequence + MFE; we track chemistry modifications, region
annotation, a feature vector, and a score vector.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ASOCandidate:
    """A single ASO candidate with all pipeline data.

    Attributes
    ----------
    sequence : str
        The ASO sequence in DNA alphabet (uppercase, ACGT).
    start : int
        Start position in the target transcript (0-based, inclusive).
    end : int
        End position in the target transcript (0-based, exclusive).
    length : int
        Length of the ASO in nucleotides.
    chemistry : str
        Chemistry type: gapmer, mixmer, 2ome, moe, morpholino.
    modification_pattern : list[str]
        Per-nucleotide modification labels (e.g. ["MOE", "MOE", "DNA", ...]).
        "DNA" = unmodified DNA, others = modified nucleotide types.
    backbone_pattern : list[str]
        Per-linkage backbone labels, length = len(sequence) - 1
        ("PS" = phosphorothioate, "PO" = natural phosphodiester).
        Empty until chemistry is applied.
    gc_content : float
        GC content as a fraction (0.0–1.0).
    region : str
        Transcript region where the candidate maps: "5UTR", "CDS",
        "3UTR", "noncoding", or "unknown".
    passes_filters : bool
        Whether the candidate passed all pre-filters (GC, homopolymer,
        end stability, ambiguous bases).
    features : dict[str, Any]
        Computed feature values (accessibility, Tm, dG duplex, etc.).
        Populated by downstream modules (accessibility, thermodynamics,
        off-target, variant).
    scores : dict[str, float]
        Sub-scores and composite score. Keys include "accessibility",
        "thermodynamics", "gc_content", "off_target", "variant",
        "self_structure", "positional", and "composite".
    rank : int | None
        Final rank among all candidates (1 = best). None until scored.
    """

    sequence: str
    start: int
    end: int
    length: int
    chemistry: str = "gapmer"
    modification_pattern: list[str] = field(default_factory=list)
    backbone_pattern: list[str] = field(default_factory=list)
    gc_content: float = 0.0
    region: str = "unknown"
    passes_filters: bool = True
    features: dict[str, Any] = field(default_factory=dict)
    scores: dict[str, float] = field(default_factory=dict)
    rank: int | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a dictionary for JSON storage."""
        return {
            "sequence": self.sequence,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "chemistry": self.chemistry,
            "modification_pattern": self.modification_pattern,
            "backbone_pattern": self.backbone_pattern,
            "gc_content": round(self.gc_content, 4),
            "region": self.region,
            "passes_filters": self.passes_filters,
            "features": self.features,
            "scores": {k: round(v, 4) for k, v in self.scores.items()},
            "rank": self.rank,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ASOCandidate:
        """Deserialize from a dictionary."""
        return cls(
            sequence=data["sequence"],
            start=data["start"],
            end=data["end"],
            length=data["length"],
            chemistry=data.get("chemistry", "gapmer"),
            modification_pattern=data.get("modification_pattern", []),
            backbone_pattern=data.get("backbone_pattern", []),
            gc_content=data.get("gc_content", 0.0),
            region=data.get("region", "unknown"),
            passes_filters=data.get("passes_filters", True),
            features=data.get("features", {}),
            scores=data.get("scores", {}),
            rank=data.get("rank"),
        )

    def __repr__(self) -> str:
        return (
            f"ASOCandidate(seq={self.sequence[:10]}..., "
            f"pos={self.start}-{self.end}, chem={self.chemistry}, "
            f"gc={self.gc_content:.2f}, region={self.region})"
        )
