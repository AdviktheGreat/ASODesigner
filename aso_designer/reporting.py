"""Report generation for ASODesigner pipeline results.

Generates structured JSON and HTML reports from pipeline outputs:
scored candidates, thermodynamic summaries, off-target results,
variant annotations, and run metadata. Designed from scratch with
pluggable sections, configurable detail levels, and round-trip
JSON serialization.

Inspired by GuideX's reporting approach but redesigned with:
- Structured ReportData container for all pipeline outputs
- Pipeline-level summary statistics (score distribution, top-N)
- Configurable detail levels (summary-only vs full off-target details)
- HTML with embedded CSS (no external dependencies)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Report summary: pipeline-level statistics
# ---------------------------------------------------------------------------


@dataclass
class ReportSummary:
    """Pipeline-level summary statistics for a report.

    Attributes
    ----------
    target_name : str
        Target gene or transcript name.
    target_id : str
        Target transcript accession (e.g., NM_001123.1).
    total_candidates : int
        Total candidates generated.
    passed_filters : int
        Candidates that passed pre-filters.
    scored_candidates : int
        Candidates that were scored.
    ranked_candidates : int
        Candidates after ranking and filtering.
    top_score : float
        Highest final score.
    median_score : float
        Median final score.
    mean_score : float
        Mean final score.
    min_score : float
        Lowest final score.
    chemistry : str
        Chemistry used for this run.
    timestamp : str
        ISO 8601 timestamp of report generation.
    """

    target_name: str = ""
    target_id: str = ""
    total_candidates: int = 0
    passed_filters: int = 0
    scored_candidates: int = 0
    ranked_candidates: int = 0
    top_score: float = 0.0
    median_score: float = 0.0
    mean_score: float = 0.0
    min_score: float = 0.0
    chemistry: str = "gapmer"
    timestamp: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "target_name": self.target_name,
            "target_id": self.target_id,
            "total_candidates": self.total_candidates,
            "passed_filters": self.passed_filters,
            "scored_candidates": self.scored_candidates,
            "ranked_candidates": self.ranked_candidates,
            "top_score": round(self.top_score, 4),
            "median_score": round(self.median_score, 4),
            "mean_score": round(self.mean_score, 4),
            "min_score": round(self.min_score, 4),
            "chemistry": self.chemistry,
            "timestamp": self.timestamp,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ReportSummary:
        """Deserialize from dictionary."""
        return cls(
            target_name=data.get("target_name", ""),
            target_id=data.get("target_id", ""),
            total_candidates=data.get("total_candidates", 0),
            passed_filters=data.get("passed_filters", 0),
            scored_candidates=data.get("scored_candidates", 0),
            ranked_candidates=data.get("ranked_candidates", 0),
            top_score=data.get("top_score", 0.0),
            median_score=data.get("median_score", 0.0),
            mean_score=data.get("mean_score", 0.0),
            min_score=data.get("min_score", 0.0),
            chemistry=data.get("chemistry", "gapmer"),
            timestamp=data.get("timestamp", ""),
        )


# ---------------------------------------------------------------------------
# Report data: full pipeline output container
# ---------------------------------------------------------------------------


@dataclass
class ReportData:
    """Complete report data from a pipeline run.

    Combines summary statistics, ranked candidates, and optional
    detail sections (thermodynamics, off-target, variant) into a
    single serializable structure.

    Attributes
    ----------
    summary : ReportSummary
        Pipeline-level summary statistics.
    candidates : list[dict]
        Ranked candidate dictionaries (from CandidateScore.to_dict()).
    thermodynamic_details : list[dict]
        Per-candidate thermodynamic summaries (optional).
    off_target_details : list[dict]
        Per-candidate off-target results (optional).
    variant_details : list[dict]
        Per-candidate variant annotations (optional).
    metadata : dict[str, Any]
        Additional run metadata (config, parameters, etc.).
    """

    summary: ReportSummary = field(default_factory=ReportSummary)
    candidates: list[dict[str, Any]] = field(default_factory=list)
    thermodynamic_details: list[dict[str, Any]] = field(default_factory=list)
    off_target_details: list[dict[str, Any]] = field(default_factory=list)
    variant_details: list[dict[str, Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a full report dictionary."""
        return {
            "summary": self.summary.to_dict(),
            "candidates": self.candidates,
            "thermodynamic_details": self.thermodynamic_details,
            "off_target_details": self.off_target_details,
            "variant_details": self.variant_details,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ReportData:
        """Deserialize from a report dictionary."""
        return cls(
            summary=ReportSummary.from_dict(data.get("summary", {})),
            candidates=data.get("candidates", []),
            thermodynamic_details=data.get("thermodynamic_details", []),
            off_target_details=data.get("off_target_details", []),
            variant_details=data.get("variant_details", []),
            metadata=data.get("metadata", {}),
        )
