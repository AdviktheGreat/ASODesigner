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

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
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


# ---------------------------------------------------------------------------
# JSON report generation
# ---------------------------------------------------------------------------


def generate_json_report(
    report: ReportData,
    output_path: str | Path,
    indent: int = 2,
) -> Path:
    """Write a ReportData to a JSON file.

    Parameters
    ----------
    report : ReportData
        Complete report data to serialize.
    output_path : str or Path
        Output JSON file path. Parent directories are created.
    indent : int
        JSON indentation (default 2).

    Returns
    -------
    Path
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = report.to_dict()

    output_path.write_text(
        json.dumps(data, indent=indent, default=str, ensure_ascii=False),
        encoding="utf-8",
    )

    logger.info("JSON report written to %s", output_path)
    return output_path


def load_json_report(input_path: str | Path) -> ReportData:
    """Load a ReportData from a JSON file.

    Parameters
    ----------
    input_path : str or Path
        Path to a JSON report file.

    Returns
    -------
    ReportData
        Deserialized report data.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Report file not found: {input_path}")

    data = json.loads(input_path.read_text(encoding="utf-8"))
    return ReportData.from_dict(data)


# ---------------------------------------------------------------------------
# Summary statistics computation
# ---------------------------------------------------------------------------


def compute_summary_statistics(
    scored_candidates: list,
    total_candidates: int = 0,
    passed_filters: int = 0,
    target_name: str = "",
    target_id: str = "",
    chemistry: str = "gapmer",
) -> ReportSummary:
    """Compute summary statistics from scored candidates.

    Parameters
    ----------
    scored_candidates : list[CandidateScore]
        Scored and ranked candidates (from rank_candidates).
    total_candidates : int
        Total candidates generated before filtering.
    passed_filters : int
        Candidates that passed pre-filters.
    target_name : str
        Target gene name.
    target_id : str
        Target transcript accession.
    chemistry : str
        Chemistry used for this run.

    Returns
    -------
    ReportSummary
        Computed summary statistics.
    """
    n = len(scored_candidates)

    if n == 0:
        return ReportSummary(
            target_name=target_name,
            target_id=target_id,
            total_candidates=total_candidates,
            passed_filters=passed_filters,
            scored_candidates=0,
            ranked_candidates=0,
            chemistry=chemistry,
            timestamp=datetime.now(timezone.utc).isoformat(),
        )

    scores = sorted(c.final_score for c in scored_candidates)
    top = scores[-1]
    min_s = scores[0]
    mean = sum(scores) / n

    mid = n // 2
    if n % 2 == 0:
        median = (scores[mid - 1] + scores[mid]) / 2
    else:
        median = scores[mid]

    return ReportSummary(
        target_name=target_name,
        target_id=target_id,
        total_candidates=total_candidates if total_candidates else n,
        passed_filters=passed_filters if passed_filters else n,
        scored_candidates=n,
        ranked_candidates=n,
        top_score=top,
        median_score=median,
        mean_score=mean,
        min_score=min_s,
        chemistry=chemistry,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )
