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

import csv
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


# ---------------------------------------------------------------------------
# Build report data from pipeline outputs
# ---------------------------------------------------------------------------


def build_report_data(
    scored_candidates: list,
    total_candidates: int = 0,
    passed_filters: int = 0,
    target_name: str = "",
    target_id: str = "",
    chemistry: str = "gapmer",
    thermodynamic_summaries: list | None = None,
    off_target_results: list | None = None,
    variant_annotations: list | None = None,
    metadata: dict[str, Any] | None = None,
    include_thermo_details: bool = True,
    include_offtarget_details: bool = True,
    include_variant_details: bool = True,
) -> ReportData:
    """Assemble ReportData from pipeline outputs.

    Combines scored candidates, thermodynamic summaries, off-target
    results, and variant annotations into a complete ReportData
    structure ready for JSON or HTML export.

    Parameters
    ----------
    scored_candidates : list[CandidateScore]
        Scored and ranked candidates.
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
    thermodynamic_summaries : list[ThermodynamicSummary] or None
        Per-candidate thermodynamic summaries.
    off_target_results : list or None
        Per-candidate off-target results.
    variant_annotations : list or None
        Per-candidate variant annotations.
    metadata : dict or None
        Additional run metadata.
    include_thermo_details : bool
        Include thermodynamic detail section.
    include_offtarget_details : bool
        Include off-target detail section.
    include_variant_details : bool
        Include variant detail section.

    Returns
    -------
    ReportData
        Complete report data.
    """
    summary = compute_summary_statistics(
        scored_candidates,
        total_candidates=total_candidates,
        passed_filters=passed_filters,
        target_name=target_name,
        target_id=target_id,
        chemistry=chemistry,
    )

    candidates = [c.to_dict() for c in scored_candidates]

    thermo_details: list[dict[str, Any]] = []
    if include_thermo_details and thermodynamic_summaries:
        for ts in thermodynamic_summaries:
            if hasattr(ts, "to_dict"):
                thermo_details.append(ts.to_dict())
            else:
                thermo_details.append(_thermo_summary_to_dict(ts))

    off_target_details: list[dict[str, Any]] = []
    if include_offtarget_details and off_target_results:
        for ot in off_target_results:
            off_target_details.append(_off_target_to_dict(ot))

    variant_details: list[dict[str, Any]] = []
    if include_variant_details and variant_annotations:
        for va in variant_annotations:
            variant_details.append(_variant_annotation_to_dict(va))

    return ReportData(
        summary=summary,
        candidates=candidates,
        thermodynamic_details=thermo_details,
        off_target_details=off_target_details,
        variant_details=variant_details,
        metadata=metadata or {},
    )


def _thermo_summary_to_dict(ts: Any) -> dict[str, Any]:
    """Convert a ThermodynamicSummary to a dict."""
    return {
        "tm": round(ts.tm, 2),
        "tm_score": round(ts.tm_score, 4),
        "dg_duplex": round(ts.dg_duplex, 2),
        "dg_self": round(ts.dg_self, 2),
        "self_structure_penalty": round(ts.self_structure_penalty, 4),
        "dg_target": round(ts.dg_target, 2),
        "dg_net": round(ts.dg_net, 2),
        "binding_class": ts.binding_class,
    }


def _off_target_to_dict(ot: Any) -> dict[str, Any]:
    """Convert an off-target result to a dict."""
    if isinstance(ot, tuple):
        result, score, specificity = ot
        return {
            "off_target_score": round(score, 4),
            "specificity_class": specificity,
            "total_hits": getattr(result, "total_hits", 0),
            "off_target_hits": getattr(result, "off_target_hits", 0),
            "perfect_hits": getattr(result, "perfect_hits", 0),
        }
    return {
        "off_target_score": 0.0,
        "total_hits": getattr(ot, "total_hits", 0),
        "off_target_hits": getattr(ot, "off_target_hits", 0),
    }


def _variant_annotation_to_dict(va: Any) -> dict[str, Any]:
    """Convert a variant annotation to a dict."""
    return {
        "penalty": round(getattr(va, "penalty", 0.0), 4),
        "overlapping_variants": getattr(va, "overlapping_variants", 0),
        "seed_variants": getattr(va, "seed_variants", 0),
        "non_seed_variants": getattr(va, "non_seed_variants", 0),
        "has_pathogenic": getattr(va, "has_pathogenic", False),
    }


# ---------------------------------------------------------------------------
# HTML report generation
# ---------------------------------------------------------------------------

_HTML_CSS = """\
/* ASODesigner report styles */
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  background: #f5f5f5; color: #333; line-height: 1.6;
}
.container { max-width: 1200px; margin: 0 auto; padding: 20px; }
h1 { color: #1a1a2e; margin-bottom: 8px; }
h2 { color: #16213e; margin: 20px 0 10px; border-bottom: 2px solid #0f3460; padding-bottom: 4px; }
h3 { color: #0f3460; margin: 15px 0 8px; }
.subtitle { color: #666; font-size: 0.9em; margin-bottom: 20px; }
.summary-grid {
  display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 12px; margin-bottom: 20px;
}
.summary-card {
  background: #fff; border-radius: 8px; padding: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}
.summary-card .label { font-size: 0.8em; color: #666; text-transform: uppercase; }
.summary-card .value { font-size: 1.4em; font-weight: bold; color: #1a1a2e; }
table {
  width: 100%; border-collapse: collapse; background: #fff;
  border-radius: 8px; overflow: hidden;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin-bottom: 20px;
}
th { background: #0f3460; color: #fff; padding: 10px 12px; text-align: left; font-size: 0.85em; }
td { padding: 8px 12px; border-bottom: 1px solid #eee; font-size: 0.9em; }
tr:hover { background: #f9f9f9; }
.score-bar {
  display: inline-block; height: 8px; border-radius: 4px;
  background: #e0e0e0; position: relative; width: 80px; vertical-align: middle;
}
.score-bar .fill { height: 100%; border-radius: 4px; }
.score-high { background: #4caf50; }
.score-mid { background: #ff9800; }
.score-low { background: #f44336; }
.badge {
  display: inline-block; padding: 2px 8px; border-radius: 12px;
  font-size: 0.75em; font-weight: bold;
}
.badge-strong { background: #c8e6c9; color: #2e7d32; }
.badge-moderate { background: #fff9c4; color: #f57f17; }
.badge-weak { background: #ffcdd2; color: #c62828; }
.detail-section { background: #fff; border-radius: 8px; padding: 16px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
.footer { text-align: center; color: #999; font-size: 0.8em; margin-top: 30px; }
"""


def _score_color(score: float) -> str:
    """Return CSS class for a score value."""
    if score >= 0.7:
        return "score-high"
    elif score >= 0.4:
        return "score-mid"
    return "score-low"


def _score_bar_html(score: float) -> str:
    """Return HTML for a score bar."""
    pct = int(score * 100)
    color = _score_color(score)
    return (
        f'<span class="score-bar">'
        f'<span class="fill {color}" style="width:{pct}%"></span>'
        f"</span> {score:.3f}"
    )


def _binding_badge(binding_class: str) -> str:
    """Return HTML badge for binding classification."""
    cls = {
        "strong": "badge-strong",
        "moderate": "badge-moderate",
        "weak": "badge-weak",
    }.get(binding_class, "")
    return f'<span class="badge {cls}">{binding_class or "N/A"}</span>'


def generate_html_report(
    report: ReportData,
    output_path: str | Path,
    title: str = "ASODesigner Report",
) -> Path:
    """Generate a self-contained HTML report from ReportData.

    Produces an HTML file with embedded CSS (no external dependencies),
    summary statistics cards, a candidate table with score bars, and
    optional detail sections.

    Parameters
    ----------
    report : ReportData
        Complete report data.
    output_path : str or Path
        Output HTML file path.
    title : str
        HTML page title.

    Returns
    -------
    Path
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    s = report.summary

    html_parts: list[str] = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="UTF-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
        f"<title>{title}</title>",
        f"<style>{_HTML_CSS}</style>",
        "</head>",
        "<body>",
        '<div class="container">',
        f"<h1>{title}</h1>",
        f'<p class="subtitle">Target: {s.target_name} ({s.target_id}) &mdash; '
        f"Chemistry: {s.chemistry} &mdash; Generated: {s.timestamp}</p>",
    ]

    # Summary section
    html_parts.append("<h2>Summary</h2>")
    html_parts.append('<div class="summary-grid">')
    summary_items = [
        ("Total Candidates", str(s.total_candidates)),
        ("Passed Filters", str(s.passed_filters)),
        ("Scored", str(s.scored_candidates)),
        ("Ranked", str(s.ranked_candidates)),
        ("Top Score", f"{s.top_score:.4f}"),
        ("Median Score", f"{s.median_score:.4f}"),
        ("Mean Score", f"{s.mean_score:.4f}"),
        ("Min Score", f"{s.min_score:.4f}"),
    ]
    for label, value in summary_items:
        html_parts.append(
            f'<div class="summary-card">'
            f'<div class="label">{label}</div>'
            f'<div class="value">{value}</div>'
            f"</div>"
        )
    html_parts.append("</div>")

    # Score distribution chart
    if report.candidates:
        html_parts.append("<h2>Score Distribution</h2>")
        html_parts.append(_generate_score_distribution_chart(report.candidates))

    # Candidate table
    html_parts.append("<h2>Ranked Candidates</h2>")
    if report.candidates:
        html_parts.append(_generate_candidate_table(report.candidates))
    else:
        html_parts.append("<p>No candidates to display.</p>")

    # Detail sections
    if report.thermodynamic_details:
        html_parts.append("<h2>Thermodynamic Details</h2>")
        html_parts.append(_generate_thermo_table(report.thermodynamic_details))

    if report.off_target_details:
        html_parts.append("<h2>Off-Target Details</h2>")
        html_parts.append(_generate_offtarget_table(report.off_target_details))

    if report.variant_details:
        html_parts.append("<h2>Variant Details</h2>")
        html_parts.append(_generate_variant_table(report.variant_details))

    # Metadata
    if report.metadata:
        html_parts.append("<h2>Run Metadata</h2>")
        html_parts.append('<div class="detail-section">')
        html_parts.append("<table><tr><th>Key</th><th>Value</th></tr>")
        for k, v in report.metadata.items():
            html_parts.append(f"<tr><td>{k}</td><td>{v}</td></tr>")
        html_parts.append("</table>")
        html_parts.append("</div>")

    # Footer
    html_parts.append(
        '<p class="footer">Generated by ASODesigner &mdash; '
        "in-silico ASO candidate design pipeline</p>"
    )
    html_parts.append("</div>")
    html_parts.append("</body>")
    html_parts.append("</html>")

    html = "\n".join(html_parts)
    output_path.write_text(html, encoding="utf-8")

    logger.info("HTML report written to %s", output_path)
    return output_path


def _generate_candidate_table(candidates: list[dict[str, Any]]) -> str:
    """Generate the candidate table HTML."""
    parts: list[str] = [
        "<table>",
        "<tr><th>Rank</th><th>Sequence</th><th>Start</th><th>End</th>"
        "<th>Chemistry</th><th>Final Score</th>"
        "<th>Accessibility</th><th>Thermo</th><th>GC</th>"
        "<th>Off-Target</th><th>Variant</th><th>Self-Struct</th>"
        "<th>Position</th></tr>",
    ]
    for c in candidates:
        rank = c.get("rank", "-")
        seq = c.get("sequence", "")
        start = c.get("start", 0)
        end = c.get("end", 0)
        chem = c.get("chemistry", "")
        final = c.get("final_score", 0.0)
        comp = c.get("component_scores", {})

        def _comp_score(key: str, comp_dict: dict = comp) -> str:
            v = comp_dict.get(key)
            return _score_bar_html(v) if v is not None else "N/A"

        parts.append(
            f"<tr><td>{rank}</td><td><code>{seq}</code></td>"
            f"<td>{start}</td><td>{end}</td><td>{chem}</td>"
            f"<td>{_score_bar_html(final)}</td>"
            f"<td>{_comp_score('accessibility')}</td>"
            f"<td>{_comp_score('thermodynamics')}</td>"
            f"<td>{_comp_score('gc_content')}</td>"
            f"<td>{_comp_score('off_target')}</td>"
            f"<td>{_comp_score('variant')}</td>"
            f"<td>{_comp_score('self_structure')}</td>"
            f"<td>{_comp_score('positional')}</td></tr>"
        )
    parts.append("</table>")
    return "\n".join(parts)


def _generate_thermo_table(details: list[dict[str, Any]]) -> str:
    """Generate the thermodynamic details table HTML."""
    parts: list[str] = [
        '<div class="detail-section">',
        "<table><tr><th>Tm (°C)</th><th>Tm Score</th><th>ΔG Duplex</th>"
        "<th>ΔG Self</th><th>Self Penalty</th><th>ΔG Target</th>"
        "<th>ΔG Net</th><th>Binding</th></tr>",
    ]
    for d in details:
        binding = _binding_badge(d.get("binding_class", ""))
        parts.append(
            f"<tr><td>{d.get('tm', 0):.1f}</td>"
            f"<td>{_score_bar_html(d.get('tm_score', 0))}</td>"
            f"<td>{d.get('dg_duplex', 0):.2f}</td>"
            f"<td>{d.get('dg_self', 0):.2f}</td>"
            f"<td>{d.get('self_structure_penalty', 0):.4f}</td>"
            f"<td>{d.get('dg_target', 0):.2f}</td>"
            f"<td>{d.get('dg_net', 0):.2f}</td>"
            f"<td>{binding}</td></tr>"
        )
    parts.append("</table>")
    parts.append("</div>")
    return "\n".join(parts)


def _generate_offtarget_table(details: list[dict[str, Any]]) -> str:
    """Generate the off-target details table HTML."""
    parts: list[str] = [
        '<div class="detail-section">',
        "<table><tr><th>Off-Target Score</th><th>Specificity</th>"
        "<th>Total Hits</th><th>Off-Target Hits</th>"
        "<th>Perfect Hits</th></tr>",
    ]
    for d in details:
        score = d.get("off_target_score", 0.0)
        parts.append(
            f"<tr><td>{_score_bar_html(score)}</td>"
            f"<td>{d.get('specificity_class', 'N/A')}</td>"
            f"<td>{d.get('total_hits', 0)}</td>"
            f"<td>{d.get('off_target_hits', 0)}</td>"
            f"<td>{d.get('perfect_hits', 0)}</td></tr>"
        )
    parts.append("</table>")
    parts.append("</div>")
    return "\n".join(parts)


def _generate_variant_table(details: list[dict[str, Any]]) -> str:
    """Generate the variant details table HTML."""
    parts: list[str] = [
        '<div class="detail-section">',
        "<table><tr><th>Penalty</th><th>Overlapping</th>"
        "<th>Seed Variants</th><th>Non-Seed Variants</th>"
        "<th>Pathogenic</th></tr>",
    ]
    for d in details:
        pathogenic = "Yes" if d.get("has_pathogenic") else "No"
        parts.append(
            f"<tr><td>{d.get('penalty', 0):.4f}</td>"
            f"<td>{d.get('overlapping_variants', 0)}</td>"
            f"<td>{d.get('seed_variants', 0)}</td>"
            f"<td>{d.get('non_seed_variants', 0)}</td>"
            f"<td>{pathogenic}</td></tr>"
        )
    parts.append("</table>")
    parts.append("</div>")
    return "\n".join(parts)


def _generate_score_distribution_chart(
    candidates: list[dict[str, Any]],
    n_bins: int = 10,
) -> str:
    """Generate an SVG histogram of final score distribution.

    Parameters
    ----------
    candidates : list[dict]
        Candidate dicts with 'final_score' key.
    n_bins : int
        Number of histogram bins (default 10).

    Returns
    -------
    str
        HTML string containing an inline SVG chart.
    """
    scores = [c.get("final_score", 0.0) for c in candidates]
    if not scores:
        return "<p>No scores to display.</p>"

    # Build histogram bins
    bin_width = 1.0 / n_bins
    bins = [0] * n_bins
    for s in scores:
        idx = min(int(s / bin_width), n_bins - 1)
        bins[idx] += 1

    max_count = max(bins) if bins else 1
    if max_count == 0:
        max_count = 1

    # SVG dimensions
    chart_w = 600
    chart_h = 200
    bar_w = chart_w / n_bins
    bar_gap = 4
    bar_actual_w = bar_w - bar_gap
    label_h = 20
    svg_h = chart_h + label_h + 10

    parts: list[str] = [
        f'<svg width="{chart_w}" height="{svg_h}" '
        f'viewBox="0 0 {chart_w} {svg_h}" '
        'style="background:#fff;border-radius:8px;'
        'box-shadow:0 1px 3px rgba(0,0,0,0.1)">',
        # Y-axis label
        '<text x="8" y="15" font-size="10" fill="#666">Count</text>',
        # X-axis label
        f'<text x="{chart_w - 60}" y="{svg_h - 4}" font-size="10" fill="#666">Score</text>',
    ]

    for i, count in enumerate(bins):
        x = i * bar_w + bar_gap / 2
        bar_h = (count / max_count) * (chart_h - 10) if count > 0 else 0
        y = chart_h - bar_h
        color = _score_color((i + 0.5) / n_bins)
        color_hex = {
            "score-high": "#4caf50",
            "score-mid": "#ff9800",
            "score-low": "#f44336",
        }.get(color, "#4caf50")

        parts.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" '
            f'width="{bar_actual_w:.1f}" height="{bar_h:.1f}" '
            f'fill="{color_hex}" rx="2" />'
        )
        if count > 0:
            parts.append(
                f'<text x="{x + bar_actual_w / 2:.1f}" y="{y - 4:.1f}" '
                f'text-anchor="middle" font-size="9" fill="#333">{count}</text>'
            )
        # X-axis labels (score ranges)
        lo = i * bin_width
        hi = (i + 1) * bin_width
        parts.append(
            f'<text x="{x + bar_actual_w / 2:.1f}" y="{chart_h + 14}" '
            f'text-anchor="middle" font-size="8" fill="#666">'
            f"{lo:.1f}-{hi:.1f}</text>"
        )

    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

# Fixed core columns for candidate CSV exports. Component-score columns
# are appended dynamically from the candidates' component_scores maps.
_CSV_CORE_COLUMNS = (
    "rank",
    "sequence",
    "start",
    "end",
    "chemistry",
    "composite_score",
    "chemistry_bonus",
    "final_score",
)


def export_candidates_csv(candidates: list[dict[str, Any]], output_path: str | Path) -> Path:
    """Export ranked candidates as a CSV table.

    Each row is one candidate. In addition to the fixed core columns,
    one column per scoring component is emitted (union across rows,
    sorted alphabetically); candidates missing a component get an empty
    cell.

    Parameters
    ----------
    candidates : list[dict]
        Ranked candidate dictionaries (``CandidateScore.to_dict()`` rows,
        as stored in ``ReportData.candidates``).
    output_path : str or Path
        Destination CSV file. Parent directories are created.

    Returns
    -------
    Path
        The path written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    component_columns: list[str] = []
    seen: set[str] = set()
    for row in candidates:
        for key in row.get("component_scores", {}):
            if key not in seen:
                seen.add(key)
                component_columns.append(key)
    component_columns.sort()

    columns = list(_CSV_CORE_COLUMNS) + component_columns

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for row in candidates:
            components = row.get("component_scores", {})
            record = []
            for col in columns:
                if col in _CSV_CORE_COLUMNS:
                    value = row.get(col, "")
                else:
                    value = components.get(col, "")
                record.append(value)
            writer.writerow(record)

    logger.info("Exported %d candidates to %s", len(candidates), output_path)
    return output_path


def export_candidates_fasta(candidates: list[dict[str, Any]], output_path: str | Path) -> Path:
    """Export candidates as FASTA (for oligo ordering / downstream tools).

    Each record's header carries the candidate index (guaranteed unique)
    plus position, chemistry, rank, and score when available:

        >aso_1 rank=1 start=10 end=30 chemistry=gapmer score=0.7521
        ATGCATGCATGCATGCATGC

    Parameters
    ----------
    candidates : list[dict]
        Ranked candidate dictionaries (``CandidateScore.to_dict()`` rows,
        as stored in ``ReportData.candidates``).
    output_path : str or Path
        Destination FASTA file. Parent directories are created.

    Returns
    -------
    Path
        The path written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    for i, row in enumerate(candidates, start=1):
        fields = [f"aso_{i}"]
        if row.get("rank") is not None:
            fields.append(f"rank={row['rank']}")
        fields.append(f"start={row.get('start', '')}")
        fields.append(f"end={row.get('end', '')}")
        fields.append(f"chemistry={row.get('chemistry', '')}")
        if row.get("final_score") is not None:
            fields.append(f"score={row['final_score']:.4f}")
        lines.append(">".join(["", " ".join(fields)]))
        lines.append(row.get("sequence", ""))

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n" if lines else "")

    logger.info("Exported %d candidates to %s", len(candidates), output_path)
    return output_path


# ---------------------------------------------------------------------------
# Combined report generation
# ---------------------------------------------------------------------------


def generate_report(
    report: ReportData,
    output_dir: str | Path,
    basename: str = "aso_report",
    generate_json: bool = True,
    generate_html: bool = True,
    generate_csv: bool = False,
    generate_fasta: bool = False,
    title: str = "ASODesigner Report",
) -> dict[str, Path]:
    """Generate reports from ReportData.

    Parameters
    ----------
    report : ReportData
        Complete report data.
    output_dir : str or Path
        Directory for output files.
    basename : str
        Base filename (without extension).
    generate_json : bool
        Write JSON report (default True).
    generate_html : bool
        Write HTML report (default True).
    generate_csv : bool
        Write a CSV table of the ranked candidates (default False).
    generate_fasta : bool
        Write a FASTA file of the candidate sequences for oligo
        ordering (default False).
    title : str
        HTML page title.

    Returns
    -------
    dict[str, Path]
        Maps format ('json', 'html', 'csv', 'fasta') to output file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}

    if generate_json:
        json_path = output_dir / f"{basename}.json"
        paths["json"] = generate_json_report(report, json_path)

    if generate_html:
        html_path = output_dir / f"{basename}.html"
        paths["html"] = generate_html_report(report, html_path, title=title)

    if generate_csv:
        csv_path = output_dir / f"{basename}.csv"
        paths["csv"] = export_candidates_csv(report.candidates, csv_path)

    if generate_fasta:
        fasta_path = output_dir / f"{basename}_candidates.fasta"
        paths["fasta"] = export_candidates_fasta(report.candidates, fasta_path)

    logger.info("Reports generated in %s: %s", output_dir, list(paths.keys()))
    return paths


def generate_report_from_ranking(
    ranking: Any,
    output_dir: str | Path,
    target_name: str = "",
    target_id: str = "",
    chemistry: str = "gapmer",
    total_candidates: int = 0,
    passed_filters: int = 0,
    basename: str = "aso_report",
    generate_json: bool = True,
    generate_html: bool = True,
    generate_csv: bool = False,
    generate_fasta: bool = False,
    thermodynamic_summaries: list | None = None,
    off_target_results: list | None = None,
    variant_annotations: list | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Path]:
    """Generate reports directly from a RankingResult.

    Convenience function that builds ReportData from a RankingResult
    (from score_rank_and_filter or rank_and_filter) and writes
    JSON + HTML reports.

    Parameters
    ----------
    ranking : RankingResult
        Ranking result with scored candidates.
    output_dir : str or Path
        Output directory.
    target_name : str
        Target gene name.
    target_id : str
        Target transcript accession.
    chemistry : str
        Chemistry used.
    total_candidates : int
        Total candidates before filtering.
    passed_filters : int
        Candidates that passed pre-filters.
    basename : str
        Base filename.
    generate_json : bool
        Write JSON report.
    generate_html : bool
        Write HTML report.
    thermodynamic_summaries : list or None
        Per-candidate thermodynamic summaries.
    off_target_results : list or None
        Per-candidate off-target results.
    variant_annotations : list or None
        Per-candidate variant annotations.
    metadata : dict or None
        Additional run metadata.

    Returns
    -------
    dict[str, Path]
        Maps format to output file path.
    """
    report = build_report_data(
        scored_candidates=ranking.candidates,
        total_candidates=total_candidates or ranking.total_candidates,
        passed_filters=passed_filters,
        target_name=target_name,
        target_id=target_id,
        chemistry=chemistry,
        thermodynamic_summaries=thermodynamic_summaries,
        off_target_results=off_target_results,
        variant_annotations=variant_annotations,
        metadata=metadata,
    )

    return generate_report(
        report,
        output_dir,
        basename=basename,
        generate_json=generate_json,
        generate_html=generate_html,
        generate_csv=generate_csv,
        generate_fasta=generate_fasta,
    )


# ---------------------------------------------------------------------------
# Report assembler: orchestrates pipeline output collection
# ---------------------------------------------------------------------------


class ReportAssembler:
    """Collects pipeline results and assembles a ReportData.

    Provides a fluent interface for incrementally adding results from
    each pipeline stage (candidates, scoring, thermodynamics, off-target,
    variant) and then assembling a complete report.

    Parameters
    ----------
    target_name : str
        Target gene name.
    target_id : str
        Target transcript accession.
    chemistry : str
        Chemistry used for this run.
    """

    def __init__(
        self,
        target_name: str = "",
        target_id: str = "",
        chemistry: str = "gapmer",
    ):
        self.target_name = target_name
        self.target_id = target_id
        self.chemistry = chemistry
        self._total_candidates = 0
        self._passed_filters = 0
        self._scored_candidates: list = []
        self._thermo_summaries: list | None = None
        self._off_target_results: list | None = None
        self._variant_annotations: list | None = None
        self._metadata: dict[str, Any] = {}

    def set_candidate_counts(
        self,
        total: int,
        passed: int = 0,
    ) -> ReportAssembler:
        """Set candidate generation counts."""
        self._total_candidates = total
        self._passed_filters = passed
        return self

    def set_scored_candidates(self, candidates: list) -> ReportAssembler:
        """Set scored and ranked candidates."""
        self._scored_candidates = candidates
        return self

    def set_thermodynamic_summaries(
        self,
        summaries: list,
    ) -> ReportAssembler:
        """Set per-candidate thermodynamic summaries."""
        self._thermo_summaries = summaries
        return self

    def set_off_target_results(
        self,
        results: list,
    ) -> ReportAssembler:
        """Set per-candidate off-target results."""
        self._off_target_results = results
        return self

    def set_variant_annotations(
        self,
        annotations: list,
    ) -> ReportAssembler:
        """Set per-candidate variant annotations."""
        self._variant_annotations = annotations
        return self

    def add_metadata(self, key: str, value: Any) -> ReportAssembler:
        """Add a metadata key-value pair."""
        self._metadata[key] = value
        return self

    def set_metadata(self, metadata: dict[str, Any]) -> ReportAssembler:
        """Replace all metadata."""
        self._metadata = dict(metadata)
        return self

    def assemble(
        self,
        include_thermo: bool = True,
        include_offtarget: bool = True,
        include_variant: bool = True,
    ) -> ReportData:
        """Assemble all collected data into a ReportData.

        Parameters
        ----------
        include_thermo : bool
            Include thermodynamic detail section.
        include_offtarget : bool
            Include off-target detail section.
        include_variant : bool
            Include variant detail section.

        Returns
        -------
        ReportData
            Complete report data.
        """
        return build_report_data(
            scored_candidates=self._scored_candidates,
            total_candidates=self._total_candidates,
            passed_filters=self._passed_filters,
            target_name=self.target_name,
            target_id=self.target_id,
            chemistry=self.chemistry,
            thermodynamic_summaries=self._thermo_summaries,
            off_target_results=self._off_target_results,
            variant_annotations=self._variant_annotations,
            metadata=self._metadata,
            include_thermo_details=include_thermo,
            include_offtarget_details=include_offtarget,
            include_variant_details=include_variant,
        )


# ---------------------------------------------------------------------------
# Top-level assemble_report convenience
# ---------------------------------------------------------------------------


def assemble_report(
    scored_candidates: list,
    output_dir: str | Path,
    target_name: str = "",
    target_id: str = "",
    chemistry: str = "gapmer",
    total_candidates: int = 0,
    passed_filters: int = 0,
    thermodynamic_summaries: list | None = None,
    off_target_results: list | None = None,
    variant_annotations: list | None = None,
    metadata: dict[str, Any] | None = None,
    basename: str = "aso_report",
    generate_json: bool = True,
    generate_html: bool = True,
    generate_csv: bool = False,
    generate_fasta: bool = False,
    include_thermo: bool = True,
    include_offtarget: bool = True,
    include_variant: bool = True,
) -> dict[str, Path]:
    """Assemble and write reports from all pipeline outputs in one call.

    This is the top-level convenience function that chains
    build_report_data → generate_report.

    Parameters
    ----------
    scored_candidates : list[CandidateScore]
        Scored and ranked candidates.
    output_dir : str or Path
        Output directory for report files.
    target_name : str
        Target gene name.
    target_id : str
        Target transcript accession.
    chemistry : str
        Chemistry used.
    total_candidates : int
        Total candidates before filtering.
    passed_filters : int
        Candidates that passed pre-filters.
    thermodynamic_summaries : list or None
        Per-candidate thermodynamic summaries.
    off_target_results : list or None
        Per-candidate off-target results.
    variant_annotations : list or None
        Per-candidate variant annotations.
    metadata : dict or None
        Additional run metadata.
    basename : str
        Base filename for reports.
    generate_json : bool
        Write JSON report.
    generate_html : bool
        Write HTML report.
    include_thermo : bool
        Include thermodynamic details.
    include_offtarget : bool
        Include off-target details.
    include_variant : bool
        Include variant details.

    Returns
    -------
    dict[str, Path]
        Maps format to output file path.
    """
    report = build_report_data(
        scored_candidates=scored_candidates,
        total_candidates=total_candidates,
        passed_filters=passed_filters,
        target_name=target_name,
        target_id=target_id,
        chemistry=chemistry,
        thermodynamic_summaries=thermodynamic_summaries,
        off_target_results=off_target_results,
        variant_annotations=variant_annotations,
        metadata=metadata,
        include_thermo_details=include_thermo,
        include_offtarget_details=include_offtarget,
        include_variant_details=include_variant,
    )

    return generate_report(
        report,
        output_dir,
        basename=basename,
        generate_json=generate_json,
        generate_html=generate_html,
        generate_csv=generate_csv,
        generate_fasta=generate_fasta,
    )


def assemble_report_from_config(
    scored_candidates: list,
    config: Any,
    target_name: str = "",
    target_id: str = "",
    chemistry: str = "gapmer",
    total_candidates: int = 0,
    passed_filters: int = 0,
    thermodynamic_summaries: list | None = None,
    off_target_results: list | None = None,
    variant_annotations: list | None = None,
    metadata: dict[str, Any] | None = None,
    basename: str = "aso_report",
) -> dict[str, Path]:
    """Assemble and write reports using ReportingConfig settings.

    Reads the ReportingConfig section of the ASO config to determine
    output directory, whether to generate HTML, and whether to include
    off-target details.

    Parameters
    ----------
    scored_candidates : list[CandidateScore]
        Scored and ranked candidates.
    config : ASOConfig
        Loaded ASO configuration (from config_loader.load_config).
    target_name : str
        Target gene name.
    target_id : str
        Target transcript accession.
    chemistry : str
        Chemistry used.
    total_candidates : int
        Total candidates before filtering.
    passed_filters : int
        Candidates that passed pre-filters.
    thermodynamic_summaries : list or None
        Per-candidate thermodynamic summaries.
    off_target_results : list or None
        Per-candidate off-target results.
    variant_annotations : list or None
        Per-candidate variant annotations.
    metadata : dict or None
        Additional run metadata.
    basename : str
        Base filename for reports.

    Returns
    -------
    dict[str, Path]
        Maps format to output file path.
    """
    reporting_cfg = config.reporting

    return assemble_report(
        scored_candidates=scored_candidates,
        output_dir=reporting_cfg.output_dir,
        target_name=target_name,
        target_id=target_id,
        chemistry=chemistry,
        total_candidates=total_candidates,
        passed_filters=passed_filters,
        thermodynamic_summaries=thermodynamic_summaries,
        off_target_results=off_target_results,
        variant_annotations=variant_annotations,
        metadata=metadata,
        basename=basename,
        generate_json=True,
        generate_html=reporting_cfg.generate_html,
        include_thermo=True,
        include_offtarget=reporting_cfg.include_offtarget_details,
        include_variant=True,
    )


# ---------------------------------------------------------------------------
# Batch report generation
# ---------------------------------------------------------------------------


def generate_reports_batch(
    targets: list[dict[str, Any]],
    output_dir: str | Path,
    generate_json: bool = True,
    generate_html: bool = True,
) -> list[dict[str, Path]]:
    """Generate reports for multiple targets in batch.

    Each target dict should contain the same keys as assemble_report
    keyword arguments. Reports are written to subdirectories named
    after the target.

    Parameters
    ----------
    targets : list[dict]
        List of target dicts. Each dict must have 'scored_candidates'
        and 'target_name'. Optional keys: 'target_id', 'chemistry',
        'total_candidates', 'passed_filters', 'thermodynamic_summaries',
        'off_target_results', 'variant_annotations', 'metadata'.
    output_dir : str or Path
        Base output directory. Each target gets a subdirectory.
    generate_json : bool
        Write JSON reports.
    generate_html : bool
        Write HTML reports.

    Returns
    -------
    list[dict[str, Path]]
        List of path dicts, one per target.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_paths: list[dict[str, Path]] = []

    for target in targets:
        name = target.get("target_name", "unknown")
        scored = target.get("scored_candidates", [])
        if not scored:
            logger.warning("Skipping target '%s': no scored candidates", name)
            continue

        target_dir = output_dir / name
        basename = f"{name}_report"

        paths = assemble_report(
            scored_candidates=scored,
            output_dir=target_dir,
            target_name=name,
            target_id=target.get("target_id", ""),
            chemistry=target.get("chemistry", "gapmer"),
            total_candidates=target.get("total_candidates", 0),
            passed_filters=target.get("passed_filters", 0),
            thermodynamic_summaries=target.get("thermodynamic_summaries"),
            off_target_results=target.get("off_target_results"),
            variant_annotations=target.get("variant_annotations"),
            metadata=target.get("metadata"),
            basename=basename,
            generate_json=generate_json,
            generate_html=generate_html,
        )
        all_paths.append(paths)

    logger.info("Batch reports generated: %d targets in %s", len(all_paths), output_dir)
    return all_paths
