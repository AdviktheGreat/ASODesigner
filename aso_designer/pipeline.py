"""Pipeline orchestration for ASODesigner.

Ties together all pipeline stages — transcript fetching, candidate
generation, pre-filtering, feature computation (accessibility,
thermodynamics, off-target, variant), scoring, ranking, and reporting
— into a single orchestrated workflow.

Designed from scratch with:
- PipelineConfig for stage-level toggles and parameters
- ASOPipeline class with stage methods and a run() entry point
- PipelineResult for structured output of the full run
- Graceful degradation when optional stages are skipped
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Pipeline configuration
# ---------------------------------------------------------------------------


@dataclass
class PipelineConfig:
    """Configuration for the ASO design pipeline.

    Controls which stages run and their parameters. Built from CLI
    arguments or constructed programmatically.

    Attributes
    ----------
    target : str
        Target gene symbol or transcript accession.
    organism : str
        Organism for NCBI lookup (default "Homo sapiens").
    chemistry : str
        ASO chemistry (default "gapmer").
    config_path : str or None
        Path to a user YAML config (overrides defaults).
    local_fasta : str or None
        Path to a local FASTA file (bypasses NCBI).
    output_dir : str
        Output directory for reports (default "results").
    top_n : int or None
        Keep only the top N candidates (None = keep all).
    min_score : float
        Minimum score threshold (default 0.0).
    run_accessibility : bool
        Compute RNA accessibility (default True).
    run_thermodynamics : bool
        Compute thermodynamic features (default True).
    run_off_target : bool
        Run off-target analysis (default True).
    run_variant : bool
        Run variant/SNP annotation (default True).
    generate_html : bool
        Generate HTML report (default True).
    """

    target: str = ""
    organism: str = "Homo sapiens"
    chemistry: str = "gapmer"
    config_path: str | None = None
    local_fasta: str | None = None
    output_dir: str = "results"
    top_n: int | None = None
    min_score: float = 0.0
    run_accessibility: bool = True
    run_thermodynamics: bool = True
    run_off_target: bool = True
    run_variant: bool = True
    generate_html: bool = True

    @classmethod
    def from_args(cls, args: Any) -> PipelineConfig:
        """Build a PipelineConfig from argparse namespace.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed CLI arguments.

        Returns
        -------
        PipelineConfig
        """
        return cls(
            target=getattr(args, "target", ""),
            organism=getattr(args, "organism", "Homo sapiens"),
            chemistry=getattr(args, "chemistry", "gapmer"),
            config_path=getattr(args, "config", None),
            local_fasta=getattr(args, "local_fasta", None),
            output_dir=getattr(args, "output_dir", "results"),
            top_n=getattr(args, "top_n", None),
            min_score=getattr(args, "min_score", 0.0),
            run_off_target=not getattr(args, "no_offtarget", False),
            run_variant=not getattr(args, "no_variant", False),
            generate_html=not getattr(args, "no_html", False),
        )


# ---------------------------------------------------------------------------
# Pipeline result
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """Structured output of a pipeline run.

    Attributes
    ----------
    target : str
        Target name used.
    transcript_id : str
        Fetched transcript accession.
    transcript_length : int
        Transcript sequence length.
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
    stages_run : list[str]
        Names of stages that were executed.
    errors : list[str]
        Non-fatal errors encountered during the run.
    report_paths : dict[str, str]
        Paths to generated report files.
    """

    target: str = ""
    transcript_id: str = ""
    transcript_length: int = 0
    total_candidates: int = 0
    passed_filters: int = 0
    scored_candidates: int = 0
    ranked_candidates: int = 0
    top_score: float = 0.0
    stages_run: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    report_paths: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "target": self.target,
            "transcript_id": self.transcript_id,
            "transcript_length": self.transcript_length,
            "total_candidates": self.total_candidates,
            "passed_filters": self.passed_filters,
            "scored_candidates": self.scored_candidates,
            "ranked_candidates": self.ranked_candidates,
            "top_score": round(self.top_score, 4),
            "stages_run": self.stages_run,
            "errors": self.errors,
            "report_paths": self.report_paths,
        }

    @property
    def success(self) -> bool:
        """True if the pipeline produced at least one ranked candidate."""
        return self.ranked_candidates > 0
