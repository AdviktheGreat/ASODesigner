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


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


class ASOPipeline:
    """Orchestrates the full ASO design pipeline.

    Runs stages in order: fetch → generate → filter → (accessibility →
    thermodynamics → off-target → variant) → score → rank → report.
    Each stage can be toggled via PipelineConfig.

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration.
    """

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.result = PipelineResult(target=config.target)
        self._config = None  # ASOConfig (loaded lazily)
        self._transcript = None  # TranscriptRecord
        self._candidates: list = []  # list[ASOCandidate]
        self._scored: list = []  # list[CandidateScore]
        self._thermo_summaries: list = []
        self._off_target_results: list = []
        self._variant_annotations: list = []

    @property
    def aso_config(self):
        """Lazily loaded ASOConfig from config system."""
        if self._config is None:
            from aso_designer.config_loader import load_config

            self._config = load_config(
                self.config.config_path,
                chemistry=self.config.chemistry,
            )
        return self._config

    def run(self) -> PipelineResult:
        """Run the full pipeline.

        Executes all enabled stages in order, collecting results and
        errors. Non-fatal errors are logged and recorded but do not
        stop the pipeline.

        Returns
        -------
        PipelineResult
            Structured output of the pipeline run.
        """
        logger.info("Starting ASO pipeline for target: %s", self.config.target)

        stages: list[tuple[str, callable]] = [
            ("fetch", self._stage_fetch),
            ("generate", self._stage_generate),
            ("filter", self._stage_filter),
        ]

        if self.config.run_accessibility:
            stages.append(("accessibility", self._stage_accessibility))
        if self.config.run_thermodynamics:
            stages.append(("thermodynamics", self._stage_thermodynamics))
        if self.config.run_off_target:
            stages.append(("off_target", self._stage_off_target))
        if self.config.run_variant:
            stages.append(("variant", self._stage_variant))

        stages.append(("score", self._stage_score))
        stages.append(("report", self._stage_report))

        for stage_name, stage_fn in stages:
            try:
                logger.info("Running stage: %s", stage_name)
                stage_fn()
                self.result.stages_run.append(stage_name)
            except Exception as e:
                msg = f"Stage '{stage_name}' failed: {e}"
                logger.error(msg)
                self.result.errors.append(msg)
                # Some stages are non-fatal
                if stage_name in ("accessibility", "off_target", "variant"):
                    continue
                break

        logger.info(
            "Pipeline complete: %d candidates, %d ranked, top score %.4f",
            self.result.total_candidates,
            self.result.ranked_candidates,
            self.result.top_score,
        )
        return self.result

    # -----------------------------------------------------------------------
    # Stage methods
    # -----------------------------------------------------------------------

    def _stage_fetch(self) -> None:
        """Fetch the target transcript."""
        from aso_designer.cli import is_accession
        from aso_designer.fetch_manager import FetchManager

        fm = FetchManager(local_fasta=self.config.local_fasta)

        if self.config.local_fasta:
            logger.info("Loading from local FASTA: %s", self.config.local_fasta)
            from aso_designer.local_loader import load_fasta

            records = load_fasta(self.config.local_fasta)
            if records:
                self._transcript = records[0]
        elif is_accession(self.config.target):
            logger.info("Fetching by accession: %s", self.config.target)
            self._transcript = fm.fetch_by_accession(self.config.target)
        else:
            logger.info("Fetching by gene symbol: %s", self.config.target)
            records = fm.fetch_by_gene(self.config.target, self.config.organism)
            if records:
                self._transcript = records[0]

        if self._transcript is None:
            raise RuntimeError(f"Could not fetch transcript for {self.config.target}")

        self.result.transcript_id = self._transcript.accession
        self.result.transcript_length = len(self._transcript.sequence)
        logger.info("Fetched %s (%d bp)", self.result.transcript_id, self.result.transcript_length)

    def _stage_generate(self) -> None:
        """Generate ASO candidates via sliding window."""
        from aso_designer.candidate_generator import generate_candidates

        if self._transcript is None:
            raise RuntimeError("No transcript available for candidate generation")

        self._candidates = generate_candidates(
            sequence=self._transcript.sequence,
            config=self.aso_config.candidate,
        )
        self.result.total_candidates = len(self._candidates)
        logger.info("Generated %d candidates", self.result.total_candidates)

    def _stage_filter(self) -> None:
        """Apply pre-filters to candidates."""
        from aso_designer.filters import run_all_filters

        passed = run_all_filters(self._candidates, self.aso_config)
        self.result.passed_filters = passed
        logger.info("Passed filters: %d / %d", passed, self.result.total_candidates)

    def _stage_accessibility(self) -> None:
        """Compute RNA accessibility features."""
        logger.info("Accessibility stage (placeholder — requires ViennaRNA)")

    def _stage_thermodynamics(self) -> None:
        """Compute thermodynamic features."""
        logger.info("Thermodynamics stage (placeholder — requires ViennaRNA)")

    def _stage_off_target(self) -> None:
        """Run off-target analysis."""
        logger.info("Off-target stage (placeholder — requires transcriptome)")

    def _stage_variant(self) -> None:
        """Run variant/SNP annotation."""
        logger.info("Variant stage (placeholder — requires variant data)")

    def _stage_score(self) -> None:
        """Score and rank candidates."""
        from aso_designer.scoring import rank_and_filter, score_candidates_batch

        if not self._candidates:
            logger.warning("No candidates to score")
            return

        # Build features dicts from candidates that passed filters
        passed = [c for c in self._candidates if c.passes_filters]
        features_list = [{} for _ in passed]

        self._scored = score_candidates_batch(passed, features_list)
        self.result.scored_candidates = len(self._scored)

        ranking = rank_and_filter(
            self._scored,
            min_score=self.config.min_score,
            top_n=self.config.top_n,
        )
        self._scored = ranking.candidates
        self.result.ranked_candidates = ranking.ranked_candidates
        self.result.top_score = ranking.top_score
        logger.info(
            "Scored %d, ranked %d, top score %.4f",
            self.result.scored_candidates,
            self.result.ranked_candidates,
            self.result.top_score,
        )

    def _stage_report(self) -> None:
        """Generate JSON and HTML reports."""
        from aso_designer.reporting import assemble_report

        if not self._scored:
            logger.warning("No scored candidates to report")
            return

        paths = assemble_report(
            scored_candidates=self._scored,
            output_dir=self.config.output_dir,
            target_name=self.config.target,
            target_id=self.result.transcript_id,
            chemistry=self.config.chemistry,
            total_candidates=self.result.total_candidates,
            passed_filters=self.result.passed_filters,
            thermodynamic_summaries=self._thermo_summaries or None,
            off_target_results=self._off_target_results or None,
            variant_annotations=self._variant_annotations or None,
            metadata={"version": "0.0.1", "stages_run": self.result.stages_run},
            basename=f"{self.config.target}_report",
            generate_html=self.config.generate_html,
        )
        self.result.report_paths = {k: str(v) for k, v in paths.items()}
        logger.info("Reports generated: %s", list(paths.keys()))
