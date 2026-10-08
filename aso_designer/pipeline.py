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

from aso_designer.config_schema import ASOConfig
from aso_designer.utils.exceptions import ASODesignError
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
    transcriptome_path : str or None
        Path to a transcriptome FASTA for off-target search (overrides
        the off_target.transcriptome_path config option).
    vcf_path : str or None
        Path to a VCF file with known variants (overrides the
        variant.vcf_path config option).
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
    generate_csv : bool
        Also write a CSV table of the ranked candidates (default False).
    generate_fasta : bool
        Also write a FASTA file of the candidate sequences (default False).
    """

    target: str = ""
    organism: str = "Homo sapiens"
    chemistry: str = "gapmer"
    config_path: str | None = None
    local_fasta: str | None = None
    transcriptome_path: str | None = None
    vcf_path: str | None = None
    output_dir: str = "results"
    top_n: int | None = None
    min_score: float = 0.0
    run_accessibility: bool = True
    run_thermodynamics: bool = True
    run_off_target: bool = True
    run_variant: bool = True
    generate_html: bool = True
    generate_csv: bool = False
    generate_fasta: bool = False

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
            transcriptome_path=getattr(args, "transcriptome", None),
            vcf_path=getattr(args, "vcf", None),
            output_dir=getattr(args, "output_dir", "results"),
            top_n=getattr(args, "top_n", None),
            min_score=getattr(args, "min_score", 0.0),
            run_off_target=not getattr(args, "no_offtarget", False),
            run_variant=not getattr(args, "no_variant", False),
            generate_html=not getattr(args, "no_html", False),
            generate_csv=getattr(args, "csv", False),
            generate_fasta=getattr(args, "export_fasta", False),
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
    def aso_config(self) -> ASOConfig:
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
        errors. Expected failures (module errors, missing inputs, I/O
        problems) are logged and recorded; the four feature stages
        (accessibility, thermodynamics, off-target, variant) are
        non-fatal — the pipeline continues and scoring falls back to
        neutral sub-scores for missing features. Programming errors
        (TypeError, AttributeError, KeyError, ...) propagate, since
        they indicate bugs rather than data problems.

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
            except (ASODesignError, RuntimeError, OSError, ValueError) as e:
                # Expected failures: module errors (ViennaRNA/network/config),
                # missing stage inputs, file I/O problems, and validation
                # errors. Programming errors (TypeError, AttributeError, ...)
                # intentionally propagate — they indicate bugs, not data
                # problems.
                msg = f"Stage '{stage_name}' failed: {e}"
                logger.error(msg)
                self.result.errors.append(msg)
                # Feature stages are non-fatal: the pipeline continues and
                # the scoring engine falls back to neutral sub-scores for
                # the missing features.
                if stage_name in ("accessibility", "thermodynamics", "off_target", "variant"):
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
        from aso_designer.fetch_manager import FetchManager
        from aso_designer.utils.sequence import is_accession

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
        """Generate ASO candidates via sliding window.

        Candidates are annotated with transcript regions (5UTR/CDS/3UTR)
        inferred from ORF detection, and each candidate's features dict
        is seeded with the transcript length for positional scoring.

        Uses the memory-managed batched generator so very long
        transcripts (e.g. DMD at ~14 kb) do not require the full
        candidate list to be materialized at once during generation.
        """
        from aso_designer.candidate_generator import generate_candidates_batched
        from aso_designer.quality import extract_metadata

        if self._transcript is None:
            raise RuntimeError("No transcript available for candidate generation")

        metadata = extract_metadata(self._transcript)
        transcript_length = len(self._transcript.sequence)
        self._candidates = []
        for batch in generate_candidates_batched(
            sequence=self._transcript.sequence,
            config=self.aso_config.candidate,
            metadata=metadata,
        ):
            for candidate in batch:
                candidate.features["transcript_length"] = transcript_length
            self._candidates.extend(batch)
        self.result.total_candidates = len(self._candidates)
        logger.info("Generated %d candidates", self.result.total_candidates)

    def _stage_filter(self) -> None:
        """Apply pre-filters to candidates."""
        from aso_designer.filters import run_all_filters

        cand_cfg = self.aso_config.candidate
        result = run_all_filters(
            self._candidates,
            gc_min=cand_cfg.gc_min,
            gc_max=cand_cfg.gc_max,
            homopolymer_max=cand_cfg.homopolymer_max,
        )
        self.result.passed_filters = result["total_passed"]
        logger.info(
            "Passed filters: %d / %d",
            self.result.passed_filters,
            self.result.total_candidates,
        )

    def _stage_accessibility(self) -> None:
        """Compute RNA accessibility features for passing candidates.

        Predicts the accessibility profile of the full target transcript
        with RNAplfold (via TargetStructurePredictor, with disk caching),
        then scores each passing candidate's binding site window with
        ``score_window_accessibility``. The per-candidate score is stored
        in ``candidate.features["accessibility_score"]``.

        Requires ViennaRNA; raises ViennaRNAError if unavailable. The
        pipeline treats this stage as non-fatal and continues with
        neutral accessibility scores.
        """
        from aso_designer.accessibility import (
            TargetStructurePredictor,
            score_window_accessibility,
        )
        from aso_designer.viennarna import verify_viennarna

        if self._transcript is None:
            raise RuntimeError("No transcript available for accessibility analysis")

        acc_cfg = self.aso_config.accessibility
        verify_viennarna()

        logger.info(
            "Predicting accessibility profile (%d nt, window=%d, max_span=%d)",
            len(self._transcript.sequence),
            acc_cfg.window_size,
            acc_cfg.max_pairing_distance,
        )
        predictor = TargetStructurePredictor(
            window_size=acc_cfg.window_size,
            max_pairing_distance=acc_cfg.max_pairing_distance,
        )
        profile = predictor.predict_accessibility(self._transcript.sequence)

        computed = 0
        for candidate in self._candidates:
            if not candidate.passes_filters:
                continue
            # score_window_accessibility uses 1-based inclusive coordinates;
            # candidates use 0-based [start, end).
            score = score_window_accessibility(profile, candidate.start + 1, candidate.end)
            candidate.features["accessibility_score"] = score
            computed += 1

        logger.info("Computed accessibility for %d candidates", computed)

    def _stage_thermodynamics(self) -> None:
        """Compute thermodynamic features for passing candidates.

        Runs ``compute_thermodynamic_summary`` for each passing
        candidate (Tm, duplex/self/target ΔG, net binding energy,
        binding classification) and stores the results in
        ``candidate.features``. Summaries are also collected in
        ``self._thermo_summaries`` for reporting.

        With the default ``method="nearest_neighbor"`` no external
        tools are needed; ViennaRNA methods raise (and the pipeline
        degrades gracefully) when ViennaRNA is unavailable.
        """
        from aso_designer.thermodynamics import (
            compute_thermodynamic_summary,
            set_thermodynamics_cache,
        )

        if self._transcript is None:
            raise RuntimeError("No transcript available for thermodynamic analysis")

        t_cfg = self.aso_config.thermodynamics
        transcript_seq = self._transcript.sequence

        if t_cfg.cache_dir:
            set_thermodynamics_cache(t_cfg.cache_dir)

        computed = 0
        for candidate in self._candidates:
            if not candidate.passes_filters:
                continue
            summary = compute_thermodynamic_summary(
                candidate.sequence,
                target_sequence=transcript_seq[candidate.start : candidate.end],
                chemistry=candidate.chemistry,
                na_conc=t_cfg.na_conc,
                oligo_conc=t_cfg.oligo_conc,
                tm_optimal=t_cfg.tm_optimal,
                tm_min=t_cfg.tm_min,
                tm_max=t_cfg.tm_max,
                self_structure_max=t_cfg.self_structure_max,
                method=t_cfg.method,
                nn_fraction=t_cfg.self_structure_nn_fraction,
            )
            candidate.features.update(
                {
                    "tm": summary.tm,
                    "tm_score": summary.tm_score,
                    "dg_duplex": summary.dg_duplex,
                    "dg_self": summary.dg_self,
                    "dg_target": summary.dg_target,
                    "dg_net": summary.dg_net,
                    "binding_class": summary.binding_class,
                    "self_structure_penalty": summary.self_structure_penalty,
                }
            )
            self._thermo_summaries.append(summary)
            computed += 1

        logger.info("Computed thermodynamics for %d candidates", computed)

    def _stage_off_target(self) -> None:
        """Run transcriptome-wide off-target analysis for passing candidates.

        Builds a k-mer index over the configured transcriptome FASTA
        (``off_target.transcriptome_path`` or the pipeline-level
        ``transcriptome_path``), searches each candidate's binding site
        with the seed-and-extend strategy, and stores the specificity
        score in ``candidate.features["off_target_score"]``. Results are
        also collected in ``self._off_target_results`` for reporting.

        If no transcriptome is configured, the stage logs a warning and
        returns, leaving candidates with neutral off-target scores.
        """
        from aso_designer.off_target import (
            OffTargetSearchEngine,
            TranscriptomeReference,
            score_off_target,
        )

        if self._transcript is None:
            raise RuntimeError("No transcript available for off-target analysis")

        ot_cfg = self.aso_config.off_target
        transcriptome_path = self.config.transcriptome_path or ot_cfg.transcriptome_path

        if not transcriptome_path:
            logger.warning(
                "No transcriptome FASTA configured — skipping off-target "
                "analysis (neutral scores will be used)"
            )
            return

        reference = TranscriptomeReference()
        try:
            reference.load_from_fasta(transcriptome_path)
        except Exception as e:
            logger.warning(
                "Could not load transcriptome %s (%s) — skipping off-target analysis",
                transcriptome_path,
                e,
            )
            return
        if reference.num_transcripts == 0:
            logger.warning(
                "Transcriptome FASTA %s contains no sequences — skipping off-target analysis",
                transcriptome_path,
            )
            return
        logger.info(
            "Built off-target reference: %d transcripts, %d nt total",
            reference.num_transcripts,
            reference.total_length,
        )

        engine = OffTargetSearchEngine(
            reference,
            max_mismatches=ot_cfg.max_mismatches,
            seed_start=ot_cfg.seed_start,
            seed_end=ot_cfg.seed_end,
        )

        scored = 0
        for candidate in self._candidates:
            if not candidate.passes_filters:
                continue
            hits = engine.search(
                candidate.sequence,
                on_target_transcript=self.result.transcript_id,
                on_target_position=candidate.start,
            )
            result, score, classification = score_off_target(
                candidate.sequence,
                hits,
                seed_mismatch_weight=ot_cfg.seed_mismatch_weight,
                max_hits=ot_cfg.max_hits,
            )
            candidate.features["off_target_score"] = score
            candidate.features["off_target_class"] = classification
            self._off_target_results.append(result)
            scored += 1

        logger.info("Scored off-targets for %d candidates", scored)

    def _stage_variant(self) -> None:
        """Run variant/SNP annotation for passing candidates.

        Variant data comes from ``variant.vcf_path`` (or the pipeline-level
        ``vcf_path``) when provided; otherwise variants are fetched from
        the dbSNP/ClinVar APIs using the target's gene name (requires
        network access). Each passing candidate's binding site is
        annotated with overlapping variants and the penalty/risk are
        stored in ``candidate.features``. Annotations are also collected
        in ``self._variant_annotations`` for reporting.

        Variant positions are interpreted in transcript coordinates —
        the same coordinate system as the fetched transcript sequence.
        Genomic-coordinate VCFs require transcript-mapped positions.

        If no variant source is available the stage logs a warning and
        returns, leaving candidates with neutral variant scores.
        """
        from aso_designer.variant import (
            VariantAnnotator,
            VariantError,
            classify_variant_risk,
        )

        if self._transcript is None:
            raise RuntimeError("No transcript available for variant annotation")

        var_cfg = self.aso_config.variant
        ot_cfg = self.aso_config.off_target
        vcf_path = self.config.vcf_path or var_cfg.vcf_path

        variant_set = None
        if vcf_path:
            from aso_designer.variant import load_vcf_filtered

            try:
                variant_set = load_vcf_filtered(vcf_path, min_maf=var_cfg.maf_threshold)
            except VariantError as e:
                logger.warning(
                    "Could not load VCF %s (%s) — skipping variant analysis",
                    vcf_path,
                    e,
                )
                return
        else:
            gene_name = self._transcript.gene_name
            if not gene_name:
                logger.warning(
                    "No VCF path and no gene name for target — skipping variant analysis"
                )
                return
            from aso_designer.variant import VariantFetchManager

            fetch_cfg = self.aso_config.fetch
            manager = VariantFetchManager(
                email=fetch_cfg.email or "",
                api_key=fetch_cfg.api_key,
            )
            try:
                variant_set = manager.fetch_variants(gene_name, organism=self.config.organism)
            except VariantError as e:
                logger.warning(
                    "Variant fetch failed for %s (%s) — skipping variant analysis",
                    gene_name,
                    e,
                )
                return

        if variant_set is None:
            logger.warning("No variant data source available — skipping variant analysis")
            return
        if variant_set.num_variants == 0:
            logger.warning(
                "No variants passed filtering — annotating candidates with an empty variant set"
            )

        # The ASO seed region definition is shared with the off-target stage
        annotator = VariantAnnotator(
            variant_set,
            seed_start=ot_cfg.seed_start,
            seed_end=ot_cfg.seed_end,
            maf_threshold=var_cfg.maf_threshold,
            seed_region_multiplier=var_cfg.seed_region_multiplier,
            flag_pathogenic=var_cfg.flag_pathogenic,
        )

        annotated = 0
        for candidate in self._candidates:
            if not candidate.passes_filters:
                continue
            # Transcript coordinates: pass chrom="" so the annotator
            # matches variant positions directly.
            annotation = annotator.annotate(candidate.start, candidate.end, chrom="")
            candidate.features["variant_penalty"] = annotation.penalty
            candidate.features["variant_risk"] = classify_variant_risk(annotation)
            self._variant_annotations.append(annotation)
            annotated += 1

        logger.info(
            "Annotated %d candidates against %d variants",
            annotated,
            variant_set.num_variants,
        )

    def _stage_score(self) -> None:
        """Score and rank candidates using their computed features.

        Each candidate's ``features`` dict — populated by the
        accessibility, thermodynamics, off-target, and variant stages —
        is passed to the scoring engine. Missing features (skipped
        stages) fall back to neutral per-component defaults. The engine
        is built from the configured scoring weights.
        """
        from aso_designer.scoring import (
            create_score_engine,
            rank_and_filter,
            score_candidates_batch,
        )

        if not self._candidates:
            logger.warning("No candidates to score")
            return

        # Features computed by earlier stages live on each candidate
        passed = [c for c in self._candidates if c.passes_filters]
        features_list = [c.features for c in passed]

        w = self.aso_config.scoring.weights
        engine = create_score_engine(
            weights={
                "accessibility": w.accessibility,
                "thermodynamics": w.thermodynamics,
                "gc_content": w.gc_content,
                "off_target": w.off_target,
                "variant": w.variant,
                "self_structure": w.self_structure,
                "positional": w.positional,
            },
            register_all=True,
        )

        self._scored = score_candidates_batch(passed, features_list, engine=engine)
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
        import aso_designer
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
            metadata={"version": aso_designer.__version__, "stages_run": self.result.stages_run},
            basename=f"{self.config.target}_report",
            generate_html=self.config.generate_html,
            generate_csv=self.config.generate_csv,
            generate_fasta=self.config.generate_fasta,
        )
        self.result.report_paths = {k: str(v) for k, v in paths.items()}
        logger.info("Reports generated: %s", list(paths.keys()))


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------


def run_pipeline(
    target: str,
    chemistry: str = "gapmer",
    config_path: str | None = None,
    local_fasta: str | None = None,
    transcriptome_path: str | None = None,
    vcf_path: str | None = None,
    output_dir: str = "results",
    top_n: int | None = None,
    min_score: float = 0.0,
    run_accessibility: bool = True,
    run_thermodynamics: bool = True,
    run_off_target: bool = True,
    run_variant: bool = True,
    generate_html: bool = True,
    generate_csv: bool = False,
    generate_fasta: bool = False,
    organism: str = "Homo sapiens",
) -> PipelineResult:
    """Run the full ASO design pipeline in one call.

    Parameters
    ----------
    target : str
        Target gene symbol or transcript accession.
    chemistry : str
        ASO chemistry (default "gapmer").
    config_path : str or None
        Path to a user YAML config.
    local_fasta : str or None
        Path to a local FASTA file (bypasses NCBI).
    transcriptome_path : str or None
        Path to a transcriptome FASTA for off-target search.
    vcf_path : str or None
        Path to a VCF file with known variants.
    output_dir : str
        Output directory for reports.
    top_n : int or None
        Keep only the top N candidates.
    min_score : float
        Minimum score threshold.
    run_accessibility : bool
        Run accessibility analysis.
    run_thermodynamics : bool
        Run thermodynamic analysis.
    run_off_target : bool
        Run off-target analysis.
    run_variant : bool
        Run variant annotation.
    generate_html : bool
        Generate HTML report.
    generate_csv : bool
        Also write a CSV table of the ranked candidates.
    generate_fasta : bool
        Also write a FASTA file of the candidate sequences.
    organism : str
        Organism for NCBI lookup.

    Returns
    -------
    PipelineResult
        Structured output of the pipeline run.
    """
    config = PipelineConfig(
        target=target,
        organism=organism,
        chemistry=chemistry,
        config_path=config_path,
        local_fasta=local_fasta,
        transcriptome_path=transcriptome_path,
        vcf_path=vcf_path,
        output_dir=output_dir,
        top_n=top_n,
        min_score=min_score,
        run_accessibility=run_accessibility,
        run_thermodynamics=run_thermodynamics,
        run_off_target=run_off_target,
        run_variant=run_variant,
        generate_html=generate_html,
        generate_csv=generate_csv,
        generate_fasta=generate_fasta,
    )
    pipeline = ASOPipeline(config)
    return pipeline.run()
