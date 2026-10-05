"""Multi-factor scoring and ranking for ASO candidates.

Combines sub-scores from accessibility, thermodynamics, GC content,
off-target analysis, variant awareness, self-structure, and positional
preference into a weighted composite score. Candidates are ranked by
composite score and filtered by configurable thresholds.

Inspired by GuideX's scoring framework but redesigned with a
pluggable component architecture, chemistry-aware bonuses, and
transparent weight documentation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from aso_designer.utils.exceptions import ScoringError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Score context: carries all computed features for a candidate
# ---------------------------------------------------------------------------


@dataclass
class ScoreContext:
    """All computed features for a single ASO candidate.

    Populated by upstream modules (accessibility, thermodynamics,
    off-target, variant) and consumed by score components.

    Attributes
    ----------
    sequence : str
        ASO sequence.
    chemistry : str
        Chemistry type.
    start : int
        Start position on the transcript (0-based).
    end : int
        End position on the transcript (0-based, exclusive).
    length : int
        ASO length.
    gc_content : float
        GC content fraction (0.0–1.0).
    region : str
        Transcript region ("5UTR", "CDS", "3UTR", etc.).

    Accessibility features:
    accessibility_score : float
        Windowed accessibility score (0.0–1.0).

    Thermodynamic features:
    tm : float
        Melting temperature in °C.
    tm_score : float
        Tm score (0.0–1.0, 1.0 = optimal).
    dg_duplex : float
        Duplex formation ΔG (kcal/mol).
    dg_self : float
        Self-structure ΔG (kcal/mol).
    dg_target : float
        Target structure ΔG (kcal/mol).
    dg_net : float
        Net binding ΔG (kcal/mol).
    binding_class : str
        Binding strength classification.

    Off-target features:
    off_target_score : float
        Off-target specificity score (0.0–1.0).

    Variant features:
    variant_penalty : float
        Variant penalty (0.0 = no penalty).
    variant_risk : str
        Variant risk classification.

    Self-structure features:
    self_structure_penalty : float
        Self-structure penalty (0.0 = no penalty).

    Positional features:
    transcript_length : int
        Total transcript length.
    """

    sequence: str = ""
    chemistry: str = "gapmer"
    start: int = 0
    end: int = 0
    length: int = 0
    gc_content: float = 0.0
    region: str = "unknown"

    # Accessibility
    accessibility_score: float = -1.0

    # Thermodynamics
    tm: float = -999.0
    tm_score: float = -1.0
    dg_duplex: float = 0.0
    dg_self: float = 0.0
    dg_target: float = 0.0
    dg_net: float = 0.0
    binding_class: str = ""

    # Off-target
    off_target_score: float = -1.0

    # Variant
    variant_penalty: float = 0.0
    variant_risk: str = ""

    # Self-structure
    self_structure_penalty: float = 0.0

    # Positional
    transcript_length: int = 0


# ---------------------------------------------------------------------------
# Score result: output of a single score component
# ---------------------------------------------------------------------------


@dataclass
class ScoreResult:
    """Result of a single score component.

    Attributes
    ----------
    component_name : str
        Name of the scoring component.
    raw_score : float
        Raw score before weighting (0.0–1.0).
    weight : float
        Weight applied to this component.
    weighted_score : float
        Raw score × weight.
    """

    component_name: str
    raw_score: float
    weight: float
    weighted_score: float


# ---------------------------------------------------------------------------
# Score component protocol
# ---------------------------------------------------------------------------


@runtime_checkable
class ScoreComponent(Protocol):
    """Protocol for pluggable score components.

    Each component computes a sub-score in [0.0, 1.0] from a
    ScoreContext. Components are registered with the ScoreEngine
    and combined into a composite score.
    """

    @property
    def name(self) -> str:
        """Component name (e.g., 'accessibility', 'thermodynamics')."""
        ...

    def score(self, context: ScoreContext) -> float:
        """Compute a sub-score in [0.0, 1.0].

        Parameters
        ----------
        context : ScoreContext
            Computed features for the candidate.

        Returns
        -------
        float
            Sub-score in [0.0, 1.0].
        """
        ...


# ---------------------------------------------------------------------------
# Composite score result
# ---------------------------------------------------------------------------


@dataclass
class CompositeScore:
    """Complete scoring result for a single candidate.

    Attributes
    ----------
    component_results : list[ScoreResult]
        Individual component results.
    composite_score : float
        Weighted sum of all component scores.
    chemistry_bonus : float
        Chemistry-specific bonus added to the composite.
    final_score : float
        Composite + chemistry bonus, clamped to [0.0, 1.0].
    """

    component_results: list[ScoreResult]
    composite_score: float
    chemistry_bonus: float
    final_score: float


# ---------------------------------------------------------------------------
# Score engine
# ---------------------------------------------------------------------------


class ScoreEngine:
    """Orchestrates multi-factor scoring for ASO candidates.

    Registers score components with their weights, computes
    sub-scores, and combines them into a weighted composite score
    with an optional chemistry bonus.

    Parameters
    ----------
    weights : dict[str, float]
        Maps component names to their weights. Weights should
        sum to 1.0 (validated by config schema).
    chemistry_bonus : float
        Bonus added to the composite score (default 0.0).
    """

    def __init__(
        self,
        weights: dict[str, float] | None = None,
        chemistry_bonus: float = 0.0,
    ):
        self._weights = weights or _default_weights()
        self.chemistry_bonus = chemistry_bonus
        self._components: dict[str, ScoreComponent] = {}

    @property
    def weights(self) -> dict[str, float]:
        """Registered component weights."""
        return self._weights

    def register_component(self, component: ScoreComponent) -> None:
        """Register a score component.

        Parameters
        ----------
        component : ScoreComponent
            Component implementing the ScoreComponent protocol.

        Raises
        ------
        ScoringError
            If the component name is not in the weights dict.
        """
        if component.name not in self._weights:
            raise ScoringError(
                f"Component '{component.name}' has no weight configured. "
                f"Known weights: {list(self._weights.keys())}"
            )
        self._components[component.name] = component
        logger.debug("Registered score component: %s", component.name)

    def score(self, context: ScoreContext) -> CompositeScore:
        """Compute the composite score for a candidate.

        Parameters
        ----------
        context : ScoreContext
            Computed features for the candidate.

        Returns
        -------
        CompositeScore
            Complete scoring result with all component scores.

        Raises
        ------
        ScoringError
            If no components are registered.
        """
        if not self._components:
            raise ScoringError("No score components registered")

        results: list[ScoreResult] = []
        composite = 0.0

        for name, component in self._components.items():
            weight = self._weights.get(name, 0.0)
            raw = component.score(context)
            raw = max(0.0, min(1.0, raw))  # Clamp to [0, 1]
            weighted = raw * weight
            composite += weighted
            results.append(
                ScoreResult(
                    component_name=name,
                    raw_score=raw,
                    weight=weight,
                    weighted_score=weighted,
                )
            )

        final = max(0.0, min(1.0, composite + self.chemistry_bonus))

        return CompositeScore(
            component_results=results,
            composite_score=composite,
            chemistry_bonus=self.chemistry_bonus,
            final_score=final,
        )

    def score_batch(self, contexts: list[ScoreContext]) -> list[CompositeScore]:
        """Score multiple candidates.

        Parameters
        ----------
        contexts : list[ScoreContext]
            List of candidate contexts.

        Returns
        -------
        list[CompositeScore]
            Scoring results for each candidate.
        """
        return [self.score(ctx) for ctx in contexts]


def _default_weights() -> dict[str, float]:
    """Return default scoring weights from config."""
    return {
        "accessibility": 0.25,
        "thermodynamics": 0.20,
        "gc_content": 0.10,
        "off_target": 0.20,
        "variant": 0.10,
        "self_structure": 0.10,
        "positional": 0.05,
    }


# ---------------------------------------------------------------------------
# Candidate score: combines composite score with candidate metadata
# ---------------------------------------------------------------------------


@dataclass
class CandidateScore:
    """A scored ASO candidate with ranking metadata.

    Attributes
    ----------
    sequence : str
        ASO sequence.
    start : int
        Start position on the transcript (0-based).
    end : int
        End position on the transcript (0-based, exclusive).
    chemistry : str
        Chemistry type.
    composite_score : CompositeScore
        Full scoring breakdown.
    final_score : float
        Final score (composite + chemistry bonus), [0.0, 1.0].
    rank : int | None
        Rank among all candidates (1 = best). None until ranked.
    """

    sequence: str
    start: int
    end: int
    chemistry: str
    composite_score: CompositeScore
    final_score: float
    rank: int | None = None

    @property
    def component_scores(self) -> dict[str, float]:
        """Map component name → raw score."""
        return {r.component_name: r.raw_score for r in self.composite_score.component_results}

    def to_dict(self) -> dict:
        """Serialize to dictionary for reporting."""
        return {
            "sequence": self.sequence,
            "start": self.start,
            "end": self.end,
            "chemistry": self.chemistry,
            "final_score": round(self.final_score, 4),
            "rank": self.rank,
            "component_scores": {
                r.component_name: round(r.raw_score, 4)
                for r in self.composite_score.component_results
            },
            "composite_score": round(self.composite_score.composite_score, 4),
            "chemistry_bonus": self.composite_score.chemistry_bonus,
        }


# ---------------------------------------------------------------------------
# Engine factory
# ---------------------------------------------------------------------------


def create_score_engine(
    weights: dict[str, float] | None = None,
    chemistry_bonus: float = 0.0,
    register_all: bool = True,
) -> ScoreEngine:
    """Create a ScoreEngine with default or custom weights.

    Parameters
    ----------
    weights : dict[str, float] or None
        Custom weights. If None, loads from default config.
    chemistry_bonus : float
        Chemistry-specific bonus (default 0.0).
    register_all : bool
        If True, register all available score components
        (accessibility, thermodynamics, gc_content, off_target,
        variant, self_structure, positional). Each component class is
        defined in this module; a component that fails to construct is
        skipped with a debug log.

    Returns
    -------
    ScoreEngine
        Configured score engine with registered components.
    """
    if weights is None:
        weights = _default_weights()

    engine = ScoreEngine(weights=weights, chemistry_bonus=chemistry_bonus)

    if register_all:
        _register_default_components(engine)

    return engine


def _register_default_components(engine: ScoreEngine) -> None:
    """Register all default score components with the engine.

    Components are referenced directly (no dynamic imports): they live
    in this module and are fully defined by the time this function is
    called at runtime. A component that raises a ScoringError during
    construction is skipped with a debug log so one broken component
    cannot prevent the engine from starting.
    """
    for component_cls in (
        AccessibilityScoreComponent,
        ThermodynamicsScoreComponent,
        GCContentScoreComponent,
        OffTargetScoreComponent,
        VariantScoreComponent,
        SelfStructureScoreComponent,
        PositionalScoreComponent,
    ):
        try:
            engine.register_component(component_cls())
        except ScoringError as e:
            logger.debug("Skipping component %s: %s", component_cls.__name__, e)


# ---------------------------------------------------------------------------
# Core score components
# ---------------------------------------------------------------------------


class AccessibilityScoreComponent:
    """Scores candidates by target RNA accessibility.

    Uses the pre-computed accessibility score from the accessibility
    module (RNAplfold-based windowed unpaired probability aggregation).
    A score of 1.0 means the target site is fully accessible (unpaired).

    If accessibility was not computed, returns a neutral 0.5 score.
    """

    @property
    def name(self) -> str:
        """Return 'accessibility'."""
        return "accessibility"

    def score(self, context: ScoreContext) -> float:
        """Compute accessibility sub-score from pre-computed accessibility_score."""
        if context.accessibility_score < 0:
            logger.debug("Accessibility not computed; using neutral 0.5")
            return 0.5
        return context.accessibility_score


class ThermodynamicsScoreComponent:
    """Scores candidates by thermodynamic profile.

    Combines Tm score and binding class into a single score:
    - Tm score (0–1) weighted by ``tm_weight`` (default 0.6)
    - Binding class score weighted by ``1 - tm_weight`` (default 0.4)

    If thermodynamics were not computed, returns a neutral 0.5 score.

    Rationale and tunability
    ------------------------
    The two sub-signals capture complementary properties: the Tm score
    measures absolute duplex stability against the configured Tm window,
    while the binding class reflects net energetics after subtracting
    self-structure and target-opening costs (see ``classify_binding``).
    The 0.6/0.4 split and the binding-class score mapping
    (strong=1.0, moderate=0.5, weak=0.0) are heuristic defaults chosen
    so that well-designed 18-22 nt gapmers score high but imperfect
    binders remain distinguishable; they are not constants from a single
    publication. Both are exposed as constructor parameters for tuning.

    Parameters
    ----------
    tm_weight : float
        Weight for the Tm sub-score in [0, 1]; the binding-class
        sub-score receives the remainder (default 0.6).
    binding_class_scores : dict[str, float] or None
        Maps binding class ("strong", "moderate", "weak") to a sub-score
        in [0, 1]. Unknown classes receive ``unknown_class_score``
        (default 0.5). If None, uses the default mapping.
    unknown_class_score : float
        Score for an unrecognized binding class (default 0.5).
    """

    _DEFAULT_BINDING_CLASS_SCORES = {"strong": 1.0, "moderate": 0.5, "weak": 0.0}

    def __init__(
        self,
        tm_weight: float = 0.6,
        binding_class_scores: dict[str, float] | None = None,
        unknown_class_score: float = 0.5,
    ):
        if not 0.0 <= tm_weight <= 1.0:
            raise ValueError(f"tm_weight must be in [0, 1], got {tm_weight}")
        self.tm_weight = tm_weight
        self.binding_class_scores = (
            dict(binding_class_scores)
            if binding_class_scores is not None
            else dict(self._DEFAULT_BINDING_CLASS_SCORES)
        )
        self.unknown_class_score = unknown_class_score

    @property
    def name(self) -> str:
        """Return 'thermodynamics'."""
        return "thermodynamics"

    def score(self, context: ScoreContext) -> float:
        """Compute thermodynamics sub-score from Tm score and binding class."""
        if context.tm_score < 0:
            logger.debug("Thermodynamics not computed; using neutral 0.5")
            return 0.5

        tm_component = context.tm_score * self.tm_weight
        binding_score = self.binding_class_scores.get(
            context.binding_class, self.unknown_class_score
        )
        binding_component = binding_score * (1.0 - self.tm_weight)

        return tm_component + binding_component


class GCContentScoreComponent:
    """Scores candidates by GC content.

    Optimal GC range is 30–60%. Scores:
    - 1.0 if GC is in [0.3, 0.6]
    - Linear decrease to 0.0 at GC=0.0 or GC=1.0

    Parameters
    ----------
    gc_min : float
        Minimum optimal GC fraction (default 0.30).
    gc_max : float
        Maximum optimal GC fraction (default 0.60).
    """

    def __init__(self, gc_min: float = 0.30, gc_max: float = 0.60):
        self.gc_min = gc_min
        self.gc_max = gc_max

    @property
    def name(self) -> str:
        """Return 'gc_content'."""
        return "gc_content"

    def score(self, context: ScoreContext) -> float:
        """Compute GC content sub-score based on optimal range."""
        gc = context.gc_content

        if self.gc_min <= gc <= self.gc_max:
            return 1.0

        if gc < self.gc_min:
            # Linear decrease from gc_min (1.0) to 0.0 (0.0)
            if self.gc_min == 0:
                return 0.0
            return gc / self.gc_min
        else:
            # Linear decrease from gc_max (1.0) to 1.0 (0.0)
            if self.gc_max == 1.0:
                return 0.0
            return (1.0 - gc) / (1.0 - self.gc_max)


# ---------------------------------------------------------------------------
# Convenience functions
# ---------------------------------------------------------------------------


def score_candidate(
    candidate,
    features: dict | None = None,
    engine: ScoreEngine | None = None,
) -> CandidateScore:
    """Score a single ASO candidate using its computed features.

    Builds a ScoreContext from the candidate and its features dict,
    then runs the score engine to produce a CandidateScore.

    Parameters
    ----------
    candidate : ASOCandidate
        The candidate to score (provides sequence, start, end, chemistry,
        gc_content, region).
    features : dict or None
        Computed features (accessibility_score, tm, tm_score, dg_duplex,
        dg_self, dg_target, dg_net, binding_class, off_target_score,
        variant_penalty, variant_risk, self_structure_penalty,
        transcript_length). Missing features use neutral defaults.
    engine : ScoreEngine or None
        Score engine to use. If None, creates a default engine.

    Returns
    -------
    CandidateScore
        Scored candidate with full breakdown.
    """
    if engine is None:
        engine = create_score_engine()
    if features is None:
        features = {}

    ctx = ScoreContext(
        sequence=candidate.sequence,
        chemistry=candidate.chemistry,
        start=candidate.start,
        end=candidate.end,
        length=candidate.length,
        gc_content=candidate.gc_content,
        region=candidate.region,
        accessibility_score=features.get("accessibility_score", -1.0),
        tm=features.get("tm", -999.0),
        tm_score=features.get("tm_score", -1.0),
        dg_duplex=features.get("dg_duplex", 0.0),
        dg_self=features.get("dg_self", 0.0),
        dg_target=features.get("dg_target", 0.0),
        dg_net=features.get("dg_net", 0.0),
        binding_class=features.get("binding_class", ""),
        off_target_score=features.get("off_target_score", -1.0),
        variant_penalty=features.get("variant_penalty", 0.0),
        variant_risk=features.get("variant_risk", ""),
        self_structure_penalty=features.get("self_structure_penalty", 0.0),
        transcript_length=features.get("transcript_length", 0),
    )

    composite = engine.score(ctx)

    return CandidateScore(
        sequence=candidate.sequence,
        start=candidate.start,
        end=candidate.end,
        chemistry=candidate.chemistry,
        composite_score=composite,
        final_score=composite.final_score,
    )


def score_candidates_batch(
    candidates: list,
    features_list: list[dict] | None = None,
    engine: ScoreEngine | None = None,
) -> list[CandidateScore]:
    """Score multiple ASO candidates in batch.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        List of candidates to score.
    features_list : list[dict] or None
        List of feature dicts (same length as candidates). If None,
        all candidates get empty feature dicts.
    engine : ScoreEngine or None
        Score engine to use. If None, creates a default engine.

    Returns
    -------
    list[CandidateScore]
        Scored candidates in the same order as input.
    """
    if engine is None:
        engine = create_score_engine()
    if features_list is None:
        features_list = [{}] * len(candidates)

    if len(features_list) != len(candidates):
        raise ScoringError(
            f"features_list length {len(features_list)} != candidates length {len(candidates)}"
        )

    return [
        score_candidate(cand, features, engine)
        for cand, features in zip(candidates, features_list, strict=True)
    ]


# ---------------------------------------------------------------------------
# Additional score components
# ---------------------------------------------------------------------------


class OffTargetScoreComponent:
    """Scores candidates by off-target specificity.

    Uses the pre-computed off-target specificity score from the
    off-target module. A score of 1.0 means highly specific (no
    off-targets), 0.0 means non-specific.

    If off-target analysis was not computed, returns a neutral 0.5.
    """

    @property
    def name(self) -> str:
        """Return 'off_target'."""
        return "off_target"

    def score(self, context: ScoreContext) -> float:
        """Compute off-target sub-score from pre-computed off_target_score."""
        if context.off_target_score < 0:
            logger.debug("Off-target not computed; using neutral 0.5")
            return 0.5
        return context.off_target_score


class VariantScoreComponent:
    """Scores candidates by variant/SNP overlap penalty.

    Converts the variant penalty (0.0 = no variants, higher = worse)
    to a score (1.0 = no variants, 0.0 = heavy penalty).

    If no variant data is available, returns 1.0 (no penalty).
    """

    @property
    def name(self) -> str:
        """Return 'variant'."""
        return "variant"

    def score(self, context: ScoreContext) -> float:
        """Compute variant sub-score from variant penalty (1 - penalty)."""
        if not context.variant_risk:
            # No variant data computed; assume no penalty
            return 1.0
        return max(0.0, 1.0 - context.variant_penalty)


class SelfStructureScoreComponent:
    """Scores candidates by self-structure penalty.

    Converts the self-structure penalty (0.0 = no self-structure,
    higher = worse) to a score (1.0 = no self-structure, 0.0 = heavy).

    If no self-structure data is available, returns 1.0.
    """

    @property
    def name(self) -> str:
        """Return 'self_structure'."""
        return "self_structure"

    def score(self, context: ScoreContext) -> float:
        """Compute self-structure sub-score from penalty (1 - penalty)."""
        if context.self_structure_penalty == 0.0 and context.dg_self == 0.0:
            # No self-structure data computed; assume no penalty
            return 1.0
        return max(0.0, 1.0 - context.self_structure_penalty)


class PositionalScoreComponent:
    """Scores candidates by transcript region and position.

    Region preferences (ASO design best practices):
    - CDS: 1.0 (most effective for RNase H recruitment)
    - 5UTR: 0.7 (good for translational blocking)
    - 3UTR: 0.5 (moderate effectiveness)
    - noncoding: 0.6 (useful for non-coding RNA targets)
    - unknown: 0.3

    Also applies a mild position preference: candidates near the
    center of the transcript get a small bonus (±10%).

    Parameters
    ----------
    region_weights : dict[str, float]
        Maps region names to base scores.
    position_bonus : float
        Maximum position-based bonus (default 0.1).
    """

    def __init__(
        self,
        region_weights: dict[str, float] | None = None,
        position_bonus: float = 0.1,
    ):
        self.region_weights = region_weights or {
            "CDS": 1.0,
            "5UTR": 0.7,
            "3UTR": 0.5,
            "noncoding": 0.6,
            "unknown": 0.3,
        }
        self.position_bonus = position_bonus

    @property
    def name(self) -> str:
        """Return 'positional'."""
        return "positional"

    def score(self, context: ScoreContext) -> float:
        """Compute positional sub-score from region weight and center bonus."""
        base = self.region_weights.get(context.region, 0.3)

        # Position bonus: candidates near center get a small boost
        if context.transcript_length > 0 and context.length > 0:
            center = context.transcript_length / 2
            candidate_center = (context.start + context.end) / 2
            # Distance from center, normalized to [0, 1]
            max_dist = context.transcript_length / 2
            if max_dist > 0:
                dist = abs(candidate_center - center) / max_dist
                # Bonus is highest at center (dist=0), zero at ends (dist=1)
                bonus = self.position_bonus * (1.0 - dist)
                base = min(1.0, base + bonus)

        return base


# ---------------------------------------------------------------------------
# Config-driven engine creation
# ---------------------------------------------------------------------------


def create_engine_from_config(
    config_path: str | None = None,
    chemistry: str | None = None,
) -> ScoreEngine:
    """Create a ScoreEngine from YAML configuration.

    Loads weights and chemistry bonus from the config system,
    then registers all available score components.

    Parameters
    ----------
    config_path : str or None
        Path to a user config YAML. If None, uses default config.
    chemistry : str or None
        Chemistry name to apply preset adjustments.

    Returns
    -------
    ScoreEngine
        Configured engine with all available components registered.
    """
    from aso_designer.config_loader import load_config

    config = load_config(config_path, chemistry=chemistry)

    weights = {
        "accessibility": config.scoring.weights.accessibility,
        "thermodynamics": config.scoring.weights.thermodynamics,
        "gc_content": config.scoring.weights.gc_content,
        "off_target": config.scoring.weights.off_target,
        "variant": config.scoring.weights.variant,
        "self_structure": config.scoring.weights.self_structure,
        "positional": config.scoring.weights.positional,
    }

    chemistry_bonus = 0.0
    # Check for chemistry bonus stored by config_loader
    raw_config = config.model_dump()
    if "scoring" in raw_config and "chemistry_bonus" in raw_config.get("scoring", {}):
        chemistry_bonus = raw_config["scoring"]["chemistry_bonus"]

    engine = create_score_engine(
        weights=weights,
        chemistry_bonus=chemistry_bonus,
        register_all=True,
    )

    logger.info(
        "Created engine from config: %d components, chemistry_bonus=%.4f",
        len(engine._components),
        chemistry_bonus,
    )
    return engine


# ---------------------------------------------------------------------------
# Integration helpers: build ScoreContext from upstream results
# ---------------------------------------------------------------------------


def build_features_dict(
    thermo_summary=None,
    off_target_result=None,
    variant_annotation=None,
    accessibility_score: float = -1.0,
    transcript_length: int = 0,
) -> dict:
    """Build a features dict from upstream module results.

    Combines results from thermodynamics, off-target, variant, and
    accessibility modules into a single dict for score_candidate.

    Parameters
    ----------
    thermo_summary : ThermodynamicSummary or None
        Result from compute_thermodynamic_summary().
    off_target_result : OffTargetResult or None
        Result from score_off_target() (the tuple's first element).
    variant_annotation : VariantAnnotation or None
        Result from VariantAnnotator.annotate().
    accessibility_score : float
        Pre-computed accessibility score (default -1.0 = not computed).
    transcript_length : int
        Total transcript length for positional scoring.

    Returns
    -------
    dict
        Features dict ready for score_candidate().
    """
    features: dict = {}

    if thermo_summary is not None:
        features["tm"] = thermo_summary.tm
        features["tm_score"] = thermo_summary.tm_score
        features["dg_duplex"] = thermo_summary.dg_duplex
        features["dg_self"] = thermo_summary.dg_self
        features["dg_target"] = thermo_summary.dg_target
        features["dg_net"] = thermo_summary.dg_net
        features["binding_class"] = thermo_summary.binding_class
        features["self_structure_penalty"] = thermo_summary.self_structure_penalty

    if off_target_result is not None:
        if isinstance(off_target_result, tuple):
            _, ot_score, _ = off_target_result
            features["off_target_score"] = ot_score
        else:
            from aso_designer.off_target import calculate_off_target_score

            features["off_target_score"] = calculate_off_target_score(off_target_result)

    if variant_annotation is not None:
        from aso_designer.variant import classify_variant_risk

        features["variant_penalty"] = variant_annotation.penalty
        features["variant_risk"] = classify_variant_risk(variant_annotation)

    features["accessibility_score"] = accessibility_score
    features["transcript_length"] = transcript_length

    return features


# ---------------------------------------------------------------------------
# Final ranking
# ---------------------------------------------------------------------------


def rank_candidates(
    scored: list[CandidateScore],
    descending: bool = True,
) -> list[CandidateScore]:
    """Rank scored candidates by final score.

    Assigns rank (1 = best) to each candidate and returns them
    sorted by score.

    Parameters
    ----------
    scored : list[CandidateScore]
        Scored candidates (from score_candidates_batch).
    descending : bool
        If True, highest score = rank 1 (default). If False,
        lowest score = rank 1.

    Returns
    -------
    list[CandidateScore]
        Sorted and ranked candidates.
    """
    sorted_scored = sorted(
        scored,
        key=lambda c: c.final_score,
        reverse=descending,
    )
    for i, cand in enumerate(sorted_scored):
        cand.rank = i + 1
    return sorted_scored


def filter_by_score(
    scored: list[CandidateScore],
    min_score: float = 0.0,
    max_score: float = 1.0,
    top_n: int | None = None,
) -> list[CandidateScore]:
    """Filter scored candidates by score thresholds.

    Parameters
    ----------
    scored : list[CandidateScore]
        Scored candidates.
    min_score : float
        Minimum final score to keep (default 0.0).
    max_score : float
        Maximum final score to keep (default 1.0).
    top_n : int or None
        If set, keep only the top N candidates by score.

    Returns
    -------
    list[CandidateScore]
        Filtered candidates (not re-ranked).
    """
    filtered = [c for c in scored if min_score <= c.final_score <= max_score]
    if top_n is not None and top_n > 0:
        filtered = sorted(filtered, key=lambda c: c.final_score, reverse=True)
        filtered = filtered[:top_n]
    return filtered


# ---------------------------------------------------------------------------
# Ranking result and convenience
# ---------------------------------------------------------------------------


@dataclass
class RankingResult:
    """Summary of a ranking operation.

    Attributes
    ----------
    total_candidates : int
        Total candidates before filtering.
    ranked_candidates : int
        Candidates after filtering.
    top_score : float
        Highest final score among ranked candidates.
    median_score : float
        Median final score among ranked candidates.
    candidates : list[CandidateScore]
        Ranked and filtered candidates.
    """

    total_candidates: int
    ranked_candidates: int
    top_score: float
    median_score: float
    candidates: list[CandidateScore]


def rank_and_filter(
    scored: list[CandidateScore],
    min_score: float = 0.0,
    top_n: int | None = None,
) -> RankingResult:
    """Rank and filter candidates in one step.

    Convenience function that chains filter_by_score → rank_candidates
    and wraps the result in a RankingResult.

    Parameters
    ----------
    scored : list[CandidateScore]
        Scored candidates.
    min_score : float
        Minimum score to keep (default 0.0).
    top_n : int or None
        Keep only top N after filtering.

    Returns
    -------
    RankingResult
        Ranking summary with filtered and ranked candidates.
    """
    total = len(scored)
    filtered = filter_by_score(scored, min_score=min_score, top_n=top_n)
    ranked = rank_candidates(filtered)

    if ranked:
        top = ranked[0].final_score
        scores = sorted(c.final_score for c in ranked)
        mid = len(scores) // 2
        if len(scores) % 2 == 0:
            median = (scores[mid - 1] + scores[mid]) / 2
        else:
            median = scores[mid]
    else:
        top = 0.0
        median = 0.0

    return RankingResult(
        total_candidates=total,
        ranked_candidates=len(ranked),
        top_score=top,
        median_score=median,
        candidates=ranked,
    )


def score_rank_and_filter(
    candidates: list,
    features_list: list[dict] | None = None,
    engine: ScoreEngine | None = None,
    min_score: float = 0.0,
    top_n: int | None = None,
) -> RankingResult:
    """End-to-end: score, rank, and filter candidates in one call.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to score.
    features_list : list[dict] or None
        Feature dicts for each candidate.
    engine : ScoreEngine or None
        Score engine (default: from config).
    min_score : float
        Minimum score to keep (default 0.0).
    top_n : int or None
        Keep only top N after filtering.

    Returns
    -------
    RankingResult
        Ranking summary with scored, ranked, and filtered candidates.
    """
    scored = score_candidates_batch(candidates, features_list, engine)
    return rank_and_filter(scored, min_score=min_score, top_n=top_n)
