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
