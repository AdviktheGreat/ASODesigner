"""Tests for the scoring module: engine, components, convenience, ranking."""

import pytest

from aso_designer.scoring import (
    AccessibilityScoreComponent,
    GCContentScoreComponent,
    OffTargetScoreComponent,
    PositionalScoreComponent,
    ScoreComponent,
    ScoreContext,
    ScoreEngine,
    ScoreResult,
    ScoringError,
    SelfStructureScoreComponent,
    ThermodynamicsScoreComponent,
    VariantScoreComponent,
    create_score_engine,
)

# ---------------------------------------------------------------------------
# ScoreContext
# ---------------------------------------------------------------------------


class TestScoreContext:
    def test_defaults(self):
        ctx = ScoreContext()
        assert ctx.sequence == ""
        assert ctx.chemistry == "gapmer"
        assert ctx.gc_content == 0.0
        assert ctx.region == "unknown"
        assert ctx.accessibility_score == -1.0
        assert ctx.tm_score == -1.0
        assert ctx.off_target_score == -1.0
        assert ctx.variant_penalty == 0.0
        assert ctx.transcript_length == 0

    def test_custom_values(self):
        ctx = ScoreContext(
            sequence="ATGC",
            chemistry="mixmer",
            gc_content=0.5,
            region="CDS",
            accessibility_score=0.8,
            tm_score=0.9,
            binding_class="strong",
        )
        assert ctx.sequence == "ATGC"
        assert ctx.chemistry == "mixmer"
        assert ctx.accessibility_score == 0.8
        assert ctx.binding_class == "strong"


# ---------------------------------------------------------------------------
# ScoreResult
# ---------------------------------------------------------------------------


class TestScoreResult:
    def test_creation(self):
        r = ScoreResult("accessibility", 0.8, 0.25, 0.2)
        assert r.component_name == "accessibility"
        assert r.raw_score == 0.8
        assert r.weight == 0.25
        assert r.weighted_score == 0.2


# ---------------------------------------------------------------------------
# ScoreComponent protocol
# ---------------------------------------------------------------------------


class TestScoreComponentProtocol:
    def test_accessibility_is_component(self):
        comp = AccessibilityScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_thermo_is_component(self):
        comp = ThermodynamicsScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_gc_is_component(self):
        comp = GCContentScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_off_target_is_component(self):
        comp = OffTargetScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_variant_is_component(self):
        comp = VariantScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_self_structure_is_component(self):
        comp = SelfStructureScoreComponent()
        assert isinstance(comp, ScoreComponent)

    def test_positional_is_component(self):
        comp = PositionalScoreComponent()
        assert isinstance(comp, ScoreComponent)


# ---------------------------------------------------------------------------
# ScoreEngine
# ---------------------------------------------------------------------------


class TestScoreEngine:
    def test_init_default_weights(self):
        engine = ScoreEngine()
        assert "accessibility" in engine.weights
        assert "thermodynamics" in engine.weights
        assert "gc_content" in engine.weights

    def test_init_custom_weights(self):
        engine = ScoreEngine(weights={"accessibility": 0.5, "gc_content": 0.5})
        assert engine.weights["accessibility"] == 0.5
        assert engine.weights["gc_content"] == 0.5

    def test_register_component(self):
        engine = ScoreEngine(weights={"accessibility": 1.0})
        comp = AccessibilityScoreComponent()
        engine.register_component(comp)
        assert "accessibility" in engine._components

    def test_register_unknown_component_raises(self):
        engine = ScoreEngine(weights={"accessibility": 1.0})
        comp = GCContentScoreComponent()
        with pytest.raises(ScoringError, match="no weight"):
            engine.register_component(comp)

    def test_score_no_components_raises(self):
        engine = ScoreEngine(weights={"accessibility": 1.0})
        ctx = ScoreContext()
        with pytest.raises(ScoringError, match="No score components"):
            engine.score(ctx)

    def test_score_single_component(self):
        engine = ScoreEngine(weights={"accessibility": 1.0})
        engine.register_component(AccessibilityScoreComponent())
        ctx = ScoreContext(accessibility_score=0.8)
        result = engine.score(ctx)
        assert result.composite_score == pytest.approx(0.8)
        assert result.final_score == pytest.approx(0.8)
        assert len(result.component_results) == 1
        assert result.component_results[0].component_name == "accessibility"
        assert result.component_results[0].raw_score == pytest.approx(0.8)

    def test_score_clamps_high(self):
        engine = ScoreEngine(weights={"accessibility": 1.0}, chemistry_bonus=0.5)
        engine.register_component(AccessibilityScoreComponent())
        ctx = ScoreContext(accessibility_score=0.8)
        result = engine.score(ctx)
        # 0.8 + 0.5 = 1.3, clamped to 1.0
        assert result.final_score == pytest.approx(1.0)
        assert result.composite_score == pytest.approx(0.8)

    def test_score_clamps_low(self):
        engine = ScoreEngine(weights={"accessibility": 1.0}, chemistry_bonus=-0.5)
        engine.register_component(AccessibilityScoreComponent())
        ctx = ScoreContext(accessibility_score=0.3)
        result = engine.score(ctx)
        # 0.3 - 0.5 = -0.2, clamped to 0.0
        assert result.final_score == pytest.approx(0.0)

    def test_score_batch(self):
        engine = ScoreEngine(weights={"accessibility": 1.0})
        engine.register_component(AccessibilityScoreComponent())
        ctxs = [
            ScoreContext(accessibility_score=0.8),
            ScoreContext(accessibility_score=0.5),
        ]
        results = engine.score_batch(ctxs)
        assert len(results) == 2
        assert results[0].final_score == pytest.approx(0.8)
        assert results[1].final_score == pytest.approx(0.5)

    def test_score_raw_clamped(self):
        """Component returning >1.0 is clamped."""
        engine = ScoreEngine(weights={"gc_content": 1.0})
        engine.register_component(GCContentScoreComponent())
        ctx = ScoreContext(gc_content=0.5)
        result = engine.score(ctx)
        assert result.component_results[0].raw_score <= 1.0


# ---------------------------------------------------------------------------
# create_score_engine factory
# ---------------------------------------------------------------------------


class TestCreateScoreEngine:
    def test_default_engine(self):
        engine = create_score_engine()
        assert len(engine._components) == 7

    def test_custom_weights(self):
        engine = create_score_engine(weights={"accessibility": 0.5, "gc_content": 0.5})
        # Only components with matching weights are registered
        assert "accessibility" in engine._components
        assert "gc_content" in engine._components
        assert "thermodynamics" not in engine._components

    def test_chemistry_bonus(self):
        engine = create_score_engine(chemistry_bonus=0.1)
        assert engine.chemistry_bonus == pytest.approx(0.1)

    def test_register_all_false(self):
        engine = create_score_engine(register_all=False)
        assert len(engine._components) == 0


# ---------------------------------------------------------------------------
# AccessibilityScoreComponent
# ---------------------------------------------------------------------------


class TestAccessibilityScoreComponent:
    def test_name(self):
        assert AccessibilityScoreComponent().name == "accessibility"

    def test_computed_score(self):
        comp = AccessibilityScoreComponent()
        ctx = ScoreContext(accessibility_score=0.8)
        assert comp.score(ctx) == pytest.approx(0.8)

    def test_not_computed_neutral(self):
        comp = AccessibilityScoreComponent()
        ctx = ScoreContext(accessibility_score=-1.0)
        assert comp.score(ctx) == pytest.approx(0.5)

    def test_zero_score(self):
        comp = AccessibilityScoreComponent()
        ctx = ScoreContext(accessibility_score=0.0)
        assert comp.score(ctx) == pytest.approx(0.0)

    def test_perfect_score(self):
        comp = AccessibilityScoreComponent()
        ctx = ScoreContext(accessibility_score=1.0)
        assert comp.score(ctx) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# ThermodynamicsScoreComponent
# ---------------------------------------------------------------------------


class TestThermodynamicsScoreComponent:
    def test_name(self):
        assert ThermodynamicsScoreComponent().name == "thermodynamics"

    def test_strong_binding(self):
        comp = ThermodynamicsScoreComponent()
        ctx = ScoreContext(tm_score=1.0, binding_class="strong")
        # 1.0 * 0.6 + 1.0 * 0.4 = 1.0
        assert comp.score(ctx) == pytest.approx(1.0)

    def test_moderate_binding(self):
        comp = ThermodynamicsScoreComponent()
        ctx = ScoreContext(tm_score=1.0, binding_class="moderate")
        # 1.0 * 0.6 + 0.5 * 0.4 = 0.8
        assert comp.score(ctx) == pytest.approx(0.8)

    def test_weak_binding(self):
        comp = ThermodynamicsScoreComponent()
        ctx = ScoreContext(tm_score=1.0, binding_class="weak")
        # 1.0 * 0.6 + 0.0 * 0.4 = 0.6
        assert comp.score(ctx) == pytest.approx(0.6)

    def test_not_computed_neutral(self):
        comp = ThermodynamicsScoreComponent()
        ctx = ScoreContext(tm_score=-1.0)
        assert comp.score(ctx) == pytest.approx(0.5)

    def test_unknown_binding_class(self):
        comp = ThermodynamicsScoreComponent()
        ctx = ScoreContext(tm_score=0.5, binding_class="unknown")
        # 0.5 * 0.6 + 0.5 * 0.4 = 0.5 (unknown gets 0.5)
        assert comp.score(ctx) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# GCContentScoreComponent
# ---------------------------------------------------------------------------


class TestGCContentScoreComponent:
    def test_name(self):
        assert GCContentScoreComponent().name == "gc_content"

    def test_optimal_range(self):
        comp = GCContentScoreComponent()
        for gc in [0.30, 0.40, 0.50, 0.60]:
            ctx = ScoreContext(gc_content=gc)
            assert comp.score(ctx) == pytest.approx(1.0)

    def test_below_range(self):
        comp = GCContentScoreComponent(gc_min=0.30, gc_max=0.60)
        ctx = ScoreContext(gc_content=0.15)
        # 0.15 / 0.30 = 0.5
        assert comp.score(ctx) == pytest.approx(0.5)

    def test_above_range(self):
        comp = GCContentScoreComponent(gc_min=0.30, gc_max=0.60)
        ctx = ScoreContext(gc_content=0.80)
        # (1.0 - 0.80) / (1.0 - 0.60) = 0.20 / 0.40 = 0.5
        assert comp.score(ctx) == pytest.approx(0.5)

    def test_zero_gc(self):
        comp = GCContentScoreComponent(gc_min=0.30, gc_max=0.60)
        ctx = ScoreContext(gc_content=0.0)
        assert comp.score(ctx) == pytest.approx(0.0)

    def test_full_gc(self):
        comp = GCContentScoreComponent(gc_min=0.30, gc_max=0.60)
        ctx = ScoreContext(gc_content=1.0)
        assert comp.score(ctx) == pytest.approx(0.0)

    def test_custom_range(self):
        comp = GCContentScoreComponent(gc_min=0.40, gc_max=0.50)
        ctx = ScoreContext(gc_content=0.45)
        assert comp.score(ctx) == pytest.approx(1.0)
