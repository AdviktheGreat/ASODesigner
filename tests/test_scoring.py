"""Tests for the scoring module: engine, components, convenience, ranking."""

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.scoring import (
    AccessibilityScoreComponent,
    CandidateScore,
    CompositeScore,
    GCContentScoreComponent,
    OffTargetScoreComponent,
    PositionalScoreComponent,
    RankingResult,
    ScoreComponent,
    ScoreContext,
    ScoreEngine,
    ScoreResult,
    ScoringError,
    SelfStructureScoreComponent,
    ThermodynamicsScoreComponent,
    VariantScoreComponent,
    build_features_dict,
    create_engine_from_config,
    create_score_engine,
    filter_by_score,
    rank_and_filter,
    rank_candidates,
    rank_candidates_pareto,
    score_candidate,
    score_candidates_batch,
    score_rank_and_filter,
)
from aso_designer.thermodynamics import ThermodynamicSummary

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

    def test_register_all_registers_every_default_component(self):
        engine = create_score_engine(register_all=True)
        assert set(engine._components) == {
            "accessibility",
            "thermodynamics",
            "gc_content",
            "off_target",
            "variant",
            "self_structure",
            "positional",
        }


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

    def test_custom_tm_weight(self):
        comp = ThermodynamicsScoreComponent(tm_weight=1.0)
        ctx = ScoreContext(tm_score=0.7, binding_class="strong")
        # Tm-only: binding class contributes nothing
        assert comp.score(ctx) == pytest.approx(0.7)

    def test_custom_binding_class_scores(self):
        comp = ThermodynamicsScoreComponent(binding_class_scores={"strong": 0.9})
        ctx = ScoreContext(tm_score=1.0, binding_class="strong")
        # 1.0 * 0.6 + 0.9 * 0.4 = 0.96
        assert comp.score(ctx) == pytest.approx(0.96)

    def test_custom_unknown_class_score(self):
        comp = ThermodynamicsScoreComponent(unknown_class_score=0.2)
        ctx = ScoreContext(tm_score=0.5, binding_class="bogus")
        # 0.5 * 0.6 + 0.2 * 0.4 = 0.38
        assert comp.score(ctx) == pytest.approx(0.38)

    def test_invalid_tm_weight_raises(self):
        with pytest.raises(ValueError, match="tm_weight"):
            ThermodynamicsScoreComponent(tm_weight=1.5)


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


# ---------------------------------------------------------------------------
# OffTargetScoreComponent
# ---------------------------------------------------------------------------


class TestOffTargetScoreComponent:
    def test_name(self):
        assert OffTargetScoreComponent().name == "off_target"

    def test_computed_score(self):
        comp = OffTargetScoreComponent()
        ctx = ScoreContext(off_target_score=0.95)
        assert comp.score(ctx) == pytest.approx(0.95)

    def test_not_computed_neutral(self):
        comp = OffTargetScoreComponent()
        ctx = ScoreContext(off_target_score=-1.0)
        assert comp.score(ctx) == pytest.approx(0.5)

    def test_zero_score(self):
        comp = OffTargetScoreComponent()
        ctx = ScoreContext(off_target_score=0.0)
        assert comp.score(ctx) == pytest.approx(0.0)

    def test_perfect_score(self):
        comp = OffTargetScoreComponent()
        ctx = ScoreContext(off_target_score=1.0)
        assert comp.score(ctx) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# VariantScoreComponent
# ---------------------------------------------------------------------------


class TestVariantScoreComponent:
    def test_name(self):
        assert VariantScoreComponent().name == "variant"

    def test_no_variants(self):
        comp = VariantScoreComponent()
        ctx = ScoreContext(variant_penalty=0.0, variant_risk="none")
        assert comp.score(ctx) == pytest.approx(1.0)

    def test_no_data(self):
        comp = VariantScoreComponent()
        ctx = ScoreContext(variant_risk="")
        assert comp.score(ctx) == pytest.approx(1.0)

    def test_moderate_penalty(self):
        comp = VariantScoreComponent()
        ctx = ScoreContext(variant_penalty=0.3, variant_risk="moderate_risk")
        assert comp.score(ctx) == pytest.approx(0.7)

    def test_heavy_penalty(self):
        comp = VariantScoreComponent()
        ctx = ScoreContext(variant_penalty=0.8, variant_risk="high_risk")
        assert comp.score(ctx) == pytest.approx(0.2)

    def test_full_penalty_clamped(self):
        comp = VariantScoreComponent()
        ctx = ScoreContext(variant_penalty=1.5, variant_risk="high_risk")
        assert comp.score(ctx) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# SelfStructureScoreComponent
# ---------------------------------------------------------------------------


class TestSelfStructureScoreComponent:
    def test_name(self):
        assert SelfStructureScoreComponent().name == "self_structure"

    def test_no_data(self):
        comp = SelfStructureScoreComponent()
        ctx = ScoreContext(self_structure_penalty=0.0, dg_self=0.0)
        assert comp.score(ctx) == pytest.approx(1.0)

    def test_moderate_penalty(self):
        comp = SelfStructureScoreComponent()
        ctx = ScoreContext(self_structure_penalty=0.4, dg_self=-2.0)
        assert comp.score(ctx) == pytest.approx(0.6)

    def test_heavy_penalty(self):
        comp = SelfStructureScoreComponent()
        ctx = ScoreContext(self_structure_penalty=0.7, dg_self=-5.0)
        assert comp.score(ctx) == pytest.approx(0.3)

    def test_full_penalty_clamped(self):
        comp = SelfStructureScoreComponent()
        ctx = ScoreContext(self_structure_penalty=1.5, dg_self=-10.0)
        assert comp.score(ctx) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# PositionalScoreComponent
# ---------------------------------------------------------------------------


class TestPositionalScoreComponent:
    def test_name(self):
        assert PositionalScoreComponent().name == "positional"

    def test_cds_region(self):
        comp = PositionalScoreComponent()
        ctx = ScoreContext(region="CDS", start=0, end=20, length=20, transcript_length=100)
        assert comp.score(ctx) >= 1.0  # CDS gets 1.0, plus position bonus

    def test_5utr_region(self):
        comp = PositionalScoreComponent()
        ctx = ScoreContext(region="5UTR", start=0, end=20, length=20, transcript_length=1000)
        score = comp.score(ctx)
        assert 0.7 <= score <= 0.8  # 0.7 base + position bonus

    def test_3utr_region(self):
        comp = PositionalScoreComponent()
        ctx = ScoreContext(region="3UTR", start=0, end=20, length=20, transcript_length=1000)
        score = comp.score(ctx)
        assert 0.5 <= score <= 0.6  # 0.5 base + position bonus

    def test_unknown_region(self):
        comp = PositionalScoreComponent()
        ctx = ScoreContext(region="unknown", start=0, end=20, length=20, transcript_length=1000)
        score = comp.score(ctx)
        assert 0.3 <= score <= 0.4  # 0.3 base + position bonus

    def test_center_bonus(self):
        comp = PositionalScoreComponent(position_bonus=0.1)
        # Candidate at center of transcript
        ctx_center = ScoreContext(
            region="CDS", start=490, end=510, length=20, transcript_length=1000
        )
        # Candidate at edge
        ctx_edge = ScoreContext(region="CDS", start=0, end=20, length=20, transcript_length=1000)
        center_score = comp.score(ctx_center)
        edge_score = comp.score(ctx_edge)
        assert center_score >= edge_score

    def test_custom_region_weights(self):
        comp = PositionalScoreComponent(
            region_weights={"CDS": 0.9, "5UTR": 0.5},
            position_bonus=0.0,
        )
        ctx = ScoreContext(region="CDS", start=0, end=20, length=20, transcript_length=100)
        assert comp.score(ctx) == pytest.approx(0.9)

    def test_no_transcript_length(self):
        comp = PositionalScoreComponent(position_bonus=0.1)
        ctx = ScoreContext(region="CDS", start=0, end=20, length=20, transcript_length=0)
        assert comp.score(ctx) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# score_candidate
# ---------------------------------------------------------------------------


class TestScoreCandidate:
    def test_basic(self):
        cand = ASOCandidate(
            sequence="ATGCATGCATGCATGCATGC",
            start=100,
            end=120,
            length=20,
            chemistry="gapmer",
            gc_content=0.5,
            region="CDS",
        )
        features = {
            "accessibility_score": 0.8,
            "tm_score": 0.9,
            "binding_class": "strong",
            "off_target_score": 0.95,
            "variant_penalty": 0.2,
            "variant_risk": "moderate",
            "self_structure_penalty": 0.3,
            "dg_self": -2.0,
            "transcript_length": 2000,
        }
        result = score_candidate(cand, features)
        assert result.sequence == "ATGCATGCATGCATGCATGC"
        assert result.start == 100
        assert result.end == 120
        assert result.chemistry == "gapmer"
        assert 0.0 < result.final_score <= 1.0
        assert result.rank is None

    def test_empty_features(self):
        cand = ASOCandidate(
            sequence="ATGCATGCATGCATGCATGC",
            start=0,
            end=20,
            length=20,
            gc_content=0.5,
            region="CDS",
        )
        result = score_candidate(cand)
        assert 0.0 < result.final_score <= 1.0

    def test_custom_engine(self):
        cand = ASOCandidate(
            sequence="ATGC",
            start=0,
            end=4,
            length=4,
            gc_content=0.5,
        )
        engine = create_score_engine(weights={"gc_content": 1.0})
        result = score_candidate(cand, engine=engine)
        assert result.final_score == pytest.approx(1.0)

    def test_component_scores(self):
        cand = ASOCandidate(
            sequence="ATGCATGCATGCATGCATGC",
            start=0,
            end=20,
            length=20,
            gc_content=0.5,
            region="CDS",
        )
        features = {"accessibility_score": 0.8, "tm_score": 0.9, "binding_class": "strong"}
        result = score_candidate(cand, features)
        assert "accessibility" in result.component_scores
        assert "thermodynamics" in result.component_scores
        assert "gc_content" in result.component_scores

    def test_to_dict(self):
        cand = ASOCandidate(
            sequence="ATGC",
            start=0,
            end=4,
            length=4,
            gc_content=0.5,
        )
        result = score_candidate(cand)
        d = result.to_dict()
        assert d["sequence"] == "ATGC"
        assert d["start"] == 0
        assert d["end"] == 4
        assert "final_score" in d
        assert "component_scores" in d
        assert d["features"] == {}

    def test_features_carried_on_score(self):
        """score_candidate copies the feature dict onto the scored row."""
        cand = ASOCandidate(
            sequence="ATGCATGCATGCATGCATGC",
            start=0,
            end=20,
            length=20,
            gc_content=0.5,
            region="CDS",
        )
        features = {"cpg_count": 3, "motif_risk": "moderate"}
        result = score_candidate(cand, features)
        assert result.features == {"cpg_count": 3, "motif_risk": "moderate"}
        assert result.to_dict()["features"] == {"cpg_count": 3, "motif_risk": "moderate"}

    def test_features_dict_is_isolated(self):
        """Mutating the input features dict later must not affect the score."""
        cand = ASOCandidate(
            sequence="ATGC",
            start=0,
            end=4,
            length=4,
            gc_content=0.5,
        )
        features = {"cpg_count": 1}
        result = score_candidate(cand, features)
        features["cpg_count"] = 99
        assert result.features["cpg_count"] == 1


# ---------------------------------------------------------------------------
# score_candidates_batch
# ---------------------------------------------------------------------------


class TestScoreCandidatesBatch:
    def test_basic(self):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="3UTR"),
        ]
        features = [
            {"accessibility_score": 0.8, "tm_score": 0.9, "binding_class": "strong"},
            {"accessibility_score": 0.3, "tm_score": 0.5, "binding_class": "moderate"},
        ]
        results = score_candidates_batch(cands, features)
        assert len(results) == 2
        assert results[0].sequence != results[1].sequence

    def test_no_features(self):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0),
        ]
        results = score_candidates_batch(cands)
        assert len(results) == 2

    def test_length_mismatch_raises(self):
        cands = [ASOCandidate("ATGC", 0, 4, 4, gc_content=0.5)]
        with pytest.raises(ScoringError, match="length"):
            score_candidates_batch(cands, [{}, {}])

    def test_custom_engine(self):
        cands = [
            ASOCandidate("ATGC", 0, 4, 4, gc_content=0.5),
        ]
        engine = create_score_engine(weights={"gc_content": 1.0})
        results = score_candidates_batch(cands, engine=engine)
        assert results[0].final_score == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# build_features_dict
# ---------------------------------------------------------------------------


class TestBuildFeaturesDict:
    def test_from_thermodynamic_summary(self):
        thermo = ThermodynamicSummary(
            tm=55.0,
            tm_score=0.8,
            dg_duplex=-15.0,
            dg_self=-2.0,
            dg_target=-3.0,
            dg_net=-10.0,
            binding_class="strong",
            self_structure_penalty=0.5,
        )
        features = build_features_dict(
            thermo_summary=thermo,
            accessibility_score=0.7,
            transcript_length=2000,
        )
        assert features["tm"] == 55.0
        assert features["tm_score"] == 0.8
        assert features["dg_duplex"] == -15.0
        assert features["binding_class"] == "strong"
        assert features["self_structure_penalty"] == 0.5
        assert features["accessibility_score"] == 0.7
        assert features["transcript_length"] == 2000

    def test_empty(self):
        features = build_features_dict()
        assert features["accessibility_score"] == -1.0
        assert features["transcript_length"] == 0

    def test_with_off_target_tuple(self):
        features = build_features_dict(off_target_result=(None, 0.85, "specific"))
        assert features["off_target_score"] == pytest.approx(0.85)

    def test_only_accessibility(self):
        features = build_features_dict(accessibility_score=0.9)
        assert features["accessibility_score"] == pytest.approx(0.9)
        assert "tm" not in features


# ---------------------------------------------------------------------------
# create_engine_from_config
# ---------------------------------------------------------------------------


class TestCreateEngineFromConfig:
    def test_default_config(self):
        engine = create_engine_from_config()
        assert len(engine._components) == 7
        assert "accessibility" in engine._components

    def test_with_chemistry(self):
        engine = create_engine_from_config(chemistry="gapmer")
        assert len(engine._components) == 7

    def test_custom_config_path(self):
        engine = create_engine_from_config(
            config_path="aso_designer/config/aso_config.yaml",
            chemistry="mixmer",
        )
        assert len(engine._components) == 7


# ---------------------------------------------------------------------------
# rank_candidates
# ---------------------------------------------------------------------------


def _make_scored(score: float, seq: str = "ATGC") -> CandidateScore:
    """Helper: build a CandidateScore with a given final score."""
    return CandidateScore(
        sequence=seq,
        start=0,
        end=20,
        chemistry="gapmer",
        composite_score=CompositeScore(
            component_results=[],
            composite_score=score,
            chemistry_bonus=0.0,
            final_score=score,
        ),
        final_score=score,
    )


class TestRankCandidates:
    def test_basic_descending(self):
        scored = [
            _make_scored(0.5, "BBB"),
            _make_scored(0.9, "AAA"),
            _make_scored(0.3, "CCC"),
        ]
        ranked = rank_candidates(scored)
        assert ranked[0].sequence == "AAA"
        assert ranked[1].sequence == "BBB"
        assert ranked[2].sequence == "CCC"
        assert ranked[0].rank == 1
        assert ranked[1].rank == 2
        assert ranked[2].rank == 3

    def test_ascending(self):
        scored = [
            _make_scored(0.5, "BBB"),
            _make_scored(0.9, "AAA"),
            _make_scored(0.3, "CCC"),
        ]
        ranked = rank_candidates(scored, descending=False)
        assert ranked[0].sequence == "CCC"
        assert ranked[0].rank == 1

    def test_empty(self):
        ranked = rank_candidates([])
        assert ranked == []

    def test_single(self):
        scored = [_make_scored(0.8, "AAA")]
        ranked = rank_candidates(scored)
        assert len(ranked) == 1
        assert ranked[0].rank == 1

    def test_ties(self):
        scored = [
            _make_scored(0.5, "AAA"),
            _make_scored(0.5, "BBB"),
        ]
        ranked = rank_candidates(scored)
        assert ranked[0].rank == 1
        assert ranked[1].rank == 2


# ---------------------------------------------------------------------------
# filter_by_score
# ---------------------------------------------------------------------------


class TestFilterByScore:
    def test_min_score(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
            _make_scored(0.3, "C"),
        ]
        filtered = filter_by_score(scored, min_score=0.4)
        assert len(filtered) == 2
        assert all(c.final_score >= 0.4 for c in filtered)

    def test_max_score(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
            _make_scored(0.3, "C"),
        ]
        filtered = filter_by_score(scored, max_score=0.6)
        assert len(filtered) == 2
        assert all(c.final_score <= 0.6 for c in filtered)

    def test_top_n(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
            _make_scored(0.3, "C"),
        ]
        filtered = filter_by_score(scored, top_n=2)
        assert len(filtered) == 2
        assert filtered[0].final_score == pytest.approx(0.9)
        assert filtered[1].final_score == pytest.approx(0.5)

    def test_no_filter(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
        ]
        filtered = filter_by_score(scored)
        assert len(filtered) == 2

    def test_empty(self):
        filtered = filter_by_score([])
        assert filtered == []


# ---------------------------------------------------------------------------
# rank_and_filter
# ---------------------------------------------------------------------------


class TestRankAndFilter:
    def test_basic(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
            _make_scored(0.3, "C"),
            _make_scored(0.7, "D"),
        ]
        result = rank_and_filter(scored, min_score=0.4, top_n=3)
        assert isinstance(result, RankingResult)
        assert result.total_candidates == 4
        assert result.ranked_candidates == 3
        assert result.top_score == pytest.approx(0.9)
        assert result.candidates[0].rank == 1
        assert result.candidates[0].sequence == "A"

    def test_no_filter(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
        ]
        result = rank_and_filter(scored)
        assert result.total_candidates == 2
        assert result.ranked_candidates == 2

    def test_empty(self):
        result = rank_and_filter([])
        assert result.total_candidates == 0
        assert result.ranked_candidates == 0
        assert result.top_score == pytest.approx(0.0)
        assert result.median_score == pytest.approx(0.0)

    def test_median_odd(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.5, "B"),
            _make_scored(0.3, "C"),
        ]
        result = rank_and_filter(scored)
        assert result.median_score == pytest.approx(0.5)

    def test_median_even(self):
        scored = [
            _make_scored(0.9, "A"),
            _make_scored(0.7, "B"),
            _make_scored(0.5, "C"),
            _make_scored(0.3, "D"),
        ]
        result = rank_and_filter(scored)
        # Sorted: 0.3, 0.5, 0.7, 0.9 → median = (0.5 + 0.7) / 2 = 0.6
        assert result.median_score == pytest.approx(0.6)


# ---------------------------------------------------------------------------
# score_rank_and_filter
# ---------------------------------------------------------------------------


class TestScoreRankAndFilter:
    def test_basic(self):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="3UTR"),
            ASOCandidate("AATTAATTAATTAATTAATT", 40, 60, 20, gc_content=0.1, region="5UTR"),
        ]
        features = [
            {
                "accessibility_score": 0.8,
                "tm_score": 0.9,
                "binding_class": "strong",
                "off_target_score": 0.95,
            },
            {
                "accessibility_score": 0.3,
                "tm_score": 0.5,
                "binding_class": "moderate",
                "off_target_score": 0.7,
            },
            {
                "accessibility_score": 0.6,
                "tm_score": 0.7,
                "binding_class": "moderate",
                "off_target_score": 0.9,
            },
        ]
        result = score_rank_and_filter(cands, features, min_score=0.0, top_n=2)
        assert result.total_candidates == 3
        assert result.ranked_candidates == 2
        assert result.candidates[0].rank == 1
        assert result.candidates[1].rank == 2

    def test_all_filtered(self):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5),
        ]
        result = score_rank_and_filter(cands, min_score=1.5)
        assert result.total_candidates == 1
        assert result.ranked_candidates == 0

    def test_no_top_n(self):
        cands = [
            ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
            ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="CDS"),
        ]
        features = [
            {"accessibility_score": 0.8, "tm_score": 0.9, "binding_class": "strong"},
            {"accessibility_score": 0.3, "tm_score": 0.5, "binding_class": "moderate"},
        ]
        result = score_rank_and_filter(cands, features)
        assert result.ranked_candidates == 2
        assert result.candidates[0].rank == 1


# ---------------------------------------------------------------------------
# Deterministic tie-breaking
# ---------------------------------------------------------------------------


class TestRankCandidatesTieBreaking:
    def _candidate(self, seq, final, thermo, access, off=0.5):
        from aso_designer.scoring import CandidateScore, CompositeScore, ScoreResult

        results = [
            ScoreResult(
                component_name="thermodynamics",
                raw_score=thermo,
                weight=0.2,
                weighted_score=thermo * 0.2,
            ),
            ScoreResult(
                component_name="accessibility",
                raw_score=access,
                weight=0.25,
                weighted_score=access * 0.25,
            ),
            ScoreResult(
                component_name="off_target", raw_score=off, weight=0.2, weighted_score=off * 0.2
            ),
        ]
        composite = CompositeScore(
            component_results=results, composite_score=final, chemistry_bonus=0.0, final_score=final
        )
        return CandidateScore(
            sequence=seq,
            start=0,
            end=len(seq),
            chemistry="gapmer",
            composite_score=composite,
            final_score=final,
        )

    def test_ties_broken_by_secondary_keys(self):
        a = self._candidate("AAA", 0.8, thermo=0.9, access=0.5)
        b = self._candidate("CCC", 0.8, thermo=0.7, access=0.9)
        ranked = rank_candidates([b, a])
        assert ranked[0].sequence == "AAA"  # higher thermo wins the tie
        assert [c.rank for c in ranked] == [1, 2]

    def test_full_determinism_regardless_of_input_order(self):
        cands = [
            self._candidate("AAA", 0.8, 0.5, 0.5),
            self._candidate("CCC", 0.8, 0.5, 0.5),
            self._candidate("GGG", 0.8, 0.5, 0.5),
        ]
        out1 = rank_candidates(list(reversed(cands)))
        out2 = rank_candidates(list(cands))
        assert [c.sequence for c in out1] == [c.sequence for c in out2]

    def test_tie_groups_marked(self):
        a = self._candidate("AAA", 0.8, 0.9, 0.5)
        b = self._candidate("CCC", 0.8, 0.7, 0.9)
        c = self._candidate("GGG", 0.5, 0.9, 0.9)
        ranked = rank_candidates([c, b, a])
        assert ranked[0].tie_group == ranked[1].tie_group == 1
        assert ranked[2].tie_group == 2
        # to_dict carries it
        assert "tie_group" in ranked[0].to_dict()

    def test_distinct_scores_distinct_groups(self):
        a = self._candidate("AAA", 0.9, 0.5, 0.5)
        b = self._candidate("CCC", 0.8, 0.5, 0.5)
        ranked = rank_candidates([b, a])
        assert ranked[0].tie_group == 1
        assert ranked[1].tie_group == 2


# ---------------------------------------------------------------------------
# Pareto ranking
# ---------------------------------------------------------------------------


class TestRankCandidatesPareto:
    def _candidate(self, seq, final, thermo, access, off, variant=0.5, self_struct=0.5):
        from aso_designer.scoring import CandidateScore, CompositeScore, ScoreResult

        results = [
            ScoreResult(
                component_name="thermodynamics", raw_score=thermo, weight=0.2, weighted_score=thermo
            ),
            ScoreResult(
                component_name="accessibility", raw_score=access, weight=0.25, weighted_score=access
            ),
            ScoreResult(component_name="off_target", raw_score=off, weight=0.2, weighted_score=off),
            ScoreResult(
                component_name="variant", raw_score=variant, weight=0.1, weighted_score=variant
            ),
            ScoreResult(
                component_name="self_structure",
                raw_score=self_struct,
                weight=0.1,
                weighted_score=self_struct,
            ),
        ]
        composite = CompositeScore(
            component_results=results, composite_score=final, chemistry_bonus=0.0, final_score=final
        )
        return CandidateScore(
            sequence=seq,
            start=0,
            end=len(seq),
            chemistry="gapmer",
            composite_score=composite,
            final_score=final,
        )

    def test_dominated_candidate_ranked_last(self):
        dominated = self._candidate("AAA", 0.9, 0.2, 0.2, 0.2, 0.2, 0.2)
        strong = self._candidate("CCC", 0.7, 0.9, 0.9, 0.9, 0.9, 0.9)
        ranked = rank_candidates_pareto([dominated, strong])
        assert ranked[0].sequence == "CCC"
        assert ranked[0].rank == 1
        assert ranked[1].rank == 2

    def test_tradeoff_candidates_share_front(self):
        # A beats B on thermo; B beats A on off-target: neither dominates
        a = self._candidate("AAA", 0.8, 0.9, 0.5, 0.3)
        b = self._candidate("CCC", 0.7, 0.5, 0.5, 0.9)
        ranked = rank_candidates_pareto([a, b])
        assert {c.sequence for c in ranked[:2]} == {"AAA", "CCC"}
        # both are rank 1 candidates of front 1
        assert ranked[0].rank == 1
        assert ranked[1].rank == 2

    def test_equal_scores_share_tie_group(self):
        a = self._candidate("AAA", 0.8, 0.9, 0.5, 0.3)
        b = self._candidate("CCC", 0.8, 0.5, 0.5, 0.9)
        ranked = rank_candidates_pareto([a, b])
        assert ranked[0].tie_group == ranked[1].tie_group

    def test_unknown_mode_rejected(self):
        import pytest as _pytest

        with _pytest.raises(ValueError, match="Unknown rank_mode"):
            rank_and_filter([], rank_mode="bogus")
