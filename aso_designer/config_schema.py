"""Pydantic schema for ASODesigner configuration validation.

Inspired by GuideX's Cas13Config but redesigned for ASO-specific parameters
with nested models, cross-field validation, and detailed field descriptions.
"""

from pydantic import BaseModel, Field, confloat, conint, model_validator


class CandidateConfig(BaseModel):
    """Candidate generation parameters."""

    length: conint(ge=10, le=30) = Field(
        default=20, description="ASO oligonucleotide length (10–30 nt)."
    )
    step: conint(ge=1, le=10) = Field(default=1, description="Sliding window step size.")
    gc_min: confloat(ge=0.0, le=1.0) = Field(
        default=0.30, description="Minimum GC content fraction."
    )
    gc_max: confloat(ge=0.0, le=1.0) = Field(
        default=0.60, description="Maximum GC content fraction."
    )
    homopolymer_max: conint(ge=2, le=8) = Field(
        default=4, description="Max consecutive identical bases."
    )
    chemistry: str = Field(
        default="gapmer",
        description="Chemistry type: gapmer, mixmer, 2ome, moe, morpholino.",
    )

    @model_validator(mode="after")
    def validate_gc_range(self) -> "CandidateConfig":
        if self.gc_min > self.gc_max:
            raise ValueError("gc_min must be <= gc_max")
        return self

    @model_validator(mode="after")
    def validate_chemistry(self) -> "CandidateConfig":
        valid = {"gapmer", "mixmer", "2ome", "moe", "morpholino"}
        if self.chemistry not in valid:
            raise ValueError(f"Invalid chemistry '{self.chemistry}'. Must be one of {valid}.")
        return self


class GapmerConfig(BaseModel):
    """Gapmer-specific parameters."""

    gap_size: conint(ge=4, le=20) = Field(
        default=10, description="DNA center length (RNase H gap)."
    )
    wing_length: conint(ge=2, le=10) = Field(
        default=5, description="Modified wing length on each side."
    )
    wing_chemistry: str = Field(default="moe", description="Wing modification: 2ome or moe.")

    @model_validator(mode="after")
    def validate_wing_chemistry(self) -> "GapmerConfig":
        valid = {"2ome", "moe"}
        if self.wing_chemistry not in valid:
            raise ValueError(
                f"Invalid wing_chemistry '{self.wing_chemistry}'. Must be one of {valid}."
            )
        return self


class MixmerConfig(BaseModel):
    """Mixmer-specific parameters."""

    pattern: list[bool] = Field(
        default_factory=list,
        description="Modification pattern (True = modified).",
    )
    modification: str = Field(default="2ome", description="Modification type for marked positions.")


class AccessibilityConfig(BaseModel):
    """RNA accessibility analysis parameters."""

    window_size: conint(ge=50, le=500) = Field(
        default=150, description="RNAplfold sliding window size (nt)."
    )
    max_pairing_distance: conint(ge=20, le=500) = Field(
        default=100, description="Max pairing distance (nt)."
    )
    accessible_threshold: confloat(ge=0.0, le=1.0) = Field(
        default=0.5, description="Accessibility score threshold for 'accessible'."
    )


class ThermodynamicsConfig(BaseModel):
    """Thermodynamic feature parameters."""

    tm_min: float = Field(default=40.0, ge=20.0, le=100.0)
    tm_max: float = Field(default=65.0, ge=20.0, le=100.0)
    tm_optimal: float = Field(default=52.0, ge=20.0, le=100.0)
    na_conc: float = Field(default=0.1, gt=0.0, description="Sodium concentration (M).")
    oligo_conc: float = Field(default=1.0e-6, gt=0.0, description="Oligo concentration (M).")
    self_structure_max: float = Field(default=-3.0, description="Max self-structure ΔG (kcal/mol).")
    binding_min: float = Field(default=-5.0, description="Min net binding ΔG (kcal/mol).")

    @model_validator(mode="after")
    def validate_tm_range(self) -> "ThermodynamicsConfig":
        if self.tm_min > self.tm_max:
            raise ValueError("tm_min must be <= tm_max")
        if not (self.tm_min <= self.tm_optimal <= self.tm_max):
            raise ValueError("tm_optimal must be between tm_min and tm_max")
        return self


class OffTargetConfig(BaseModel):
    """Off-target analysis parameters."""

    max_mismatches: conint(ge=0, le=5) = Field(
        default=3, description="Max mismatches for off-target hits."
    )
    seed_start: conint(ge=0, le=10) = Field(default=2, description="Seed region start (0-based).")
    seed_end: conint(ge=1, le=15) = Field(default=8, description="Seed region end (exclusive).")
    seed_mismatch_weight: confloat(ge=0.0, le=10.0) = Field(
        default=2.0, description="Weight for seed-region mismatches."
    )
    max_hits: conint(ge=1, le=10000) = Field(
        default=50, description="Max off-target hits before filtering."
    )

    @model_validator(mode="after")
    def validate_seed_range(self) -> "OffTargetConfig":
        if self.seed_start >= self.seed_end:
            raise ValueError("seed_start must be < seed_end")
        return self


class VariantConfig(BaseModel):
    """Variant/SNP awareness parameters."""

    maf_threshold: confloat(ge=0.0, le=1.0) = Field(
        default=0.01, description="Min MAF to flag a variant."
    )
    seed_region_multiplier: confloat(ge=1.0, le=10.0) = Field(
        default=3.0, description="Penalty multiplier for seed-region variants."
    )
    flag_pathogenic: bool = Field(
        default=True, description="Flag candidates overlapping pathogenic variants."
    )


class ScoringWeightsConfig(BaseModel):
    """Scoring weight definitions. Weights should sum to 1.0."""

    accessibility: confloat(ge=0.0, le=1.0) = 0.25
    thermodynamics: confloat(ge=0.0, le=1.0) = 0.20
    gc_content: confloat(ge=0.0, le=1.0) = 0.10
    off_target: confloat(ge=0.0, le=1.0) = 0.20
    variant: confloat(ge=0.0, le=1.0) = 0.10
    self_structure: confloat(ge=0.0, le=1.0) = 0.10
    positional: confloat(ge=0.0, le=1.0) = 0.05

    @model_validator(mode="after")
    def validate_weights_sum(self) -> "ScoringWeightsConfig":
        total = (
            self.accessibility
            + self.thermodynamics
            + self.gc_content
            + self.off_target
            + self.variant
            + self.self_structure
            + self.positional
        )
        if abs(total - 1.0) > 0.001:
            raise ValueError(f"Scoring weights must sum to 1.0, got {total:.4f}.")
        return self


class ScoringConfig(BaseModel):
    """Scoring configuration."""

    weights: ScoringWeightsConfig = Field(default_factory=ScoringWeightsConfig)


class ReportingConfig(BaseModel):
    """Reporting parameters."""

    output_dir: str = Field(default="results", description="Output directory.")
    generate_html: bool = Field(default=True, description="Generate HTML report.")
    include_offtarget_details: bool = Field(
        default=True, description="Include full off-target details."
    )


class ASOConfig(BaseModel):
    """Top-level ASODesigner configuration."""

    candidate: CandidateConfig = Field(default_factory=CandidateConfig)
    gapmer: GapmerConfig = Field(default_factory=GapmerConfig)
    mixmer: MixmerConfig = Field(default_factory=MixmerConfig)
    accessibility: AccessibilityConfig = Field(default_factory=AccessibilityConfig)
    thermodynamics: ThermodynamicsConfig = Field(default_factory=ThermodynamicsConfig)
    off_target: OffTargetConfig = Field(default_factory=OffTargetConfig)
    variant: VariantConfig = Field(default_factory=VariantConfig)
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)
    reporting: ReportingConfig = Field(default_factory=ReportingConfig)
