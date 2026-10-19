"""Pydantic schema for ASODesigner configuration validation.

Inspired by GuideX's Cas13Config but redesigned for ASO-specific parameters
with nested models, cross-field validation, and detailed field descriptions.
"""

from pydantic import BaseModel, Field, confloat, conint, field_validator, model_validator


class FetchConfig(BaseModel):
    """NCBI E-utilities access parameters."""

    email: str | None = Field(
        default=None,
        description=(
            "Email address for NCBI API identification (recommended by NCBI). "
            "Falls back to the NCBI_EMAIL environment variable."
        ),
    )
    api_key: str | None = Field(
        default=None,
        description=(
            "NCBI API key for higher rate limits (10 req/s vs 3 req/s). "
            "Falls back to the NCBI_API_KEY environment variable."
        ),
    )


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
    rnase_h_eligible: bool | None = Field(
        default=None,
        description=(
            "Override for RNase H eligibility (set by chemistry presets, "
            "e.g. morpholino: false). None derives it from the chemistry design."
        ),
    )
    self_dimer_min_stretch: conint(ge=0, le=20) = Field(
        default=0,
        description=(
            "Minimum 3'-anchored self-complementary stretch (nt) that REJECTS "
            "an ASO as a homodimer risk. 0 (default) records the stretch only "
            "as a report flag; e.g. 4 turns it into a hard filter."
        ),
    )
    chemistry: str = Field(
        default="gapmer",
        description="Chemistry type: gapmer, mixmer, 2ome, moe, morpholino.",
    )

    @model_validator(mode="after")
    def validate_gc_range(self) -> "CandidateConfig":
        """Validate that gc_min does not exceed gc_max."""
        if self.gc_min > self.gc_max:
            raise ValueError("gc_min must be <= gc_max")
        return self

    @model_validator(mode="after")
    def validate_chemistry(self) -> "CandidateConfig":
        """Validate that chemistry is a supported type."""
        valid = {"gapmer", "mixmer", "2ome", "moe", "morpholino", "lna"}
        if self.chemistry not in valid:
            raise ValueError(f"Invalid chemistry '{self.chemistry}'. Must be one of {valid}.")
        return self


class LnaConfig(BaseModel):
    """LNA-specific parameters."""

    pattern: list[bool] = Field(
        default_factory=list,
        description="Per-position LNA flags (True = LNA, False = DNA). Empty = default every-third pattern.",
    )


class GapmerConfig(BaseModel):
    """Gapmer-specific parameters."""

    gap_size: conint(ge=4, le=20) = Field(
        default=10, description="DNA center length (RNase H gap)."
    )
    wing_length: conint(ge=2, le=10) = Field(
        default=5, description="Modified wing length on each side."
    )
    backbone: str = Field(
        default="ps",
        description="Internucleotide backbone: 'ps' (full phosphorothioate, default) or 'po'.",
    )
    wing_chemistry: str = Field(default="moe", description="Wing modification: 2ome or moe.")

    @model_validator(mode="after")
    def validate_wing_chemistry(self) -> "GapmerConfig":
        """Validate that wing chemistry is 2ome or moe."""
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
    modification: str = Field(
        default="2ome",
        description="Modification type for marked positions: '2ome' or 'moe'.",
    )

    @field_validator("modification")
    @classmethod
    def validate_modification(cls, v: str) -> str:
        """Restrict mixmer modifications to the supported enum."""
        allowed = {"2ome", "moe"}
        if v not in allowed:
            raise ValueError(f"mixmer.modification must be one of {sorted(allowed)}, got {v!r}")
        return v


class AccessibilityConfig(BaseModel):
    """RNA accessibility analysis parameters."""

    method: str = Field(
        default="rnaplfold",
        description=(
            "Accessibility tool: 'rnaplfold' (sliding-window unpaired "
            "probability) or 'rnaup' (per-site opening energy)."
        ),
    )
    unpaired_region: conint(ge=1, le=15) = Field(
        default=4,
        description="Unpaired-region length for the RNAup method (-u).",
    )
    unpaired_lengths: list[conint(ge=1, le=15)] = Field(
        default_factory=lambda: [1],
        description=(
            "Unpaired-region lengths for the RNAplfold method (-u 1,2,...). "
            "The default [1] keeps per-nucleotide accessibility; e.g. "
            "[1, 2, 3, 4] adds multi-nucleotide stretch probabilities."
        ),
    )
    window_size: conint(ge=50, le=500) = Field(
        default=150, description="RNAplfold sliding window size (nt)."
    )
    max_pairing_distance: conint(ge=20, le=500) = Field(
        default=100, description="Max pairing distance (nt)."
    )
    accessible_threshold: confloat(ge=0.0, le=1.0) = Field(
        default=0.5, description="Accessibility score threshold for 'accessible'."
    )

    @field_validator("method")
    @classmethod
    def validate_method(cls, v: str) -> str:
        """Restrict the accessibility method to supported tools."""
        allowed = {"rnaplfold", "rnaup"}
        if v not in allowed:
            raise ValueError(f"accessibility.method must be one of {sorted(allowed)}, got {v!r}")
        return v

    @field_validator("unpaired_lengths")
    @classmethod
    def validate_unpaired_lengths(cls, v: list[int]) -> list[int]:
        """Require at least one unique, positive unpaired-region length."""
        if not v:
            raise ValueError("accessibility.unpaired_lengths must not be empty")
        if len(set(v)) != len(v):
            raise ValueError("accessibility.unpaired_lengths must not contain duplicates")
        return sorted(v)


class ThermodynamicsConfig(BaseModel):
    """Thermodynamic feature parameters."""

    method: str = Field(
        default="nearest_neighbor",
        description=(
            "Computation method: 'nearest_neighbor' (fast, no ViennaRNA), "
            "'rnafold' (ViennaRNA for self/target structure), or "
            "'rnaduplex' (ViennaRNA RNAduplex for duplex ΔG)."
        ),
    )
    tm_min: float = Field(default=40.0, ge=20.0, le=100.0)
    tm_max: float = Field(default=65.0, ge=20.0, le=100.0)
    tm_optimal: float = Field(default=52.0, ge=20.0, le=100.0)
    na_conc: float = Field(default=0.1, gt=0.0, description="Sodium concentration (M).")
    oligo_conc: float = Field(default=1.0e-6, gt=0.0, description="Oligo concentration (M).")
    mg_conc: float = Field(
        default=0.0,
        ge=0.0,
        description=(
            "Magnesium concentration (M). When > 0, Tm uses the Owczarzy 2008 "
            "divalent-ion correction (which replaces the monovalent-only "
            "SantaLucia correction)."
        ),
    )
    dntp_conc: float = Field(
        default=0.0,
        ge=0.0,
        description=(
            "dNTP concentration (M). dNTPs chelate Mg2+; free Mg2+ is used "
            "in the Owczarzy 2008 correction."
        ),
    )
    self_structure_max: float = Field(default=-3.0, description="Max self-structure ΔG (kcal/mol).")
    self_structure_nn_fraction: confloat(gt=0.0, le=1.0) = Field(
        default=0.3,
        description=(
            "Fraction of the perfect-duplex ΔG attributed to intramolecular "
            "self-structure in the nearest-neighbor approximation "
            "(used when method is 'nearest_neighbor')."
        ),
    )
    binding_min: float = Field(default=-5.0, description="Min net binding ΔG (kcal/mol).")
    cache_dir: str | None = Field(
        default=None,
        description=(
            "Optional directory for caching ViennaRNA duplex/self-structure "
            "results on disk (reduces repeat invocations across runs)."
        ),
    )

    @model_validator(mode="after")
    def validate_method(self) -> "ThermodynamicsConfig":
        """Validate that the computation method is supported."""
        valid = {"nearest_neighbor", "rnafold", "rnaduplex"}
        if self.method not in valid:
            raise ValueError(f"Invalid method '{self.method}'. Must be one of {valid}.")
        return self

    @model_validator(mode="after")
    def validate_tm_range(self) -> "ThermodynamicsConfig":
        """Validate that Tm bounds are consistent and optimal is within range."""
        if self.tm_min > self.tm_max:
            raise ValueError("tm_min must be <= tm_max")
        if not (self.tm_min <= self.tm_optimal <= self.tm_max):
            raise ValueError("tm_optimal must be between tm_min and tm_max")
        return self


class OffTargetConfig(BaseModel):
    """Off-target analysis parameters."""

    transcriptome_path: str | None = Field(
        default=None,
        description=(
            "Path to a transcriptome FASTA file for off-target search. "
            "If None, off-target analysis is skipped and neutral scores are used."
        ),
    )
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
        """Validate that seed_start is less than seed_end."""
        if self.seed_start >= self.seed_end:
            raise ValueError("seed_start must be < seed_end")
        return self


class VariantConfig(BaseModel):
    """Variant/SNP awareness parameters."""

    vcf_path: str | None = Field(
        default=None,
        description=(
            "Path to a VCF file with known variants for the target. "
            "If None, variants are fetched from dbSNP/ClinVar APIs (when a "
            "gene name is available); if that fails, variant analysis is skipped."
        ),
    )
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
        """Validate that all scoring weights sum to 1.0."""
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


# Named weight profiles: quick presets over the seven component weights
# so users don't hand-edit floats. "balanced" equals the schema defaults.
WEIGHT_PROFILES: dict[str, dict[str, float]] = {
    "balanced": {
        "accessibility": 0.25,
        "thermodynamics": 0.20,
        "gc_content": 0.10,
        "off_target": 0.20,
        "variant": 0.10,
        "self_structure": 0.10,
        "positional": 0.05,
    },
    "efficacy": {
        "accessibility": 0.35,
        "thermodynamics": 0.30,
        "gc_content": 0.05,
        "off_target": 0.10,
        "variant": 0.05,
        "self_structure": 0.10,
        "positional": 0.05,
    },
    "specificity": {
        "accessibility": 0.15,
        "thermodynamics": 0.10,
        "gc_content": 0.10,
        "off_target": 0.45,
        "variant": 0.15,
        "self_structure": 0.05,
        "positional": 0.00,
    },
    "safety": {
        "accessibility": 0.15,
        "thermodynamics": 0.15,
        "gc_content": 0.10,
        "off_target": 0.15,
        "variant": 0.35,
        "self_structure": 0.10,
        "positional": 0.00,
    },
}


class ScoringConfig(BaseModel):
    """Scoring configuration."""

    weights: ScoringWeightsConfig = Field(default_factory=ScoringWeightsConfig)
    profile: str = Field(
        default="balanced",
        description=(
            "Named weight profile: balanced, efficacy, specificity, or "
            "safety. Applied by the config loader unless the user set "
            "explicit weights; 'balanced' keeps the schema defaults."
        ),
    )
    ignore_neutral_components: bool = Field(
        default=False,
        description=(
            "Recompute final scores using only components backed by real "
            "data (weights of neutral components are dropped and the rest "
            "renormalized) before ranking."
        ),
    )
    chemistry_bonus: confloat(ge=0.0, le=0.2) = Field(
        default=0.0,
        description="Flat bonus added to the composite score (set by chemistry presets).",
    )
    rank_mode: str = Field(
        default="weighted",
        description="Ranking mode: 'weighted' (weighted sum) or 'pareto' (non-dominated sort).",
    )

    @field_validator("profile")
    @classmethod
    def validate_profile(cls, v: str) -> str:
        """Restrict the weight-profile name."""
        allowed = set(WEIGHT_PROFILES)
        if v not in allowed:
            raise ValueError(f"scoring.profile must be one of {sorted(allowed)}, got {v!r}")
        return v

    @field_validator("rank_mode")
    @classmethod
    def validate_rank_mode(cls, v: str) -> str:
        """Restrict the ranking mode."""
        allowed = {"weighted", "pareto"}
        if v not in allowed:
            raise ValueError(f"scoring.rank_mode must be one of {sorted(allowed)}, got {v!r}")
        return v


class ReportingConfig(BaseModel):
    """Reporting parameters."""

    output_dir: str = Field(default="results", description="Output directory.")
    generate_html: bool = Field(default=True, description="Generate HTML report.")
    include_offtarget_details: bool = Field(
        default=True, description="Include full off-target details."
    )


class ASOConfig(BaseModel):
    """Top-level ASODesigner configuration."""

    fetch: FetchConfig = Field(default_factory=FetchConfig)
    candidate: CandidateConfig = Field(default_factory=CandidateConfig)
    gapmer: GapmerConfig = Field(default_factory=GapmerConfig)
    lna: LnaConfig = Field(default_factory=LnaConfig)
    mixmer: MixmerConfig = Field(default_factory=MixmerConfig)
    accessibility: AccessibilityConfig = Field(default_factory=AccessibilityConfig)
    thermodynamics: ThermodynamicsConfig = Field(default_factory=ThermodynamicsConfig)
    off_target: OffTargetConfig = Field(default_factory=OffTargetConfig)
    variant: VariantConfig = Field(default_factory=VariantConfig)
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)
    reporting: ReportingConfig = Field(default_factory=ReportingConfig)

    @model_validator(mode="after")
    def validate_gapmer_geometry(self) -> "ASOConfig":
        """Reject gapmer configs whose length does not match wings + gap."""
        if self.candidate.chemistry == "gapmer":
            expected = 2 * self.gapmer.wing_length + self.gapmer.gap_size
            if self.candidate.length != expected:
                raise ValueError(
                    f"Invalid gapmer geometry: candidate.length {self.candidate.length} != "
                    f"2 * gapmer.wing_length {self.gapmer.wing_length} + gapmer.gap_size "
                    f"{self.gapmer.gap_size} (expected {expected}). Adjust candidate.length, "
                    "gapmer.gap_size, or gapmer.wing_length."
                )
        return self
