"""
ASODesigner — In-silico antisense oligonucleotide (ASO) candidate design and ranking.

A computational pipeline for designing, scoring, and ranking ASO candidates
against a target transcript. Provides RNA accessibility analysis,
thermodynamic feature computation, transcriptome-wide off-target screening,
variant/SNP awareness, and interpretable multi-factor candidate scoring.

This is a research/candidate-ranking tool. It does not claim experimental
efficacy or safety for any candidate produced.
"""

from aso_designer.candidate import ASOCandidate
from aso_designer.chemistry import (
    ChemistryDesign,
    CleavageProfile,
    compute_cleavage_profile,
    create_chemistry_design,
)
from aso_designer.config_loader import load_config
from aso_designer.config_schema import ASOConfig
from aso_designer.motifs import MotifReport, assess_motifs
from aso_designer.pipeline import ASOPipeline, PipelineConfig
from aso_designer.thermodynamics import (
    DefectDuplexResult,
    calculate_defect_duplex,
    reset_chemistry_presets_cache,
    set_thermodynamics_cache,
)

__version__ = "1.0.0"
__author__ = "AdviktheGreat"
__email__ = "advikanand7@gmail.com"

__all__ = [
    "ASOConfig",
    "ASOCandidate",
    "ASOPipeline",
    "ChemistryDesign",
    "CleavageProfile",
    "DefectDuplexResult",
    "MotifReport",
    "PipelineConfig",
    "assess_motifs",
    "calculate_defect_duplex",
    "compute_cleavage_profile",
    "create_chemistry_design",
    "load_config",
    "reset_chemistry_presets_cache",
    "set_thermodynamics_cache",
    "__author__",
    "__email__",
    "__version__",
]
