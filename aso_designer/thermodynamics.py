"""Thermodynamic calculations for ASO candidate evaluation.

Implements nearest-neighbor (NN) model for melting temperature (Tm)
prediction, duplex ΔG, self-structure ΔG, target structure ΔG, and
net binding energy. Uses SantaLucia 1998 unified NN parameters for
DNA/DNA, Turner 2004 parameters for RNA/RNA, and the Banerjee et al.
2020 (corrigendum 2021) RNA/DNA hybrid parameters for DNA-containing
ASOs binding RNA targets. ViennaRNA tools (RNAfold, RNAduplex,
RNAcofold) are used for structure-based ΔG calculations where the NN
model is insufficient.

Scientific note: Tm and ΔG values are computational predictions for
candidate ranking, not experimental measurements. The NN model
assumes 1 M NaCl monovalent salt and non-self-complementary sequences
unless otherwise specified.
"""

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

from aso_designer.config_loader import CHEMISTRY_PRESETS_PATH
from aso_designer.utils.cache_key import cache_key
from aso_designer.utils.exceptions import ThermodynamicsError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Optional on-disk cache for ViennaRNA results (see set_thermodynamics_cache).
_VIENNA_CACHE_DIR: Path | None = None
_VIENNA_CACHE_ENABLED: bool = False


def set_thermodynamics_cache(cache_dir: str | Path | None) -> None:
    """Enable or disable the on-disk cache for ViennaRNA results.

    Parameters
    ----------
    cache_dir : str | Path | None
        Directory for cached duplex/self/target structure results, or
        None to disable caching. Results are stored as small JSON files
        keyed by a hash of the tool and input sequence, mirroring the
        accessibility-stage cache.
    """
    global _VIENNA_CACHE_DIR, _VIENNA_CACHE_ENABLED
    if cache_dir is None:
        _VIENNA_CACHE_DIR = None
        _VIENNA_CACHE_ENABLED = False
        return
    _VIENNA_CACHE_DIR = Path(cache_dir)
    try:
        _VIENNA_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logger.warning("Could not create thermodynamics cache dir: %s", e)
        _VIENNA_CACHE_ENABLED = False
        return
    _VIENNA_CACHE_ENABLED = True


def _vienna_cache_path(tool: str, key: str) -> Path | None:
    """Return the cache file path for a ViennaRNA call, or None if disabled."""
    if not _VIENNA_CACHE_ENABLED or _VIENNA_CACHE_DIR is None:
        return None
    digest = cache_key(tool, key)
    return _VIENNA_CACHE_DIR / f"{tool}_{digest}.json"


def _load_vienna_cache(path: Path | None) -> dict | None:
    """Load a cached ViennaRNA result; None on miss or corruption."""
    if path is None:
        return None
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _save_vienna_cache(path: Path | None, data: dict) -> None:
    """Persist a ViennaRNA result; failures are non-fatal."""
    if path is None:
        return
    try:
        with open(path, "w") as f:
            json.dump(data, f)
    except OSError as e:
        logger.warning("Failed to write thermodynamics cache %s: %s", path, e)


# Per-process cache for the bundled chemistry presets, loaded from
# CHEMISTRY_PRESETS_PATH (owned by config_loader).
_CHEMISTRY_PRESETS_CACHE: dict | None = None


def _load_chemistry_presets() -> dict:
    """Load the chemistry presets YAML (cached after the first load).

    Returns
    -------
    dict
        Chemistry presets keyed by chemistry name.

    Raises
    ------
    OSError, yaml.YAMLError
        If the bundled presets file is missing, unreadable, or invalid.
    """
    global _CHEMISTRY_PRESETS_CACHE
    if _CHEMISTRY_PRESETS_CACHE is None:
        with open(CHEMISTRY_PRESETS_PATH) as f:
            _CHEMISTRY_PRESETS_CACHE = yaml.safe_load(f) or {}
    return _CHEMISTRY_PRESETS_CACHE


def reset_chemistry_presets_cache() -> None:
    """Drop the cached chemistry presets (mainly for tests)."""
    global _CHEMISTRY_PRESETS_CACHE
    _CHEMISTRY_PRESETS_CACHE = None


# Default reaction temperature for Tm/ΔG calculations (°C)
DEFAULT_TEMPERATURE_C = 37.0

# Default fraction of the perfect-duplex ΔG attributed to intramolecular
# self-structure in the nearest-neighbor approximation.
_DEFAULT_NN_FRACTION = 0.3

# Gas constant in cal/(mol·K)
R = 1.987

# Reference temperature for ΔG calculations (37 °C in Kelvin)
T_REF = 273.15 + DEFAULT_TEMPERATURE_C

# ---------------------------------------------------------------------------
# SantaLucia 1998 unified NN parameters for DNA/DNA
# Reference: SantaLucia (1998) PNAS 95:1460-1465
# Units: ΔH in kcal/mol, ΔS in cal/(mol·K)
# Key format: "XY" = 5'-X Y-3' / 3'-Y' X'-5' (top strand 5'→3')
# ---------------------------------------------------------------------------

_DNA_NN: dict[str, tuple[float, float]] = {
    # NN pair (5'-3'/3'-5'): (ΔH, ΔS)
    # 10 unique parameters from SantaLucia 1998, expanded to all 16
    # dinucleotides using reverse-complement symmetry.
    "AA": (-7.9, -22.2),
    "AC": (-8.4, -22.4),  # = GT
    "AG": (-7.8, -21.0),  # = CT
    "AT": (-7.2, -20.4),
    "CA": (-8.5, -22.7),
    "CC": (-8.0, -19.9),  # = GG
    "CG": (-10.6, -27.2),
    "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "GT": (-8.4, -22.4),
    "TA": (-7.2, -21.3),
    "TC": (-8.2, -22.2),  # = GA
    "TG": (-8.5, -22.7),  # = CA
    "TT": (-7.9, -22.2),  # = AA
}

# Initiation parameters (SantaLucia 1998)
# Terminal G·C: ΔH = 0.1, ΔS = -2.8
# Terminal A·T: ΔH = 2.3, ΔS = 4.1
# Symmetry correction (self-complementary): ΔH = 0, ΔS = -1.4
_DNA_INIT_GC = (0.1, -2.8)
_DNA_INIT_AT = (2.3, 4.1)
_DNA_SYMMETRY = (0.0, -1.4)

# ---------------------------------------------------------------------------
# Turner 2004 NN parameters for RNA/RNA
# Reference: Xia et al. (1998) Biochemistry 37:14719-14735, as compiled in
# the Turner 2004 Nearest Neighbor Database (NNDB),
# https://rna.urmc.rochester.edu/NNDB/ (Turner & Mathews 2010).
# Units: ΔH in kcal/mol, ΔS in cal/(mol·K). ΔS values are derived from the
# published ΔG°37 and ΔH° via ΔS = (ΔH − ΔG°37) × 1000 / 310.15.
# Key format: "XY" = 5'-X Y-3' / 3'-revcomp(XY)-5' (top strand 5'→3');
# entries marked "= ZW" are equal to ZW by reverse-complement symmetry.
# ---------------------------------------------------------------------------

_RNA_NN: dict[str, tuple[float, float]] = {
    # NN pair (5'-3'/3'-5'): (ΔH, ΔS)
    "AA": (-6.82, -19.0),  # = UU
    "AC": (-11.40, -29.5),  # = GU
    "AG": (-10.48, -27.1),  # = CU
    "AU": (-9.38, -26.7),
    "CA": (-10.44, -26.9),  # = UG
    "CC": (-13.39, -32.7),  # = GG
    "CG": (-10.64, -26.7),
    "CU": (-10.48, -27.1),
    "GA": (-12.44, -32.5),  # = UC
    "GC": (-14.88, -36.9),
    "GG": (-13.39, -32.7),
    "GU": (-11.40, -29.5),
    "UA": (-7.69, -20.5),
    "UC": (-12.44, -32.5),
    "UG": (-10.44, -26.9),
    "UU": (-6.82, -19.0),
}

# RNA initiation parameters (Turner 2004 / Xia et al. 1998, from the NNDB)
# Intermolecular initiation: ΔH = +3.61 kcal/mol, ΔG°37 = +4.09 kcal/mol
#   → ΔS = (3.61 − 4.09) × 1000 / 310.15 = −1.5 cal/(mol·K)
# Applied once per bimolecular duplex.
_RNA_INIT_INTERMOLECULAR = (3.61, -1.5)
# AU end penalty: ΔH = +3.72 kcal/mol, ΔG°37 = +0.45 kcal/mol
#   → ΔS = (3.72 − 0.45) × 1000 / 310.15 = +10.5 cal/(mol·K)
# Applied once per terminal A·U base pair.
_RNA_AU_END = (3.72, 10.5)
# Symmetry correction for self-complementary RNA duplexes:
# ΔG°37 = +0.43 kcal/mol, ΔH = 0 → ΔS = −1.4 cal/(mol·K)
_RNA_SYMMETRY = (0.0, -1.4)

# ---------------------------------------------------------------------------
# RNA/DNA hybrid NN parameters (RNA target strand paired with DNA ASO strand)
# Reference: Banerjee et al. (2020) Nucleic Acids Res 48:12042-12055,
# "Improved nearest-neighbor parameters for the stability of RNA/DNA hybrids
# under a physiological condition", with the corrigendum (NAR 2021, 49:10796,
# gkab780) that recalculated ΔS values. Values are the re-derived 1 M NaCl
# parameter set (Supplementary Table S3, "parameters calculated by new
# program" column, corrigendum revision); ΔS values follow the same
# Gibbs-relation convention as the RNA table above.
# Units: ΔH in kcal/mol, ΔS in cal/(mol·K).
# Key format: "XY" = RNA strand 5'-X Y-3' paired with the antiparallel DNA
# strand (e.g. "AG" = rAG/dCT). Unlike the DNA and RNA tables this set is
# NOT reverse-complement symmetric, so all 16 keys are independent entries.
# Model notes (per the publication):
# - A single helix-initiation term applies to every duplex (no separate
#   terminal GC/AT initiation and no terminal AU penalty; the paper found an
#   explicit terminal rA-dT term did not improve predictions).
# - No symmetry correction: an RNA/DNA hybrid is never self-complementary
#   (the two strands are chemically different), consistent with the paper's
#   Tm fits using ln(Ct/4).
# ---------------------------------------------------------------------------

_HYBRID_NN: dict[str, tuple[float, float]] = {
    # NN pair on the RNA strand (5'-3'), paired with DNA (3'-5'): (ΔH, ΔS)
    "AA": (-7.5, -21.0),  # rAA/dTT
    "AC": (-8.5, -21.6),  # rAC/dGT
    "AG": (-8.6, -22.6),  # rAG/dCT
    "AU": (-7.5, -21.9),  # rAU/dAT
    "CA": (-8.4, -22.9),  # rCA/dTG
    "CC": (-7.9, -19.0),  # rCC/dGG
    "CG": (-14.8, -41.6),  # rCG/dCG
    "CU": (-4.4, -10.6),  # rCU/dAG
    "GA": (-9.4, -24.5),  # rGA/dTC
    "GC": (-10.5, -25.5),  # rGC/dGC
    "GG": (-12.5, -31.6),  # rGG/dCC
    "GU": (-9.5, -25.8),  # rGU/dAC
    "UA": (-9.5, -28.4),  # rUA/dTA
    "UC": (-8.1, -22.2),  # rUC/dGA
    "UG": (-8.6, -23.2),  # rUG/dCA
    "UU": (-5.6, -16.8),  # rUU/dAA
}

# Helix initiation for RNA/DNA hybrids (1 M NaCl):
# ΔH = +1.9 kcal/mol, ΔG°37 = +3.1 kcal/mol
#   → ΔS = (1.9 − 3.1) × 1000 / 310.15 = −3.9 cal/(mol·K)
_HYBRID_INIT = (1.9, -3.9)

# Complement map used to derive the RNA target strand from a DNA ASO
# sequence: the RNA strand antiparallel to a DNA 5'→3' strand, read in its
# own 5'→3' direction, is the DNA strand's reverse complement (with T→U).
_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# ---------------------------------------------------------------------------
# DNA single-mismatch / dangling-end defect parameters (1 M NaCl)
# ---------------------------------------------------------------------------
# Source: Nearest Neighbor Database (Mittal, Turner & Mathews 2024), DNA
# parameter tables downloaded from
#   https://rna.urmc.rochester.edu/NNDB/parameter_tables/dna/
# Underlying optical-melting data: 1x1 internal loops (single internal
# mismatches; Allawi & SantaLucia 1998 and related compilations), terminal
# mismatches (SantaLucia-compilation tables), and dangling ends
# (Bommarito, Peyret & SantaLucia 2000, Biochemistry 39:3423).
#
# Units: dG37 and dH in kcal/mol at 37 C, 1 M NaCl; dS is derived at import
# via dS = (dH - dG37) * 1000 / 310.15. A value of None (the '.' entries in
# the source tables) marks a disallowed motif. Blocks of the 1x1 internal-
# loop tables whose every entry equals 3.0 are the published fallback
# penalties for wobble-closed motifs with insufficient experimental data
# (dS = 0 for those entries).
#
# Key conventions (bottom strand written 3'->5' under the top strand):
# - _INT11_DG37/DH: key "LT/LB|RT/RB" is expanded as (left closing pair
#   "LT/LB", right closing pair "RT/RB"); grid key "XY" = (top mismatch
#   base X, bottom mismatch base Y). The tabulated value REPLACES the two
#   flanking Watson-Crick stacks (NNDB internal-loop composition).
# - _TERMINAL_MISMATCH_DG37/DH: key "TB" = closing pair (top penultimate,
#   bottom penultimate); grid key "XY" = (top terminal base X, bottom
#   terminal base Y). The value REPLACES the terminal stack (NNDB
#   "WC helix + terminal-mismatch" composition).
# - _DANGLE_DG37/DH: key "TB|side" = closing pair and side ("top_3p" =
#   3'-terminal dangle on the top strand; "bottom_5p" = 5'-terminal dangle
#   on the bottom strand); value keyed by the dangling base. The value
#   ADDS to the paired-region helix (NNDB dangle composition). Only WC and
#   G.T-wobble closing pairs are tabulated.
# ---------------------------------------------------------------------------

_INT11_DG37: dict[str, dict[str, float | None]] = {
    "AT/AT": {
        "AA": 1.3,
        "AC": 2.2,
        "AG": 0.9,
        "AT": 1.0,
        "CA": 1.7,
        "CC": 2.4,
        "CG": 1.0,
        "CT": 1.4,
        "GA": 0.4,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.7,
        "TG": 2.0,
        "TT": 1.4,
    },
    "AT/CG": {
        "AA": 0.8,
        "AC": 1.4,
        "AG": -0.4,
        "AT": 1.0,
        "CA": 1.6,
        "CC": 2.1,
        "CG": 1.0,
        "CT": 1.6,
        "GA": -0.2,
        "GC": 1.0,
        "GG": -1.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 1.1,
    },
    "AT/GC": {
        "AA": 1.0,
        "AC": 1.7,
        "AG": 0.3,
        "AT": 1.0,
        "CA": 1.5,
        "CC": 2.0,
        "CG": 1.0,
        "CT": 1.0,
        "GA": 0.1,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 0.6,
    },
    "AT/GT": {
        "AA": 2.0,
        "AC": 2.7,
        "AG": 1.3,
        "AT": 1.0,
        "CA": 2.5,
        "CC": 3.0,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 1.1,
        "GC": 1.0,
        "GG": 0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.4,
        "TG": 2.0,
        "TT": 1.6,
    },
    "AT/TA": {
        "AA": 1.2,
        "AC": 1.7,
        "AG": 0.2,
        "AT": 1.0,
        "CA": 1.7,
        "CC": 2.7,
        "CG": 1.0,
        "CT": 1.4,
        "GA": 0.2,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 1.4,
    },
    "AT/TG": {
        "AA": 1.8,
        "AC": 2.4,
        "AG": 0.6,
        "AT": 1.0,
        "CA": 2.6,
        "CC": 3.1,
        "CG": 1.0,
        "CT": 2.6,
        "GA": 0.8,
        "GC": 1.0,
        "GG": -0.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.4,
        "TG": 2.0,
        "TT": 2.1,
    },
    "CG/AT": {
        "AA": 1.1,
        "AC": 2.1,
        "AG": 0.8,
        "AT": 1.0,
        "CA": 1.7,
        "CC": 1.8,
        "CG": 1.0,
        "CT": 1.4,
        "GA": 0.5,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 0.6,
    },
    "CG/CG": {
        "AA": 0.6,
        "AC": 1.2,
        "AG": -0.5,
        "AT": 1.0,
        "CA": 1.6,
        "CC": 1.5,
        "CG": 1.0,
        "CT": 1.6,
        "GA": -0.1,
        "GC": 1.0,
        "GG": -1.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.0,
        "TG": 2.0,
        "TT": 0.3,
    },
    "CG/GC": {
        "AA": 0.9,
        "AC": 1.5,
        "AG": 0.1,
        "AT": 1.0,
        "CA": 1.5,
        "CC": 1.4,
        "CG": 1.0,
        "CT": 1.0,
        "GA": 0.1,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.0,
        "TG": 2.0,
        "TT": -0.2,
    },
    "CG/GT": {
        "AA": 1.9,
        "AC": 2.5,
        "AG": 1.1,
        "AT": 1.0,
        "CA": 2.5,
        "CC": 2.4,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 1.1,
        "GC": 1.0,
        "GG": 0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 0.8,
    },
    "CG/TA": {
        "AA": 1.0,
        "AC": 1.5,
        "AG": 0.1,
        "AT": 1.0,
        "CA": 1.7,
        "CC": 2.0,
        "CG": 1.0,
        "CT": 1.4,
        "GA": 0.3,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.0,
        "TG": 2.0,
        "TT": 0.6,
    },
    "CG/TG": {
        "AA": 1.6,
        "AC": 2.2,
        "AG": 0.5,
        "AT": 1.0,
        "CA": 2.6,
        "CC": 2.5,
        "CG": 1.0,
        "CT": 2.6,
        "GA": 0.9,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 1.3,
    },
    "GC/AT": {
        "AA": 0.9,
        "AC": 2.1,
        "AG": 0.5,
        "AT": 1.0,
        "CA": 1.4,
        "CC": 1.8,
        "CG": 1.0,
        "CT": 1.4,
        "GA": -0.1,
        "GC": 1.0,
        "GG": -0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 1.1,
    },
    "GC/CG": {
        "AA": 0.3,
        "AC": 1.3,
        "AG": -0.8,
        "AT": 1.0,
        "CA": 1.3,
        "CC": 1.6,
        "CG": 1.0,
        "CT": 1.6,
        "GA": -0.8,
        "GC": 1.0,
        "GG": -2.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.6,
        "TG": 2.0,
        "TT": 0.9,
    },
    "GC/GC": {
        "AA": 0.6,
        "AC": 1.6,
        "AG": -0.1,
        "AT": 1.0,
        "CA": 1.2,
        "CC": 1.5,
        "CG": 1.0,
        "CT": 1.0,
        "GA": -0.5,
        "GC": 1.0,
        "GG": -1.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.6,
        "TG": 2.0,
        "TT": 0.3,
    },
    "GC/GT": {
        "AA": 1.6,
        "AC": 2.6,
        "AG": 0.9,
        "AT": 1.0,
        "CA": 2.2,
        "CC": 2.5,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 0.5,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.6,
        "TG": 2.0,
        "TT": 1.3,
    },
    "GC/TA": {
        "AA": 0.8,
        "AC": 1.6,
        "AG": -0.2,
        "AT": 1.0,
        "CA": 1.4,
        "CC": 2.1,
        "CG": 1.0,
        "CT": 1.4,
        "GA": -0.4,
        "GC": 1.0,
        "GG": -1.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.6,
        "TG": 2.0,
        "TT": 1.1,
    },
    "GC/TG": {
        "AA": 1.3,
        "AC": 2.3,
        "AG": 0.2,
        "AT": 1.0,
        "CA": 2.3,
        "CC": 2.6,
        "CG": 1.0,
        "CT": 2.6,
        "GA": 0.2,
        "GC": 1.0,
        "GG": -1.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.6,
        "TG": 2.0,
        "TT": 1.9,
    },
    "GT/AT": {
        "AA": 1.9,
        "AC": 3.1,
        "AG": 1.5,
        "AT": 1.0,
        "CA": 2.4,
        "CC": 2.8,
        "CG": 1.0,
        "CT": 2.4,
        "GA": 0.9,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 3.0,
        "TG": 2.0,
        "TT": 2.1,
    },
    "GT/CG": {
        "AA": 1.3,
        "AC": 2.3,
        "AG": 0.2,
        "AT": 1.0,
        "CA": 2.3,
        "CC": 2.6,
        "CG": 1.0,
        "CT": 2.6,
        "GA": 0.2,
        "GC": 1.0,
        "GG": -1.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.6,
        "TG": 2.0,
        "TT": 1.9,
    },
    "GT/GC": {
        "AA": 1.6,
        "AC": 2.6,
        "AG": 0.9,
        "AT": 1.0,
        "CA": 2.2,
        "CC": 2.5,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 0.5,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.6,
        "TG": 2.0,
        "TT": 1.3,
    },
    "GT/GT": {
        "AA": 3.0,
        "AC": 3.0,
        "AG": 3.0,
        "AT": 3.0,
        "CA": 3.0,
        "CC": 3.0,
        "CG": 3.0,
        "CT": 3.0,
        "GA": 3.0,
        "GC": 3.0,
        "GG": 3.0,
        "GT": 3.0,
        "TA": 3.0,
        "TC": 3.0,
        "TG": 3.0,
        "TT": 3.0,
    },
    "GT/TA": {
        "AA": 1.8,
        "AC": 2.6,
        "AG": 0.8,
        "AT": 1.0,
        "CA": 2.4,
        "CC": 3.1,
        "CG": 1.0,
        "CT": 2.4,
        "GA": 0.6,
        "GC": 1.0,
        "GG": -0.2,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.6,
        "TG": 2.0,
        "TT": 2.1,
    },
    "GT/TG": {
        "AA": 3.0,
        "AC": 3.0,
        "AG": 3.0,
        "AT": 3.0,
        "CA": 3.0,
        "CC": 3.0,
        "CG": 3.0,
        "CT": 3.0,
        "GA": 3.0,
        "GC": 3.0,
        "GG": 3.0,
        "GT": 3.0,
        "TA": 3.0,
        "TC": 3.0,
        "TG": 3.0,
        "TT": 3.0,
    },
    "TA/AT": {
        "AA": 1.4,
        "AC": 2.3,
        "AG": 1.2,
        "AT": 1.0,
        "CA": 2.3,
        "CC": 2.1,
        "CG": 1.0,
        "CT": 1.7,
        "GA": 1.2,
        "GC": 1.0,
        "GG": 0.9,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.7,
        "TG": 2.0,
        "TT": 1.4,
    },
    "TA/CG": {
        "AA": 0.9,
        "AC": 1.4,
        "AG": -0.1,
        "AT": 1.0,
        "CA": 2.1,
        "CC": 1.8,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 0.5,
        "GC": 1.0,
        "GG": -0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 1.1,
    },
    "TA/GC": {
        "AA": 1.1,
        "AC": 1.7,
        "AG": 0.5,
        "AT": 1.0,
        "CA": 2.1,
        "CC": 1.8,
        "CG": 1.0,
        "CT": 1.4,
        "GA": 0.8,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 0.6,
    },
    "TA/GT": {
        "AA": 2.1,
        "AC": 2.7,
        "AG": 1.5,
        "AT": 1.0,
        "CA": 3.1,
        "CC": 2.8,
        "CG": 1.0,
        "CT": 2.4,
        "GA": 1.8,
        "GC": 1.0,
        "GG": 1.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.4,
        "TG": 2.0,
        "TT": 1.6,
    },
    "TA/TA": {
        "AA": 1.3,
        "AC": 1.7,
        "AG": 0.4,
        "AT": 1.0,
        "CA": 2.2,
        "CC": 2.4,
        "CG": 1.0,
        "CT": 1.7,
        "GA": 0.9,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 1.4,
        "TG": 2.0,
        "TT": 1.4,
    },
    "TA/TG": {
        "AA": 1.9,
        "AC": 2.4,
        "AG": 0.9,
        "AT": 1.0,
        "CA": 3.1,
        "CC": 2.8,
        "CG": 1.0,
        "CT": 3.0,
        "GA": 1.5,
        "GC": 1.0,
        "GG": 0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.4,
        "TG": 2.0,
        "TT": 2.1,
    },
    "TG/AT": {
        "AA": 2.1,
        "AC": 3.1,
        "AG": 1.8,
        "AT": 1.0,
        "CA": 2.7,
        "CC": 2.8,
        "CG": 1.0,
        "CT": 2.4,
        "GA": 1.5,
        "GC": 1.0,
        "GG": 1.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.4,
        "TG": 2.0,
        "TT": 1.6,
    },
    "TG/CG": {
        "AA": 1.6,
        "AC": 2.2,
        "AG": 0.5,
        "AT": 1.0,
        "CA": 2.6,
        "CC": 2.5,
        "CG": 1.0,
        "CT": 2.6,
        "GA": 0.9,
        "GC": 1.0,
        "GG": -0.3,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 1.3,
    },
    "TG/GC": {
        "AA": 1.9,
        "AC": 2.5,
        "AG": 1.1,
        "AT": 1.0,
        "CA": 2.5,
        "CC": 2.4,
        "CG": 1.0,
        "CT": 2.0,
        "GA": 1.1,
        "GC": 1.0,
        "GG": 0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 0.8,
    },
    "TG/GT": {
        "AA": 3.0,
        "AC": 3.0,
        "AG": 3.0,
        "AT": 3.0,
        "CA": 3.0,
        "CC": 3.0,
        "CG": 3.0,
        "CT": 3.0,
        "GA": 3.0,
        "GC": 3.0,
        "GG": 3.0,
        "GT": 3.0,
        "TA": 3.0,
        "TC": 3.0,
        "TG": 3.0,
        "TT": 3.0,
    },
    "TG/TA": {
        "AA": 2.0,
        "AC": 2.5,
        "AG": 1.1,
        "AT": 1.0,
        "CA": 2.7,
        "CC": 3.0,
        "CG": 1.0,
        "CT": 2.4,
        "GA": 1.3,
        "GC": 1.0,
        "GG": 0.7,
        "GT": 2.0,
        "TA": 1.0,
        "TC": 2.0,
        "TG": 2.0,
        "TT": 1.6,
    },
    "TG/TG": {
        "AA": 3.0,
        "AC": 3.0,
        "AG": 3.0,
        "AT": 3.0,
        "CA": 3.0,
        "CC": 3.0,
        "CG": 3.0,
        "CT": 3.0,
        "GA": 3.0,
        "GC": 3.0,
        "GG": 3.0,
        "GT": 3.0,
        "TA": 3.0,
        "TC": 3.0,
        "TG": 3.0,
        "TT": 3.0,
    },
}

_INT11_DH: dict[str, dict[str, float | None]] = {
    "AT/AT": {
        "AA": 12.9,
        "AC": 1.5,
        "AG": 10.7,
        "AT": 1.5,
        "CA": 10.4,
        "CC": 14.7,
        "CG": 10.4,
        "CT": -1.6,
        "GA": 1.8,
        "GC": 1.5,
        "GG": 1.6,
        "GT": 1.5,
        "TA": 10.4,
        "TC": 3.5,
        "TG": 10.4,
        "TT": -4.8,
    },
    "AT/CG": {
        "AA": 10.0,
        "AC": 10.9,
        "AG": 8.9,
        "AT": 10.9,
        "CA": 4.1,
        "CC": 12.8,
        "CG": 4.1,
        "CT": 1.3,
        "GA": 4.5,
        "GC": 10.9,
        "GG": -5.0,
        "GT": 10.9,
        "TA": 4.1,
        "TC": 11.6,
        "TG": 4.1,
        "TT": 0.5,
    },
    "AT/GC": {
        "AA": 5.0,
        "AC": 3.9,
        "AG": -6.3,
        "AT": 3.9,
        "CA": 5.3,
        "CC": 2.6,
        "CG": 5.3,
        "CT": -4.9,
        "GA": -2.2,
        "GC": 3.9,
        "GG": -9.9,
        "GT": 3.9,
        "TA": 5.3,
        "TC": -12.3,
        "TG": 5.3,
        "TT": -7.5,
    },
    "AT/GT": {
        "AA": 5.0,
        "AC": 3.9,
        "AG": -6.3,
        "AT": 3.9,
        "CA": 5.3,
        "CC": 2.6,
        "CG": 5.3,
        "CT": -4.9,
        "GA": -2.2,
        "GC": 3.9,
        "GG": -9.9,
        "GT": 3.9,
        "TA": 5.3,
        "TC": -12.3,
        "TG": 5.3,
        "TT": -7.5,
    },
    "AT/TA": {
        "AA": 4.9,
        "AC": 9.8,
        "AG": 7.1,
        "AT": 9.8,
        "CA": 9.8,
        "CC": 7.5,
        "CG": 9.8,
        "CT": 4.3,
        "GA": 7.1,
        "GC": 9.8,
        "GG": 4.1,
        "GT": 9.8,
        "TA": 9.8,
        "TC": 4.3,
        "TG": 9.8,
        "TT": 5.7,
    },
    "AT/TG": {
        "AA": 10.0,
        "AC": 10.9,
        "AG": 8.9,
        "AT": 10.9,
        "CA": 4.1,
        "CC": 12.8,
        "CG": 4.1,
        "CT": 1.3,
        "GA": 4.5,
        "GC": 10.9,
        "GG": -5.0,
        "GT": 10.9,
        "TA": 4.1,
        "TC": 11.6,
        "TG": 4.1,
        "TT": 0.5,
    },
    "CG/AT": {
        "AA": 10.0,
        "AC": 13.2,
        "AG": 11.3,
        "AT": 13.2,
        "CA": 9.6,
        "CC": 5.9,
        "CG": 9.6,
        "CT": 0.1,
        "GA": 6.9,
        "GC": 13.2,
        "GG": -1.1,
        "GT": 13.2,
        "TA": 9.6,
        "TC": 3.0,
        "TG": 9.6,
        "TT": 0.7,
    },
    "CG/CG": {
        "AA": -0.6,
        "AC": 4.6,
        "AG": 1.3,
        "AT": 4.6,
        "CA": 9.3,
        "CC": 7.5,
        "CG": 9.3,
        "CT": 1.9,
        "GA": 3.9,
        "GC": 4.6,
        "GG": -0.9,
        "GT": 4.6,
        "TA": 9.3,
        "TC": 0.5,
        "TG": 9.3,
        "TT": 2.1,
    },
    "CG/GC": {
        "AA": 4.0,
        "AC": 5.1,
        "AG": 6.1,
        "AT": 5.1,
        "CA": 5.1,
        "CC": 4.1,
        "CG": 5.1,
        "CT": 1.6,
        "GA": 6.1,
        "GC": 5.1,
        "GG": 0.1,
        "GT": 5.1,
        "TA": 5.1,
        "TC": 1.6,
        "TG": 5.1,
        "TT": -2.7,
    },
    "CG/GT": {
        "AA": 4.0,
        "AC": 5.1,
        "AG": 6.1,
        "AT": 5.1,
        "CA": 5.1,
        "CC": 4.1,
        "CG": 5.1,
        "CT": 1.6,
        "GA": 6.1,
        "GC": 5.1,
        "GG": 0.1,
        "GT": 5.1,
        "TA": 5.1,
        "TC": 1.6,
        "TG": 5.1,
        "TT": -2.7,
    },
    "CG/TA": {
        "AA": 5.0,
        "AC": 5.3,
        "AG": -2.2,
        "AT": 5.3,
        "CA": 3.9,
        "CC": 2.6,
        "CG": 3.9,
        "CT": -12.3,
        "GA": -6.3,
        "GC": 5.3,
        "GG": -9.9,
        "GT": 5.3,
        "TA": 3.9,
        "TC": -4.9,
        "TG": 3.9,
        "TT": -7.5,
    },
    "CG/TG": {
        "AA": -0.6,
        "AC": 4.6,
        "AG": 1.3,
        "AT": 4.6,
        "CA": 9.3,
        "CC": 7.5,
        "CG": 9.3,
        "CT": 1.9,
        "GA": 3.9,
        "GC": 4.6,
        "GG": -0.9,
        "GT": 4.6,
        "TA": 9.3,
        "TC": 0.5,
        "TG": 9.3,
        "TT": 2.1,
    },
    "GC/AT": {
        "AA": -1.4,
        "AC": 13.2,
        "AG": 11.3,
        "AT": 13.2,
        "CA": 4.2,
        "CC": 12.9,
        "CG": 4.2,
        "CT": -0.4,
        "GA": 6.9,
        "GC": 13.2,
        "GG": 3.4,
        "GT": 13.2,
        "TA": 4.2,
        "TC": 3.7,
        "TG": 4.2,
        "TT": -10.0,
    },
    "GC/CG": {
        "AA": 2.6,
        "AC": 3.9,
        "AG": 2.4,
        "AT": 3.9,
        "CA": 3.9,
        "CC": 6.0,
        "CG": 3.9,
        "CT": -0.6,
        "GA": 2.4,
        "GC": 3.9,
        "GG": -2.8,
        "GT": 3.9,
        "TA": 3.9,
        "TC": -0.6,
        "TG": 3.9,
        "TT": 0.3,
    },
    "GC/GC": {
        "AA": -0.6,
        "AC": 9.3,
        "AG": 3.9,
        "AT": 9.3,
        "CA": 4.6,
        "CC": 7.5,
        "CG": 4.6,
        "CT": 0.5,
        "GA": 1.3,
        "GC": 9.3,
        "GG": -0.9,
        "GT": 9.3,
        "TA": 4.6,
        "TC": 1.9,
        "TG": 4.6,
        "TT": 2.1,
    },
    "GC/GT": {
        "AA": -0.6,
        "AC": 9.3,
        "AG": 3.9,
        "AT": 9.3,
        "CA": 4.6,
        "CC": 7.5,
        "CG": 4.6,
        "CT": 0.5,
        "GA": 1.3,
        "GC": 9.3,
        "GG": -0.9,
        "GT": 9.3,
        "TA": 4.6,
        "TC": 1.9,
        "TG": 4.6,
        "TT": 2.1,
    },
    "GC/TA": {
        "AA": 10.0,
        "AC": 4.1,
        "AG": 4.5,
        "AT": 4.1,
        "CA": 10.9,
        "CC": 12.8,
        "CG": 10.9,
        "CT": 11.6,
        "GA": 8.9,
        "GC": 4.1,
        "GG": -5.0,
        "GT": 4.1,
        "TA": 10.9,
        "TC": 1.3,
        "TG": 10.9,
        "TT": 0.5,
    },
    "GC/TG": {
        "AA": 2.6,
        "AC": 3.9,
        "AG": 2.4,
        "AT": 3.9,
        "CA": 3.9,
        "CC": 6.0,
        "CG": 3.9,
        "CT": -0.6,
        "GA": 2.4,
        "GC": 3.9,
        "GG": -2.8,
        "GT": 3.9,
        "TA": 3.9,
        "TC": -0.6,
        "TG": 3.9,
        "TT": 0.3,
    },
    "GT/AT": {
        "AA": -1.4,
        "AC": 13.2,
        "AG": 11.3,
        "AT": 13.2,
        "CA": 4.2,
        "CC": 12.9,
        "CG": 4.2,
        "CT": -0.4,
        "GA": 6.9,
        "GC": 13.2,
        "GG": 3.4,
        "GT": 13.2,
        "TA": 4.2,
        "TC": 3.7,
        "TG": 4.2,
        "TT": -10.0,
    },
    "GT/CG": {
        "AA": 2.6,
        "AC": 3.9,
        "AG": 2.4,
        "AT": 3.9,
        "CA": 3.9,
        "CC": 6.0,
        "CG": 3.9,
        "CT": -0.6,
        "GA": 2.4,
        "GC": 3.9,
        "GG": -2.8,
        "GT": 3.9,
        "TA": 3.9,
        "TC": -0.6,
        "TG": 3.9,
        "TT": 0.3,
    },
    "GT/GC": {
        "AA": -0.6,
        "AC": 9.3,
        "AG": 3.9,
        "AT": 9.3,
        "CA": 4.6,
        "CC": 7.5,
        "CG": 4.6,
        "CT": 0.5,
        "GA": 1.3,
        "GC": 9.3,
        "GG": -0.9,
        "GT": 9.3,
        "TA": 4.6,
        "TC": 1.9,
        "TG": 4.6,
        "TT": 2.1,
    },
    "GT/GT": {
        "AA": -0.6,
        "AC": 9.3,
        "AG": 3.9,
        "AT": 9.3,
        "CA": 4.6,
        "CC": 7.5,
        "CG": 4.6,
        "CT": 0.5,
        "GA": 1.3,
        "GC": 9.3,
        "GG": -0.9,
        "GT": 9.3,
        "TA": 4.6,
        "TC": 1.9,
        "TG": 4.6,
        "TT": 2.1,
    },
    "GT/TA": {
        "AA": 10.0,
        "AC": 4.1,
        "AG": 4.5,
        "AT": 4.1,
        "CA": 10.9,
        "CC": 12.8,
        "CG": 10.9,
        "CT": 11.6,
        "GA": 8.9,
        "GC": 4.1,
        "GG": -5.0,
        "GT": 4.1,
        "TA": 10.9,
        "TC": 1.3,
        "TG": 10.9,
        "TT": 0.5,
    },
    "GT/TG": {
        "AA": 2.6,
        "AC": 3.9,
        "AG": 2.4,
        "AT": 3.9,
        "CA": 3.9,
        "CC": 6.0,
        "CG": 3.9,
        "CT": -0.6,
        "GA": 2.4,
        "GC": 3.9,
        "GG": -2.8,
        "GT": 3.9,
        "TA": 3.9,
        "TC": -0.6,
        "TG": 3.9,
        "TT": 0.3,
    },
    "TA/AT": {
        "AA": 12.1,
        "AC": 14.3,
        "AG": 11.9,
        "AT": 14.3,
        "CA": 14.3,
        "CC": 17.4,
        "CG": 14.3,
        "CT": 2.5,
        "GA": 11.9,
        "GC": 14.3,
        "GG": 8.4,
        "GT": 14.3,
        "TA": 14.3,
        "TC": 2.5,
        "TG": 14.3,
        "TT": 6.2,
    },
    "TA/CG": {
        "AA": -1.4,
        "AC": 4.2,
        "AG": 6.9,
        "AT": 4.2,
        "CA": 13.2,
        "CC": 12.9,
        "CG": 13.2,
        "CT": 3.7,
        "GA": 11.3,
        "GC": 4.2,
        "GG": 3.4,
        "GT": 4.2,
        "TA": 13.2,
        "TC": -0.4,
        "TG": 13.2,
        "TT": -10.0,
    },
    "TA/GC": {
        "AA": 10.0,
        "AC": 9.6,
        "AG": 6.9,
        "AT": 9.6,
        "CA": 13.2,
        "CC": 5.9,
        "CG": 13.2,
        "CT": 3.0,
        "GA": 11.3,
        "GC": 9.6,
        "GG": -1.1,
        "GT": 9.6,
        "TA": 13.2,
        "TC": 0.1,
        "TG": 13.2,
        "TT": 0.7,
    },
    "TA/GT": {
        "AA": 10.0,
        "AC": 9.6,
        "AG": 6.9,
        "AT": 9.6,
        "CA": 13.2,
        "CC": 5.9,
        "CG": 13.2,
        "CT": 3.0,
        "GA": 11.3,
        "GC": 9.6,
        "GG": -1.1,
        "GT": 9.6,
        "TA": 13.2,
        "TC": 0.1,
        "TG": 13.2,
        "TT": 0.7,
    },
    "TA/TA": {
        "AA": 12.9,
        "AC": 10.4,
        "AG": 1.8,
        "AT": 10.4,
        "CA": 1.5,
        "CC": 14.7,
        "CG": 1.5,
        "CT": 3.5,
        "GA": 10.7,
        "GC": 10.4,
        "GG": 1.6,
        "GT": 10.4,
        "TA": 1.5,
        "TC": -1.6,
        "TG": 1.5,
        "TT": -4.8,
    },
    "TA/TG": {
        "AA": -1.4,
        "AC": 4.2,
        "AG": 6.9,
        "AT": 4.2,
        "CA": 13.2,
        "CC": 12.9,
        "CG": 13.2,
        "CT": 3.7,
        "GA": 11.3,
        "GC": 4.2,
        "GG": 3.4,
        "GT": 4.2,
        "TA": 13.2,
        "TC": -0.4,
        "TG": 13.2,
        "TT": -10.0,
    },
    "TG/AT": {
        "AA": 10.0,
        "AC": 13.2,
        "AG": 11.3,
        "AT": 13.2,
        "CA": 9.6,
        "CC": 5.9,
        "CG": 9.6,
        "CT": 0.1,
        "GA": 6.9,
        "GC": 13.2,
        "GG": -1.1,
        "GT": 13.2,
        "TA": 9.6,
        "TC": 3.0,
        "TG": 9.6,
        "TT": 0.7,
    },
    "TG/CG": {
        "AA": -0.6,
        "AC": 4.6,
        "AG": 1.3,
        "AT": 4.6,
        "CA": 9.3,
        "CC": 7.5,
        "CG": 9.3,
        "CT": 1.9,
        "GA": 3.9,
        "GC": 4.6,
        "GG": -0.9,
        "GT": 4.6,
        "TA": 9.3,
        "TC": 0.5,
        "TG": 9.3,
        "TT": 2.1,
    },
    "TG/GC": {
        "AA": 4.0,
        "AC": 5.1,
        "AG": 6.1,
        "AT": 5.1,
        "CA": 5.1,
        "CC": 4.1,
        "CG": 5.1,
        "CT": 1.6,
        "GA": 6.1,
        "GC": 5.1,
        "GG": 0.1,
        "GT": 5.1,
        "TA": 5.1,
        "TC": 1.6,
        "TG": 5.1,
        "TT": -2.7,
    },
    "TG/GT": {
        "AA": 4.0,
        "AC": 5.1,
        "AG": 6.1,
        "AT": 5.1,
        "CA": 5.1,
        "CC": 4.1,
        "CG": 5.1,
        "CT": 1.6,
        "GA": 6.1,
        "GC": 5.1,
        "GG": 0.1,
        "GT": 5.1,
        "TA": 5.1,
        "TC": 1.6,
        "TG": 5.1,
        "TT": -2.7,
    },
    "TG/TA": {
        "AA": 5.0,
        "AC": 5.3,
        "AG": -2.2,
        "AT": 5.3,
        "CA": 3.9,
        "CC": 2.6,
        "CG": 3.9,
        "CT": -12.3,
        "GA": -6.3,
        "GC": 5.3,
        "GG": -9.9,
        "GT": 5.3,
        "TA": 3.9,
        "TC": -4.9,
        "TG": 3.9,
        "TT": -7.5,
    },
    "TG/TG": {
        "AA": -0.6,
        "AC": 4.6,
        "AG": 1.3,
        "AT": 4.6,
        "CA": 9.3,
        "CC": 7.5,
        "CG": 9.3,
        "CT": 1.9,
        "GA": 3.9,
        "GC": 4.6,
        "GG": -0.9,
        "GT": 4.6,
        "TA": 9.3,
        "TC": 0.5,
        "TG": 9.3,
        "TT": 2.1,
    },
}

_TERMINAL_MISMATCH_DG37: dict[str, dict[str, float | None]] = {
    "AT": {
        "AA": -0.7,
        "AC": -0.3,
        "AG": -0.5,
        "AT": 0.0,
        "CA": -0.6,
        "CC": -0.2,
        "CG": 0.0,
        "CT": -0.3,
        "GA": -0.6,
        "GC": 0.0,
        "GG": -0.4,
        "GT": -0.5,
        "TA": 0.0,
        "TC": -0.3,
        "TG": -0.5,
        "TT": -0.4,
    },
    "CG": {
        "AA": -1.0,
        "AC": -0.8,
        "AG": -0.9,
        "AT": 0.0,
        "CA": -0.8,
        "CC": -0.5,
        "CG": 0.0,
        "CT": -0.7,
        "GA": -1.0,
        "GC": 0.0,
        "GG": -0.9,
        "GT": -1.0,
        "TA": 0.0,
        "TC": -0.6,
        "TG": -0.9,
        "TT": -0.9,
    },
    "GC": {
        "AA": -1.0,
        "AC": -0.7,
        "AG": -0.8,
        "AT": 0.0,
        "CA": -1.0,
        "CC": -0.6,
        "CG": 0.0,
        "CT": -0.7,
        "GA": -1.0,
        "GC": 0.0,
        "GG": -1.0,
        "GT": -0.8,
        "TA": 0.0,
        "TC": -0.6,
        "TG": -0.9,
        "TT": -0.9,
    },
    "GT": {
        "AA": -0.5,
        "AC": -0.2,
        "AG": -0.5,
        "AT": -0.6,
        "CA": -0.2,
        "CC": -0.2,
        "CG": -0.9,
        "CT": -0.2,
        "GA": -0.5,
        "GC": -0.9,
        "GG": -0.5,
        "GT": -0.5,
        "TA": -0.5,
        "TC": -0.2,
        "TG": -0.5,
        "TT": -0.2,
    },
    "TA": {
        "AA": -0.6,
        "AC": -0.4,
        "AG": -0.5,
        "AT": 0.0,
        "CA": -0.5,
        "CC": -0.2,
        "CG": 0.0,
        "CT": -0.5,
        "GA": -0.6,
        "GC": -0.0,
        "GG": -0.4,
        "GT": -0.5,
        "TA": 0.0,
        "TC": -0.3,
        "TG": -0.6,
        "TT": -0.3,
    },
    "TG": {
        "AA": -0.5,
        "AC": -0.2,
        "AG": -0.5,
        "AT": -0.5,
        "CA": -0.2,
        "CC": -0.2,
        "CG": -0.8,
        "CT": -0.2,
        "GA": -0.5,
        "GC": -1.0,
        "GG": -0.5,
        "GT": -0.5,
        "TA": -0.5,
        "TC": -0.2,
        "TG": -0.5,
        "TT": -0.2,
    },
}

_TERMINAL_MISMATCH_DH: dict[str, dict[str, float | None]] = {
    "AT": {
        "AA": 4.0,
        "AC": 5.2,
        "AG": 4.7,
        "AT": 4.3,
        "CA": 4.0,
        "CC": 5.2,
        "CG": 4.7,
        "CT": 4.3,
        "GA": 4.0,
        "GC": 5.2,
        "GG": 4.7,
        "GT": 4.3,
        "TA": 4.0,
        "TC": 5.2,
        "TG": 4.7,
        "TT": 4.3,
    },
    "CG": {
        "AA": -4.6,
        "AC": -4.2,
        "AG": -4.7,
        "AT": -5.0,
        "CA": -4.6,
        "CC": -4.2,
        "CG": -4.7,
        "CT": -5.0,
        "GA": -4.6,
        "GC": -4.2,
        "GG": -4.7,
        "GT": -5.0,
        "TA": -4.6,
        "TC": -4.2,
        "TG": -4.7,
        "TT": -5.0,
    },
    "GC": {
        "AA": -6.2,
        "AC": -3.3,
        "AG": -5.1,
        "AT": -0.9,
        "CA": -6.2,
        "CC": -3.3,
        "CG": -5.1,
        "CT": -0.9,
        "GA": -6.2,
        "GC": -3.3,
        "GG": -5.1,
        "GT": -0.9,
        "TA": -6.2,
        "TC": -3.3,
        "TG": -5.1,
        "TT": -0.9,
    },
    "GT": {
        "AA": -3.0,
        "AC": -0.1,
        "AG": -1.9,
        "AT": 2.3,
        "CA": -3.0,
        "CC": -0.1,
        "CG": -1.9,
        "CT": 2.3,
        "GA": -3.0,
        "GC": -0.1,
        "GG": -1.9,
        "GT": 2.3,
        "TA": -3.0,
        "TC": -0.1,
        "TG": -1.9,
        "TT": 2.3,
    },
    "TA": {
        "AA": 2.3,
        "AC": 1.7,
        "AG": 0.7,
        "AT": -5.8,
        "CA": 2.3,
        "CC": 1.7,
        "CG": 0.7,
        "CT": -5.8,
        "GA": 2.3,
        "GC": 1.7,
        "GG": 0.7,
        "GT": -5.8,
        "TA": 2.3,
        "TC": 1.7,
        "TG": 0.7,
        "TT": -5.8,
    },
    "TG": {
        "AA": 2.3,
        "AC": 1.7,
        "AG": 0.7,
        "AT": -5.8,
        "CA": 2.3,
        "CC": 1.7,
        "CG": 0.7,
        "CT": -5.8,
        "GA": 2.3,
        "GC": 1.7,
        "GG": 0.7,
        "GT": -5.8,
        "TA": 2.3,
        "TC": 1.7,
        "TG": 0.7,
        "TT": -5.8,
    },
}

_DANGLE_DG37: dict[str, dict[str, float | None]] = {
    "AT|bottom_5p": {"A": -0.5, "C": -0.2, "G": -0.5, "T": -0.3},
    "AT|top_3p": {"A": -0.2, "C": -0.2, "G": -0.2, "T": -0.2},
    "CG|bottom_5p": {"A": -0.8, "C": -0.5, "G": -0.8, "T": -0.7},
    "CG|top_3p": {"A": -1.1, "C": -0.4, "G": -1.1, "T": -0.8},
    "GC|bottom_5p": {"A": -0.9, "C": -0.5, "G": -0.7, "T": -0.6},
    "GC|top_3p": {"A": -0.4, "C": -0.2, "G": -0.4, "T": -0.3},
    "GT|bottom_5p": {"A": -0.1, "C": -0.1, "G": -0.1, "T": -0.1},
    "GT|top_3p": {"A": -0.2, "C": -0.2, "G": -0.2, "T": -0.2},
    "TA|bottom_5p": {"A": -0.5, "C": -0.2, "G": -0.6, "T": -0.3},
    "TA|top_3p": {"A": -0.2, "C": -0.2, "G": -0.2, "T": -0.2},
    "TG|bottom_5p": {"A": -0.1, "C": -0.1, "G": -0.1, "T": -0.1},
    "TG|top_3p": {"A": -0.2, "C": -0.2, "G": -0.2, "T": -0.2},
}

_DANGLE_DH: dict[str, dict[str, float | None]] = {
    "AT|bottom_5p": {"A": 0.8, "C": 2.0, "G": 1.5, "T": 1.1},
    "AT|top_3p": {"A": 4.3, "C": 6.3, "G": 2.5, "T": 2.6},
    "CG|bottom_5p": {"A": -4.6, "C": -4.2, "G": -4.7, "T": -5.0},
    "CG|top_3p": {"A": -8.0, "C": -2.5, "G": -3.1, "T": -3.9},
    "GC|bottom_5p": {"A": -6.2, "C": -3.3, "G": -5.1, "T": -0.9},
    "GC|top_3p": {"A": -1.9, "C": -0.2, "G": -4.0, "T": -5.5},
    "GT|bottom_5p": {"A": -6.2, "C": -3.3, "G": -5.1, "T": -0.9},
    "GT|top_3p": {"A": -1.9, "C": -0.2, "G": -4.0, "T": -5.5},
    "TA|bottom_5p": {"A": -0.9, "C": -0.6, "G": -2.5, "T": -9.0},
    "TA|top_3p": {"A": -5.3, "C": 2.0, "G": -4.8, "T": 0.9},
    "TG|bottom_5p": {"A": -4.6, "C": -4.2, "G": -4.7, "T": -5.0},
    "TG|top_3p": {"A": -8.0, "C": -2.5, "G": -3.1, "T": -3.9},
}


def _get_nn_params(rna: bool = False, hybrid: bool = False) -> dict[str, tuple[float, float]]:
    """Return the appropriate NN parameter table.

    Parameters
    ----------
    rna : bool
        RNA/RNA table (Turner 2004).
    hybrid : bool
        RNA/DNA hybrid table (Banerjee et al. 2020). Takes precedence
        over ``rna`` when True.
    """
    if hybrid:
        return _HYBRID_NN
    return _RNA_NN if rna else _DNA_NN


def _get_init_params(
    rna: bool = False,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Return (init_GC, init_AT, symmetry) parameters (DNA model only)."""
    return _DNA_INIT_GC, _DNA_INIT_AT, _DNA_SYMMETRY


def _is_self_complementary(sequence: str, rna: bool = False) -> bool:
    """Check if a sequence is self-complementary (palindromic)."""
    seq = sequence.upper()
    if rna:
        seq = seq.replace("T", "U")
    comp_map = str.maketrans("ACGU", "UGCA") if rna else str.maketrans("ACGT", "TGCA")
    complement = seq.translate(comp_map)
    return seq == complement[::-1]


@lru_cache(maxsize=65536)
def calculate_nn_thermodynamics(
    sequence: str,
    rna: bool = False,
    hybrid: bool = False,
) -> tuple[float, float]:
    """Calculate total ΔH and ΔS using the nearest-neighbor model.

    Sums nearest-neighbor stacking enthalpies and entropies, plus
    initiation and symmetry corrections. Results are memoized per
    (sequence, rna, hybrid) triple — batch scoring of sliding-window
    candidates calls this repeatedly with overlapping sequences, so
    caching eliminates redundant computation.

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence. DNA alphabet for DNA/DNA and hybrid
        models; DNA or RNA alphabet for the RNA model (T → U converted).
    rna : bool
        If True, use RNA/RNA NN parameters (Turner 2004).
        If False, use DNA/DNA NN parameters (SantaLucia 1998).
    hybrid : bool
        If True, use RNA/DNA hybrid NN parameters (Banerjee et al.
        2020) for an ASO:target-RNA duplex. Takes precedence over
        ``rna``. The table is keyed on the RNA target strand, so the
        RNA strand 5'→3' (the reverse complement of the ASO, T → U)
        is what the stacking sum runs over.

    Returns
    -------
    tuple[float, float]
        (Total ΔH in kcal/mol, Total ΔS in cal/(mol·K))

    Raises
    ------
    ThermodynamicsError
        If the sequence is too short or contains invalid NN pairs.
    """
    seq = sequence.upper()
    if rna and not hybrid:
        seq = seq.replace("T", "U")

    if len(seq) < 2:
        raise ThermodynamicsError(f"Sequence too short for NN calculation: {len(seq)} nt")

    nn_params = _get_nn_params(rna, hybrid)

    total_h = 0.0
    total_s = 0.0

    if hybrid:
        # RNA/DNA hybrid duplex. The published parameters are defined on
        # the RNA (target) strand, so convert the DNA ASO 5'→3' sequence
        # to its antiparallel RNA target strand read 5'→3' (reverse
        # complement, T→U) and sum stacks over that strand.
        rna_strand = seq.translate(_DNA_COMPLEMENT)[::-1].replace("T", "U")
        for i in range(len(rna_strand) - 1):
            nn = rna_strand[i : i + 2]
            if nn not in nn_params:
                raise ThermodynamicsError(
                    f"Invalid nearest-neighbor pair: {nn} in sequence {rna_strand}"
                )
            h, s = nn_params[nn]
            total_h += h
            total_s += s
        # Single helix-initiation term; no terminal or symmetry corrections
        # (see the table notes above).
        total_h += _HYBRID_INIT[0]
        total_s += _HYBRID_INIT[1]
        return total_h, total_s

    # Sum NN stacking parameters
    for i in range(len(seq) - 1):
        nn = seq[i : i + 2]
        if nn not in nn_params:
            raise ThermodynamicsError(f"Invalid nearest-neighbor pair: {nn} in sequence {seq}")
        h, s = nn_params[nn]
        total_h += h
        total_s += s

    if rna:
        # Turner 2004 initiation model (NNDB):
        # intermolecular initiation (once) + AU end penalty (per terminal
        # A·U base pair) + symmetry correction (self-complementary duplexes)
        total_h += _RNA_INIT_INTERMOLECULAR[0]
        total_s += _RNA_INIT_INTERMOLECULAR[1]
        for terminal_base in (seq[0], seq[-1]):
            if terminal_base in ("A", "U"):
                total_h += _RNA_AU_END[0]
                total_s += _RNA_AU_END[1]
        if _is_self_complementary(seq, rna=True):
            total_h += _RNA_SYMMETRY[0]
            total_s += _RNA_SYMMETRY[1]
    else:
        # SantaLucia 1998 initiation model: per-terminal base pair
        # corrections (init GC or AT) + symmetry correction
        init_gc, init_at, symmetry = _get_init_params(rna=False)

        # 5' terminal
        first_base = seq[0]
        if first_base in ("G", "C"):
            total_h += init_gc[0]
            total_s += init_gc[1]
        else:
            total_h += init_at[0]
            total_s += init_at[1]

        # 3' terminal
        last_base = seq[-1]
        if last_base in ("G", "C"):
            total_h += init_gc[0]
            total_s += init_gc[1]
        else:
            total_h += init_at[0]
            total_s += init_at[1]

        # Symmetry correction
        if _is_self_complementary(seq, rna=False):
            total_h += symmetry[0]
            total_s += symmetry[1]

    return total_h, total_s


# ---------------------------------------------------------------------------
# Divalent-ion (Mg2+) salt correction
# Owczarzy et al. (2008) Biochemistry 47:5336-5353, method 7: a correction
# added to 1/Tm (Kelvin) that accounts for both Mg2+ and the monovalent
# background. Cross-validated against the reference implementation in
# Biopython's Bio.SeqUtils.MeltingTemp.salt_correction (method 7).
# ---------------------------------------------------------------------------

# Mg:dNTP dissociation constant (M^-1) used to compute free Mg2+ when dNTPs
# are present (Owczarzy 2008).
_MG_DNTP_KA = 3.0e4


def calculate_mg_correction(
    sequence: str,
    na_conc: float,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
) -> float:
    """Compute the Owczarzy 2008 divalent-ion correction for 1/Tm.

    The returned term ``corr`` is added to ``1/Tm`` with Tm in Kelvin:

        Tm(corrected) = 1 / (1 / Tm(1 M standard state) + corr)

    The correction replaces (not supplements) the monovalent ΔS salt
    correction: its low-Mg branch is the Owczarzy 2004 monovalent
    formula, and the Mg branches fold the monovalent background into
    the empirical coefficients via the sqrt([Mg2+])/[monovalent] ratio.

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence (DNA alphabet); used for the GC
        fraction and length terms.
    na_conc : float
        Total monovalent cation concentration in M (Na+ + K+ + Tris/2).
        Must be > 0.
    mg_conc : float
        Total magnesium concentration in M (0 disables the correction).
    dntp_conc : float
        Total dNTP concentration in M. dNTPs chelate Mg2+; free Mg2+ is
        computed with a 1:1 binding model (Ka = 3e4 M^-1). When dNTPs
        overwhelm Mg2+, free Mg2+ becomes negligible and the correction
        reduces to the monovalent (2004) branch.

    Returns
    -------
    float
        Correction term in K^-1 (0.0 when ``mg_conc`` <= 0).

    Raises
    ------
    ThermodynamicsError
        If ``na_conc`` is not positive when Mg2+ is present.
    """
    import math

    if mg_conc <= 0:
        return 0.0
    if na_conc <= 0:
        raise ThermodynamicsError(
            "na_conc must be > 0 when Mg2+ is present (Owczarzy correction "
            "requires a monovalent background)"
        )

    seq = sequence.upper()
    gc_fraction = (seq.count("G") + seq.count("C")) / len(seq)
    mon = na_conc  # M
    mg = mg_conc  # M

    if dntp_conc > 0:
        # Free Mg2+ after dNTP chelation (quadratic solution of the
        # 1:1 binding equilibrium).
        x = _MG_DNTP_KA * dntp_conc - _MG_DNTP_KA * mg + 1.0
        mg = (-x + math.sqrt(x * x + 4.0 * _MG_DNTP_KA * mg)) / (2.0 * _MG_DNTP_KA)

    # Decision tree on R = sqrt([free Mg2+]) / [monovalent]
    ratio = math.sqrt(mg) / mon
    if ratio < 0.22:
        # Monovalent-dominated regime: Owczarzy 2004 monovalent formula
        # (Biochemistry 43:3537-3554).
        return (4.29 * gc_fraction - 3.95) * 1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2

    a, b, c, d = 3.92, -0.911, 6.26, 1.42
    e, f, g = -48.2, 52.5, 8.31
    if ratio < 6.0:
        # Intermediate regime: monovalent-adjusted coefficients.
        a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
        d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
        g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    # Else (ratio >= 6.0): Mg2+-dominated regime with the default coefficients.

    corr = (
        a
        + b * math.log(mg)
        + gc_fraction * (c + d * math.log(mg))
        + (1.0 / (2.0 * (len(seq) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5
    return corr


# ---------------------------------------------------------------------------
# Melting temperature (Tm) calculation
# ---------------------------------------------------------------------------


def calculate_tm(
    sequence: str,
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
    rna: bool = False,
    salt_correction: str = "santalucia",
    hybrid: bool = False,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
) -> float:
    """Calculate the melting temperature (Tm) using the NN model.

    Tm = (ΔH * 1000) / (ΔS + R * ln(C_T / x)) - 273.15

    where:
    - ΔH = total enthalpy (kcal/mol)
    - ΔS = total entropy (cal/(mol·K))
    - R = 1.987 cal/(mol·K)
    - C_T = total strand concentration (M)
    - x = 4 for non-self-complementary, 1 for self-complementary

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence (DNA or RNA alphabet).
    na_conc : float
        Sodium concentration in M (default 0.1 M = 100 mM). When Mg2+
        is present this is the total monovalent background
        (Na+ + K+ + Tris/2).
    oligo_conc : float
        Total oligo strand concentration in M (default 1 µM).
    rna : bool
        If True, use RNA/RNA NN parameters.
    salt_correction : str
        Salt correction method: "santalucia" (1998) or "log10" (simple).
        Ignored when ``mg_conc`` > 0 (the Owczarzy 2008 correction
        handles monovalent + divalent ions together).
    hybrid : bool
        If True, use RNA/DNA hybrid NN parameters. The concentration
        factor x is 4 (hybrid duplexes are never self-complementary).
    mg_conc : float
        Magnesium concentration in M (0 disables Mg2+ correction).
    dntp_conc : float
        dNTP concentration in M (chelates Mg2+; see
        ``calculate_mg_correction``).

    Returns
    -------
    float
        Melting temperature in °C.

    Raises
    ------
    ThermodynamicsError
        If NN calculation fails or salt correction is unknown.

    Notes
    -----
    The SantaLucia ΔS salt correction was derived for DNA/DNA duplexes;
    applying it to hybrid or RNA parameters is an extrapolation kept for
    consistency until a hybrid-specific salt model is adopted (see the
    RNA/DNA hybrid table notes).
    """
    import math

    dh, ds = calculate_nn_thermodynamics(sequence, rna=rna, hybrid=hybrid)

    # Symmetry factor: 1 for self-complementary, 4 otherwise
    seq = sequence.upper()
    if rna and not hybrid:
        seq = seq.replace("T", "U")
    is_self_comp = False if hybrid else _is_self_complementary(seq, rna=rna)
    x = 1 if is_self_comp else 4

    # Tm at 1 M NaCl
    # ΔH in kcal/mol → * 1000 to get cal/mol
    # ΔS in cal/(mol·K)
    # R * ln(C_T / x) in cal/(mol·K)
    denominator = ds + R * math.log(oligo_conc / x)
    if denominator == 0:
        raise ThermodynamicsError("Denominator is zero in Tm calculation")

    tm_1m = (dh * 1000) / denominator - 273.15

    # Mg2+ present: the Owczarzy 2008 correction replaces the monovalent
    # salt corrections entirely (it corrects 1/Tm from the 1 M standard
    # state and folds the monovalent background into its own branches).
    if mg_conc > 0:
        corr = calculate_mg_correction(seq, na_conc, mg_conc, dntp_conc)
        return 1.0 / (1.0 / (tm_1m + 273.15) + corr) - 273.15

    # Apply salt correction
    if salt_correction == "santalucia":
        # SantaLucia 1998: adjust ΔS for salt concentration
        # ΔS_corrected = ΔS + 0.368 * (N-1) * ln([Na+])
        # Then recompute Tm with corrected ΔS
        n_phosphates = len(seq)
        ds_corrected = ds + 0.368 * (n_phosphates - 1) * math.log(na_conc)
        if ds_corrected + R * math.log(oligo_conc / x) == 0:
            raise ThermodynamicsError("Denominator is zero in Tm calculation (salt-corrected)")
        tm = (dh * 1000) / (ds_corrected + R * math.log(oligo_conc / x)) - 273.15
    elif salt_correction == "log10":
        # Simple: Tm(salt) = Tm(1M) + 16.6 * log10([Na+])
        tm = tm_1m + 16.6 * math.log10(na_conc)
    else:
        raise ThermodynamicsError(
            f"Unknown salt correction method: {salt_correction}. Use 'santalucia' or 'log10'."
        )

    return tm


def calculate_dg(
    sequence: str,
    temperature: float = DEFAULT_TEMPERATURE_C,
    rna: bool = False,
    hybrid: bool = False,
) -> float:
    """Calculate ΔG at a given temperature using the NN model.

    ΔG = ΔH - T * ΔS / 1000

    where T is in Kelvin and ΔS is converted from cal to kcal.

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence.
    temperature : float
        Temperature in °C (default 37 °C).
    rna : bool
        If True, use RNA/RNA NN parameters.
    hybrid : bool
        If True, use RNA/DNA hybrid NN parameters (takes precedence
        over ``rna``).

    Returns
    -------
    float
        ΔG in kcal/mol (negative = favorable).
    """
    dh, ds = calculate_nn_thermodynamics(sequence, rna=rna, hybrid=hybrid)
    t_kelvin = temperature + 273.15
    # ΔS is in cal/(mol·K), ΔH is in kcal/mol → convert ΔS to kcal
    dg = dh - (t_kelvin * ds / 1000.0)
    return dg


# ---------------------------------------------------------------------------
# Chemistry-aware Tm calculation and scoring
# ---------------------------------------------------------------------------

# Chemistry types that use RNA NN parameters (target is RNA, ASO is RNA-like)
_RNA_CHEMISTRIES = {"2ome", "moe"}
# Chemistry types that form an RNA/DNA hybrid duplex with the RNA target
# (DNA-containing ASOs: DNA gap gapmers, mixed DNA/modified mixmers, and the
# DNA-like morpholino backbone). See _HYBRID_NN for the parameter source.
_HYBRID_CHEMISTRIES = {"gapmer", "mixmer", "morpholino", "lna"}

# Per-position LNA affinity adjustment (heuristic, from the LNA literature:
# each LNA raises duplex Tm by ~2-8 degC and stabilizes the duplex by
# roughly -1 to -2 kcal/mol; we use mid-range values per LNA site).
_LNA_TM_PER_SITE = 3.5
_LNA_DG_PER_SITE = -1.4
_LNA_MOD_LABEL = "LNA"


def get_chemistry_end_gc_bounds(chemistry: str | None) -> tuple[float, float]:
    """Get per-chemistry end-stability GC bounds from the presets.

    Parameters
    ----------
    chemistry : str or None
        Chemistry name. None or an unknown name falls back to the
        chemistry-agnostic default bounds.

    Returns
    -------
    tuple[float, float]
        ``(min_end_gc, max_end_gc)`` for the end-stability filter.
    """
    default = (0.0, 0.67)
    if not chemistry:
        return default
    try:
        presets = _load_chemistry_presets()
    except (OSError, yaml.YAMLError) as e:
        logger.warning("Could not load chemistry presets; using default bounds: %s", e)
        return default
    preset = presets.get(chemistry, {})
    return (
        float(preset.get("end_gc_min", default[0])),
        float(preset.get("end_gc_max", default[1])),
    )


def get_chemistry_tm_adjustment(chemistry: str) -> float:
    """Get the Tm adjustment for a given chemistry from presets.

    Parameters
    ----------
    chemistry : str
        Chemistry name: gapmer, mixmer, 2ome, moe, morpholino.

    Returns
    -------
    float
        Tm adjustment in °C (added to computed Tm).
    """
    try:
        presets = _load_chemistry_presets()
    except (OSError, yaml.YAMLError) as e:
        logger.warning("Could not load chemistry presets; using 0.0 adjustment: %s", e)
        return 0.0

    preset = presets.get(chemistry, {})
    return preset.get("tm_adjustment", 0.0)


def calculate_tm_with_chemistry(
    sequence: str,
    chemistry: str = "gapmer",
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    modification_pattern: list[str] | None = None,
) -> float:
    """Calculate Tm with chemistry-specific NN parameters and adjustment.

    Selects the appropriate NN parameter set based on chemistry:
    - 2'OMe, MOE: RNA/RNA NN parameters (ASO is RNA-like, binds RNA target)
    - Gapmer, mixmer, morpholino: RNA/DNA hybrid NN parameters (the ASO
      carries a DNA-like backbone and pairs with the RNA target)
    - RNA/DNA hybrid calculations use the Banerjee et al. (2020) 1 M NaCl
      parameter set; see _HYBRID_NN.

    Then applies the chemistry-specific Tm adjustment from presets.

    Parameters
    ----------
    sequence : str
        ASO sequence in DNA alphabet (ACGT, will be converted to RNA if needed).
    chemistry : str
        Chemistry name (default "gapmer").
    na_conc : float
        Sodium concentration in M (default 0.1).
    oligo_conc : float
        Oligo concentration in M (default 1 µM).
    mg_conc : float
        Magnesium concentration in M (0 disables Mg2+ correction).
    dntp_conc : float
        dNTP concentration in M (chelates Mg2+).
    modification_pattern : list of str or None
        Per-nucleotide modification labels; with the "lna" chemistry each
        LNA position raises the Tm by ``_LNA_TM_PER_SITE``.

    Returns
    -------
    float
        Adjusted Tm in °C.
    """
    use_rna = chemistry in _RNA_CHEMISTRIES
    use_hybrid = chemistry in _HYBRID_CHEMISTRIES
    base_tm = calculate_tm(
        sequence,
        na_conc=na_conc,
        oligo_conc=oligo_conc,
        rna=use_rna,
        hybrid=use_hybrid,
        mg_conc=mg_conc,
        dntp_conc=dntp_conc,
    )
    adjustment = get_chemistry_tm_adjustment(chemistry)
    if chemistry == "lna" and modification_pattern:
        n_lna = sum(1 for mod in modification_pattern if mod == _LNA_MOD_LABEL)
        adjustment += _LNA_TM_PER_SITE * n_lna
    adjusted_tm = base_tm + adjustment
    logger.debug(
        "Tm for %s (%s): base=%.2f, adj=%+.2f, final=%.2f",
        sequence[:10],
        chemistry,
        base_tm,
        adjustment,
        adjusted_tm,
    )
    return adjusted_tm


def score_tm(
    tm: float,
    tm_optimal: float = 52.0,
    tm_min: float = 40.0,
    tm_max: float = 65.0,
) -> float:
    """Score a Tm value based on distance from optimal.

    Returns a score in [0.0, 1.0]:
    - 1.0 at tm_optimal
    - Linear decrease to 0.5 at tm_min or tm_max
    - 0.0 outside [tm_min, tm_max]

    Parameters
    ----------
    tm : float
        Computed Tm in °C.
    tm_optimal : float
        Optimal Tm in °C (default 52.0).
    tm_min : float
        Minimum acceptable Tm (default 40.0).
    tm_max : float
        Maximum acceptable Tm (default 65.0).

    Returns
    -------
    float
        Tm score in [0.0, 1.0].
    """
    if tm < tm_min or tm > tm_max:
        return 0.0

    if tm == tm_optimal:
        return 1.0

    # Linear interpolation: 1.0 at optimal, 0.5 at boundaries
    if tm < tm_optimal:
        # Below optimal: interpolate between tm_min (0.5) and tm_optimal (1.0)
        if tm_min == tm_optimal:
            return 1.0
        return 0.5 + 0.5 * (tm - tm_min) / (tm_optimal - tm_min)
    else:
        # Above optimal: interpolate between tm_optimal (1.0) and tm_max (0.5)
        if tm_max == tm_optimal:
            return 1.0
        return 1.0 - 0.5 * (tm - tm_optimal) / (tm_max - tm_optimal)


# ---------------------------------------------------------------------------
# Duplex ΔG computation
# ---------------------------------------------------------------------------


@dataclass
class DuplexResult:
    """Result of duplex ΔG computation.

    Attributes
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str
        Target subsequence (DNA alphabet).
    method : str
        Method used: "rnaduplex" or "nearest_neighbor".
    dg_duplex : float
        Duplex formation ΔG in kcal/mol (negative = favorable).
    start_target : int
        Start position of hybrid in target (1-based, 0 if NN method).
    end_target : int
        End position of hybrid in target (1-based, 0 if NN method).
    """

    aso_sequence: str
    target_sequence: str
    method: str
    dg_duplex: float
    start_target: int = 0
    end_target: int = 0


def calculate_duplex_dg_rnaduplex(
    aso_sequence: str,
    target_sequence: str,
) -> DuplexResult:
    """Compute duplex ΔG using ViennaRNA RNAduplex.

    RNAduplex computes the minimum free energy of inter-molecular
    hybridization, ignoring intramolecular structure. This is the
    primary method for ASO:target binding energy.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA or RNA, T → U conversion applied).
    target_sequence : str
        Target subsequence (DNA or RNA, T → U conversion applied).

    Returns
    -------
    DuplexResult
        Duplex ΔG and hybridization positions.

    Raises
    ------
    ThermodynamicsError
        If RNAduplex fails.
    """
    cache_path = _vienna_cache_path("rnaduplex", f"{aso_sequence}|{target_sequence}")
    cached = _load_vienna_cache(cache_path)
    if cached is not None:
        return DuplexResult(
            aso_sequence=aso_sequence.upper(),
            target_sequence=target_sequence.upper(),
            method="rnaduplex",
            dg_duplex=cached["mfe"],
            start_target=cached["start2"],
            end_target=cached["end2"],
        )

    from aso_designer.viennarna import run_rnaduplex

    try:
        result = run_rnaduplex(aso_sequence, target_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAduplex failed for ASO:target duplex: {e}") from e

    _save_vienna_cache(
        cache_path,
        {"mfe": result.mfe, "start2": result.start2, "end2": result.end2},
    )

    return DuplexResult(
        aso_sequence=aso_sequence.upper(),
        target_sequence=target_sequence.upper(),
        method="rnaduplex",
        dg_duplex=result.mfe,
        start_target=result.start2,
        end_target=result.end2,
    )


def calculate_duplex_dg_nn(
    aso_sequence: str,
    chemistry: str = "gapmer",
    temperature: float = DEFAULT_TEMPERATURE_C,
    modification_pattern: list[str] | None = None,
) -> DuplexResult:
    """Compute duplex ΔG using the nearest-neighbor model.

    The NN model computes ΔG for the ASO binding to its perfect
    complement. This is faster than RNAduplex but does not account
    for imperfect matches or bulges.

    The duplex modeled here is ASO:target-RNA, so chemistry routing
    uses the RNA/DNA hybrid parameters for DNA-like ASO chemistries
    (see ``_HYBRID_CHEMISTRIES``) and RNA parameters for RNA-like ones.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (determines hybrid vs RNA vs DNA NN parameters).
    temperature : float
        Temperature in °C (default 37 °C).
    modification_pattern : list of str or None
        Per-nucleotide modification labels; with the "lna" chemistry each
        LNA position stabilizes the duplex by ``_LNA_DG_PER_SITE``.

    Returns
    -------
    DuplexResult
        Duplex ΔG (NN-based).
    """
    use_rna = chemistry in _RNA_CHEMISTRIES
    use_hybrid = chemistry in _HYBRID_CHEMISTRIES
    dg = calculate_dg(aso_sequence, temperature=temperature, rna=use_rna, hybrid=use_hybrid)
    if chemistry == "lna" and modification_pattern:
        n_lna = sum(1 for mod in modification_pattern if mod == _LNA_MOD_LABEL)
        dg += _LNA_DG_PER_SITE * n_lna
    return DuplexResult(
        aso_sequence=aso_sequence.upper(),
        target_sequence="",
        method="nearest_neighbor",
        dg_duplex=dg,
    )


def calculate_duplex_dg(
    aso_sequence: str,
    target_sequence: str | None = None,
    chemistry: str = "gapmer",
    method: str = "nearest_neighbor",
    temperature: float = DEFAULT_TEMPERATURE_C,
    modification_pattern: list[str] | None = None,
) -> DuplexResult:
    """Compute duplex ΔG using the specified method.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str or None
        Target subsequence. Required for "rnaduplex" method.
    chemistry : str
        Chemistry name (for NN method parameter selection).
    method : str
        "nearest_neighbor" (fast, perfect match only) or
        "rnaduplex" (ViennaRNA, handles mismatches).
    temperature : float
        Temperature for NN calculation (default 37 °C).
    modification_pattern : list of str or None
        Per-nucleotide modification labels; with the "lna" chemistry each
        LNA position stabilizes the duplex.

    Returns
    -------
    DuplexResult
        Duplex ΔG result.
    """
    if method == "rnaduplex":
        if target_sequence is None:
            raise ThermodynamicsError("target_sequence is required for rnaduplex method")
        return calculate_duplex_dg_rnaduplex(aso_sequence, target_sequence)
    elif method == "nearest_neighbor":
        return calculate_duplex_dg_nn(
            aso_sequence,
            chemistry=chemistry,
            temperature=temperature,
            modification_pattern=modification_pattern,
        )
    else:
        raise ThermodynamicsError(
            f"Unknown duplex method: {method}. Use 'nearest_neighbor' or 'rnaduplex'."
        )


# ---------------------------------------------------------------------------
# Defective-duplex thermodynamics (single mismatch / dangling end)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DefectDuplexResult:
    """Thermodynamics of an ASO duplexed with an imperfect target site.

    Attributes
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet, 5'→3').
    target_site : str
        Target-site sequence exactly as evaluated (DNA alphabet).
    dh : float
        Duplex ΔH in kcal/mol.
    ds : float
        Duplex ΔS in cal/(mol·K).
    dg : float
        Duplex ΔG in kcal/mol at ``temperature``.
    defect : str
        One of "none", "internal_mismatch", "terminal_mismatch",
        "dangling_end".
    aso_position : int
        0-based ASO index of the defect (mismatched pair, or the paired
        position the dangling base stacks onto); −1 when defect is "none".
    target_position : int
        0-based index into ``target_site`` of the defect; −1 when "none".
    detail : str
        Human-readable description of the defect and parameters used.
    """

    aso_sequence: str
    target_site: str
    dh: float
    ds: float
    dg: float
    defect: str
    aso_position: int
    target_position: int
    detail: str


# Closing base pairs allowed by the defect tables: the four Watson-Crick
# pairs plus the two G·T wobble orientations (bottom base written second).
_DNA_CLOSURES: frozenset[tuple[str, str]] = frozenset(
    {("A", "T"), ("C", "G"), ("G", "C"), ("T", "A"), ("G", "T"), ("T", "G")}
)


def _combine_defect_tables(
    dg_table: dict[str, dict[str, float | None]],
    dh_table: dict[str, dict[str, float | None]],
) -> dict[str, dict[str, tuple[float, float]]]:
    """Combine published (ΔG°37, ΔH) defect tables into (ΔH, ΔS) tables.

    ΔS is derived from the published pair via the Gibbs relation
    ΔS = (ΔH − ΔG°37) × 1000 / 310.15, the same convention used for the
    DNA/RNA/hybrid NN tables above. Disallowed motifs ("." in the source
    tables, stored as None) are dropped; lookups of missing motifs raise
    :class:`ThermodynamicsError`.
    """
    combined: dict[str, dict[str, tuple[float, float]]] = {}
    for block, grid_dg in dg_table.items():
        grid_dh = dh_table[block]
        entries: dict[str, tuple[float, float]] = {}
        for motif, dg in grid_dg.items():
            dh = grid_dh[motif]
            if dg is None or dh is None:
                continue
            entries[motif] = (dh, (dh - dg) * 1000.0 / 310.15)
        combined[block] = entries
    return combined


_INT11_PARAMS = _combine_defect_tables(_INT11_DG37, _INT11_DH)
_TERMINAL_MISMATCH_PARAMS = _combine_defect_tables(_TERMINAL_MISMATCH_DG37, _TERMINAL_MISMATCH_DH)
_DANGLE_PARAMS = _combine_defect_tables(_DANGLE_DG37, _DANGLE_DH)


def _nn_duplex_params(top: str) -> tuple[float, float]:
    """(ΔH, ΔS) of a perfectly paired DNA duplex region (SantaLucia 1998).

    Sums the stacking terms over the top-strand sequence plus the
    per-terminal initiation corrections at the region's own two ends.
    No symmetry correction: this helper is used for baselines of
    defective duplexes, which are never self-complementary.
    """
    total_h = 0.0
    total_s = 0.0
    for i in range(len(top) - 1):
        h, s = _DNA_NN[top[i : i + 2]]
        total_h += h
        total_s += s
    for base in (top[0], top[-1]):
        h, s = _DNA_INIT_GC if base in ("G", "C") else _DNA_INIT_AT
        total_h += h
        total_s += s
    return total_h, total_s


def _lookup_defect(
    table: dict[str, dict[str, tuple[float, float]]],
    block: str,
    motif: str,
    description: str,
) -> tuple[float, float]:
    """Look up a (ΔH, ΔS) defect entry, raising on disallowed motifs."""
    try:
        return table[block][motif]
    except KeyError:
        raise ThermodynamicsError(
            f"Defect motif not tabulated ({description}): block {block}, motif {motif}"
        ) from None


def calculate_defect_duplex(
    aso_sequence: str,
    target_site: str,
    temperature: float = DEFAULT_TEMPERATURE_C,
) -> DefectDuplexResult:
    """Thermodynamics of an ASO duplexed with a (possibly imperfect) target.

    Models an ASO (top strand, DNA, 5'→3') bound to a target site using
    the SantaLucia 1998 DNA/DNA NN model extended with the NNDB defect
    tables (see the defect-parameter section above). The pairing
    convention matches :mod:`aso_designer.off_target`: the ASO binds the
    reverse complement of ``target_site``, i.e. ``aso[i]`` pairs with
    ``target_site[n - 1 - i]`` (antiparallel, complementary pairing).

    Exactly one defect is supported:

    - ``target_site`` same length as the ASO with zero mismatches: the
      perfect duplex (equivalent to :func:`calculate_dg`).
    - one internal mismatch: a 1×1 internal loop. Requires Watson-Crick
      or G·T-wobble pairs on both flanking sides (the tabulated
      closures). The tabulated 1×1 value replaces the two flanking
      stacks.
    - one terminal mismatch (at either end): the tabulated terminal-
      mismatch value replaces the terminal stack, added to the WC helix
      of the remaining pairs (NNDB composition; verified against the
      RNA worked examples, which the repo's RNA model reproduces
      exactly).
    - ``target_site`` longer by exactly one base: a single dangling
      nucleotide on the target strand. The aligned region must pair
      perfectly. The dangling value adds to the paired-region helix;
      5'-side overhangs use the bottom-5' table directly, 3'-side
      overhangs are evaluated in the rotated (strand-swapped) frame
      with the top-3' table, which is exact for this parameter set.

    Duplexes with more than one defect (two mismatches, mismatch
    adjacent to a dangle, mismatch next to a non-tabulated closing
    pair, overhangs on both strands) are outside the tabulated NN
    domain and raise :class:`ThermodynamicsError`.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet, 5'→3'). U is converted to T.
    target_site : str
        Target-site sequence, 5'→3', of length equal to the ASO
        (mismatch case) or longer by one base (dangling-end case).
        U is converted to T.
    temperature : float
        Temperature in °C for ΔG (default 37 °C).

    Returns
    -------
    DefectDuplexResult

    Raises
    ------
    ThermodynamicsError
        On invalid input, unsupported defect geometry, or motifs
        missing from the published tables.
    """
    aso = aso_sequence.upper().replace("U", "T")
    site = target_site.upper().replace("U", "T")
    for label, seq in (("ASO", aso), ("target site", site)):
        invalid = sorted(set(seq) - set("ACGT"))
        if invalid:
            raise ThermodynamicsError(f"Invalid base(s) {invalid} in {label} sequence: {seq}")

    n = len(aso)
    m = len(site)
    if n < 3:
        raise ThermodynamicsError(f"ASO too short for defect-duplex modeling: {n} nt (< 3)")
    if m not in (n, n + 1):
        raise ThermodynamicsError(
            f"Target site length {m} must equal the ASO length {n} "
            f"(mismatch) or exceed it by one base (dangling end): "
            f"{aso} vs {site}"
        )

    comp_map = str.maketrans("ACGT", "TGCA")
    revcomp_aso = aso.translate(comp_map)[::-1]

    def dg_at(dh: float, ds: float) -> float:
        return dh - (273.15 + temperature) * ds / 1000.0

    def finish(
        defect: str,
        aso_pos: int,
        target_pos: int,
        dh: float,
        ds: float,
        detail: str,
    ) -> DefectDuplexResult:
        return DefectDuplexResult(
            aso_sequence=aso,
            target_site=site,
            dh=dh,
            ds=ds,
            dg=dg_at(dh, ds),
            defect=defect,
            aso_position=aso_pos,
            target_position=target_pos,
            detail=detail,
        )

    if m == n:
        if site == revcomp_aso:
            dh, ds = calculate_nn_thermodynamics(aso)
            return finish(
                "none",
                -1,
                -1,
                dh,
                ds,
                "perfectly paired duplex (SantaLucia 1998 NN model)",
            )

        # under[i] = target-site base paired opposite aso[i] (bottom strand
        # read 3'→5' under the top strand). A position pairs (Watson-Crick
        # or G·T wobble) iff the pair is in the tabulated closure set.
        under = [site[n - 1 - i] for i in range(n)]
        mismatches = [i for i in range(n) if (aso[i], under[i]) not in _DNA_CLOSURES]
        if not mismatches:
            # Paired throughout via WC pairs and/or G·T wobbles (the site
            # need not be the strict reverse complement). Stacking and
            # initiation depend only on the top strand, so the energetics
            # equal the perfect duplex; the symmetry correction is skipped
            # because a duplex containing wobble pairs is never
            # self-complementary.
            dh, ds = _nn_duplex_params(aso)
            return finish(
                "none",
                -1,
                -1,
                dh,
                ds,
                "paired duplex with Watson-Crick and/or G·T wobble pairs "
                "(SantaLucia 1998 NN model)",
            )
        if len(mismatches) > 1:
            raise ThermodynamicsError(
                f"Unsupported duplex with {len(mismatches)} mismatches "
                f"(at most one is tabulated): {aso} vs {site}"
            )

        i = mismatches[0]
        top_base, bottom_base = aso[i], under[i]

        if 0 < i < n - 1:
            # 1x1 internal loop: requires tabulated (WC/wobble) flanking
            # closures; the tabulated value replaces the two flanking stacks.
            left = (aso[i - 1], under[i - 1])
            right = (aso[i + 1], under[i + 1])
            for side, closure in (("left", left), ("right", right)):
                if closure not in _DNA_CLOSURES:
                    raise ThermodynamicsError(
                        f"Unsupported internal mismatch at ASO position {i}: "
                        f"{side} flanking pair {closure[0]}·{closure[1]} is not "
                        f"a Watson-Crick or G·T-wobble pair"
                    )
            block = f"{left[0]}{left[1]}/{right[0]}{right[1]}"
            motif = f"{top_base}{bottom_base}"
            dh_defect, ds_defect = _lookup_defect(_INT11_PARAMS, block, motif, "1x1 internal loop")
            h1, s1 = _DNA_NN[aso[i - 1 : i + 1]]
            h2, s2 = _DNA_NN[aso[i : i + 2]]
            base_h, base_s = _nn_duplex_params(aso)
            dh = base_h - h1 - h2 + dh_defect
            ds = base_s - s1 - s2 + ds_defect
            return finish(
                "internal_mismatch",
                i,
                n - 1 - i,
                dh,
                ds,
                f"1x1 internal mismatch {top_base}·{bottom_base} at ASO "
                f"position {i} between closures {block} "
                f"(ΔG°37 {dh_defect - 310.15 * ds_defect / 1000.0:+.2f} "
                f"kcal/mol replacing two stacks)",
            )

        if i == n - 1:
            # 3'-terminal mismatch: WC helix through the penultimate pair +
            # the tabulated terminal-mismatch value (replaces the terminal
            # stack).
            closing = (aso[n - 2], under[n - 2])
            if closing not in _DNA_CLOSURES:
                raise ThermodynamicsError(
                    f"Unsupported terminal mismatch at ASO position {i}: "
                    f"penultimate pair {closing[0]}·{closing[1]} is not a "
                    f"Watson-Crick or G·T-wobble pair"
                )
            dh_base, ds_base = _nn_duplex_params(aso[: n - 1])
            dh_defect, ds_defect = _lookup_defect(
                _TERMINAL_MISMATCH_PARAMS,
                f"{closing[0]}{closing[1]}",
                f"{top_base}{bottom_base}",
                "terminal mismatch",
            )
            return finish(
                "terminal_mismatch",
                i,
                n - 1 - i,
                dh_base + dh_defect,
                ds_base + ds_defect,
                f"3'-terminal mismatch {top_base}·{bottom_base} adjacent to "
                f"closing pair {closing[0]}·{closing[1]}",
            )

        # i == 0: 5'-terminal mismatch, evaluated in the rotated frame
        # (strands swapped and reversed) where the mismatched base sits on
        # the top strand at the 3' end, matching the table orientation.
        closing = (under[1], aso[1])
        if closing not in _DNA_CLOSURES:
            raise ThermodynamicsError(
                f"Unsupported terminal mismatch at ASO position {i}: "
                f"penultimate pair {closing[0]}·{closing[1]} is not a "
                f"Watson-Crick or G·T-wobble pair"
            )
        dh_base, ds_base = _nn_duplex_params(aso[1:])
        dh_defect, ds_defect = _lookup_defect(
            _TERMINAL_MISMATCH_PARAMS,
            f"{closing[0]}{closing[1]}",
            f"{under[0]}{aso[0]}",
            "terminal mismatch",
        )
        return finish(
            "terminal_mismatch",
            0,
            n - 1,
            dh_base + dh_defect,
            ds_base + ds_defect,
            f"5'-terminal mismatch {top_base}·{bottom_base} adjacent to "
            f"closing pair {aso[1]}·{under[1]} (evaluated on the rotated "
            f"strand)",
        )

    # m == n + 1: single dangling nucleotide on the target strand. Exactly
    # one alignment must pair via tabulated (Watson-Crick/wobble) closures.
    five_prime = all((aso[i], site[n - i]) in _DNA_CLOSURES for i in range(n))
    three_prime = all((aso[i], site[n - 1 - i]) in _DNA_CLOSURES for i in range(n))
    if five_prime and three_prime:
        raise ThermodynamicsError(
            f"Ambiguous dangling-end alignment for {aso} vs {site}: the "
            f"unpaired base could sit at either end; disambiguate the "
            f"overhang side before calling"
        )
    if not five_prime and not three_prime:
        raise ThermodynamicsError(
            f"Dangling-end modeling requires the aligned region to pair "
            f"throughout (Watson-Crick or G·T wobble): {aso} vs {site}"
        )

    if five_prime:
        # The target's 5'-most base dangles beyond the ASO's 3' end:
        # bottom-5' dangle, closing pair (aso[n-1], site[1]).
        closing = f"{aso[n - 1]}{site[1]}"
        dh_defect, ds_defect = _lookup_defect(
            _DANGLE_PARAMS, f"{closing}|bottom_5p", site[0], "dangling end"
        )
        dh, ds = _nn_duplex_params(aso)
        return finish(
            "dangling_end",
            n - 1,
            0,
            dh + dh_defect,
            ds + ds_defect,
            f"5' dangling {site[0]} on the target strand adjacent to "
            f"closing pair {closing[0]}·{closing[1]}",
        )

    # The target's 3'-most base dangles beyond the ASO's 5' end: bottom-3'
    # dangle, evaluated in the rotated frame as a top-3' dangle with
    # closing pair (site[n-1], aso[0]).
    closing = f"{site[n - 1]}{aso[0]}"
    dh_defect, ds_defect = _lookup_defect(
        _DANGLE_PARAMS, f"{closing}|top_3p", site[n], "dangling end"
    )
    dh, ds = _nn_duplex_params(aso)
    return finish(
        "dangling_end",
        0,
        n,
        dh + dh_defect,
        ds + ds_defect,
        f"3' dangling {site[n]} on the target strand adjacent to closing "
        f"pair {closing[0]}·{closing[1]} (evaluated on the rotated strand)",
    )


# ---------------------------------------------------------------------------
# Self-structure ΔG computation
# ---------------------------------------------------------------------------


@dataclass
class SelfStructureResult:
    """Result of self-structure ΔG computation.

    Attributes
    ----------
    sequence : str
        ASO sequence.
    dg_self : float
        Self-folding ΔG in kcal/mol (more negative = more self-structure).
    structure : str
        Dot-bracket secondary structure (empty if NN method).
    method : str
        Method used: "rnafold" or "nearest_neighbor".
    penalty : float
        Self-structure penalty (0.0 if no significant self-structure).
    """

    sequence: str
    dg_self: float
    structure: str
    method: str
    penalty: float = 0.0


def calculate_self_structure_dg_rnafold(
    aso_sequence: str,
) -> SelfStructureResult:
    """Compute self-structure ΔG using ViennaRNA RNAfold.

    RNAfold computes the MFE secondary structure of the ASO alone.
    More negative ΔG indicates stronger self-folding, which reduces
    the effective concentration available for target binding.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA or RNA, T → U conversion applied).

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG and dot-bracket structure.

    Raises
    ------
    ThermodynamicsError
        If RNAfold fails.
    """
    cache_path = _vienna_cache_path("rnafold", aso_sequence)
    cached = _load_vienna_cache(cache_path)
    if cached is not None:
        return SelfStructureResult(
            sequence=aso_sequence.upper(),
            dg_self=cached["mfe"],
            structure=cached["structure"],
            method="rnafold",
        )

    from aso_designer.viennarna import run_rnafold

    try:
        result = run_rnafold(aso_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAfold failed for self-structure: {e}") from e

    _save_vienna_cache(cache_path, {"mfe": result.mfe, "structure": result.structure})

    return SelfStructureResult(
        sequence=aso_sequence.upper(),
        dg_self=result.mfe,
        structure=result.structure,
        method="rnafold",
    )


def calculate_self_structure_dg_nn(
    aso_sequence: str,
    chemistry: str = "gapmer",
    temperature: float = DEFAULT_TEMPERATURE_C,
    nn_fraction: float = _DEFAULT_NN_FRACTION,
) -> SelfStructureResult:
    """Compute self-structure ΔG using the nearest-neighbor model.

    This computes the ΔG of the ASO folding back on itself (hairpin).
    For short ASOs (15-25 nt), self-structure is typically minimal, so
    the NN model provides a quick approximation. For more accurate
    results, use the RNAfold method.

    The NN calculation gives the ΔG of the ASO paired with its perfect
    complement; since only part of the oligo can fold back on itself,
    the result is scaled by ``nn_fraction`` (default 0.3). This is a
    rough approximation; RNAfold is preferred.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (determines NN parameter set).
    temperature : float
        Temperature in °C (default 37 °C).
    nn_fraction : float
        Fraction of the perfect-duplex ΔG attributed to intramolecular
        self-structure (default 0.3, in (0, 1]).

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG (NN-based, no structure).
    """
    # Self-structure models the ASO duplexing with *itself* (ASO:ASO), so
    # same-strand parameters apply: DNA/DNA for DNA-like chemistries,
    # RNA/RNA for RNA-like ones — NOT the RNA/DNA hybrid table, which
    # describes ASO:target-RNA duplexes only.
    use_rna = chemistry in _RNA_CHEMISTRIES
    dg = calculate_dg(aso_sequence, temperature=temperature, rna=use_rna)
    dg_self = dg * nn_fraction
    return SelfStructureResult(
        sequence=aso_sequence.upper(),
        dg_self=dg_self,
        structure="",
        method="nearest_neighbor",
    )


def calculate_self_structure_penalty(
    dg_self: float,
    self_structure_max: float = -3.0,
) -> float:
    """Compute a penalty for self-structure.

    The penalty is 0.0 when ΔG ≥ 0 (no self-structure) and increases
    linearly as ΔG becomes more negative. At the threshold
    (self_structure_max), the penalty is 1.0.

    Parameters
    ----------
    dg_self : float
        Self-structure ΔG in kcal/mol.
    self_structure_max : float
        Threshold ΔG for maximum penalty (default -3.0 kcal/mol).

    Returns
    -------
    float
        Penalty in [0.0, 1.0+]. Values >1.0 indicate ΔG worse than threshold.
    """
    if dg_self >= 0:
        return 0.0
    if self_structure_max >= 0:
        return 0.0
    return abs(dg_self) / abs(self_structure_max)


def get_chemistry_self_structure_weight(chemistry: str) -> float:
    """Get the self-structure weight adjustment for a chemistry.

    Some chemistries (e.g., gapmers with DNA center) are less prone
    to self-structure issues, while others (e.g., 2'OMe) are more.

    Parameters
    ----------
    chemistry : str
        Chemistry name.

    Returns
    -------
    float
        Weight multiplier (1.0 = default, <1.0 = less sensitive,
        >1.0 = more sensitive).
    """
    try:
        presets = _load_chemistry_presets()
    except (OSError, yaml.YAMLError) as e:
        logger.warning("Could not load chemistry presets; using 1.0 weight: %s", e)
        return 1.0

    preset = presets.get(chemistry, {})
    return preset.get("self_structure_weight_adjust", 1.0)


def calculate_self_structure_dg(
    aso_sequence: str,
    chemistry: str = "gapmer",
    method: str = "nearest_neighbor",
    self_structure_max: float = -3.0,
    temperature: float = DEFAULT_TEMPERATURE_C,
    nn_fraction: float = _DEFAULT_NN_FRACTION,
) -> SelfStructureResult:
    """Compute self-structure ΔG and penalty with chemistry awareness.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (for weight adjustment and NN parameter selection).
    method : str
        "nearest_neighbor" (fast) or "rnafold" (ViennaRNA, accurate).
    self_structure_max : float
        Threshold for maximum penalty (default -3.0 kcal/mol).
    temperature : float
        Temperature for NN calculation (default 37 °C).
    nn_fraction : float
        Fraction of the perfect-duplex ΔG attributed to
        intramolecular self-structure in the NN approximation
        (default 0.3, in (0, 1]). Ignored for the "rnafold" method.

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG, structure, and chemistry-adjusted penalty.
    """
    if method == "rnafold":
        result = calculate_self_structure_dg_rnafold(aso_sequence)
    elif method == "nearest_neighbor":
        result = calculate_self_structure_dg_nn(
            aso_sequence,
            chemistry=chemistry,
            temperature=temperature,
            nn_fraction=nn_fraction,
        )
    else:
        raise ThermodynamicsError(
            f"Unknown self-structure method: {method}. Use 'nearest_neighbor' or 'rnafold'."
        )

    # Compute base penalty
    base_penalty = calculate_self_structure_penalty(result.dg_self, self_structure_max)

    # Apply chemistry-specific weight adjustment
    weight = get_chemistry_self_structure_weight(chemistry)
    result.penalty = base_penalty * weight

    return result


# ---------------------------------------------------------------------------
# Target structure ΔG and net binding energy
# ---------------------------------------------------------------------------


@dataclass
class TargetStructureResult:
    """Result of target structure ΔG computation at a binding site.

    Attributes
    ----------
    target_sequence : str
        Target subsequence at the binding site.
    dg_target : float
        Target folding ΔG in kcal/mol (more negative = more structured).
    structure : str
        Dot-bracket secondary structure (empty if NN method).
    method : str
        Method used: "rnafold" or "nearest_neighbor".
    """

    target_sequence: str
    dg_target: float
    structure: str
    method: str


def calculate_target_structure_dg_rnafold(
    target_sequence: str,
) -> TargetStructureResult:
    """Compute target structure ΔG at a binding site using RNAfold.

    The target subsequence at the ASO binding site may have local
    secondary structure that must be opened for ASO binding. RNAfold
    computes the MFE of this subsequence.

    Parameters
    ----------
    target_sequence : str
        Target subsequence at the binding site (DNA or RNA).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG and dot-bracket structure.
    """
    cache_path = _vienna_cache_path("rnafold", target_sequence)
    cached = _load_vienna_cache(cache_path)
    if cached is not None:
        return TargetStructureResult(
            target_sequence=target_sequence.upper(),
            dg_target=cached["mfe"],
            structure=cached["structure"],
            method="rnafold",
        )

    from aso_designer.viennarna import run_rnafold

    try:
        result = run_rnafold(target_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAfold failed for target structure: {e}") from e

    _save_vienna_cache(cache_path, {"mfe": result.mfe, "structure": result.structure})

    return TargetStructureResult(
        target_sequence=target_sequence.upper(),
        dg_target=result.mfe,
        structure=result.structure,
        method="rnafold",
    )


def calculate_target_structure_dg_nn(
    target_sequence: str,
    temperature: float = DEFAULT_TEMPERATURE_C,
) -> TargetStructureResult:
    """Compute target structure ΔG using the NN model.

    Uses RNA NN parameters since the target is mRNA.

    Parameters
    ----------
    target_sequence : str
        Target subsequence (DNA alphabet, converted to RNA internally).
    temperature : float
        Temperature in °C (default 37 °C).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG (NN-based).
    """
    dg = calculate_dg(target_sequence, temperature=temperature, rna=True)
    return TargetStructureResult(
        target_sequence=target_sequence.upper(),
        dg_target=dg,
        structure="",
        method="nearest_neighbor",
    )


def calculate_target_structure_dg(
    target_sequence: str,
    method: str = "nearest_neighbor",
    temperature: float = DEFAULT_TEMPERATURE_C,
) -> TargetStructureResult:
    """Compute target structure ΔG using the specified method.

    Parameters
    ----------
    target_sequence : str
        Target subsequence at the binding site.
    method : str
        "nearest_neighbor" (fast) or "rnafold" (ViennaRNA, accurate).
    temperature : float
        Temperature for NN calculation (default 37 °C).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG result.
    """
    if method == "rnafold":
        return calculate_target_structure_dg_rnafold(target_sequence)
    elif method == "nearest_neighbor":
        return calculate_target_structure_dg_nn(target_sequence, temperature=temperature)
    else:
        raise ThermodynamicsError(
            f"Unknown target structure method: {method}. Use 'nearest_neighbor' or 'rnafold'."
        )


# ---------------------------------------------------------------------------
# Net binding energy
# ---------------------------------------------------------------------------


@dataclass
class ThermodynamicSummary:
    """Complete thermodynamic profile of an ASO candidate.

    Combines Tm, duplex ΔG, self-structure ΔG, target structure ΔG,
    and net binding energy into a single summary object.

    Attributes
    ----------
    tm : float
        Melting temperature in °C.
    tm_score : float
        Tm score (0–1, 1 = optimal).
    dg_duplex : float
        Duplex formation ΔG in kcal/mol.
    dg_self : float
        Self-structure ΔG in kcal/mol.
    self_structure_penalty : float
        Chemistry-adjusted self-structure penalty.
    dg_target : float
        Target structure ΔG at binding site in kcal/mol.
    dg_net : float
        Net binding ΔG = ΔG_duplex - ΔG_self - ΔG_target.
    binding_class : str
        Binding strength classification: "strong", "moderate", or "weak".
    """

    tm: float
    tm_score: float
    dg_duplex: float
    dg_self: float
    self_structure_penalty: float
    dg_target: float
    dg_net: float
    binding_class: str


def calculate_net_binding_dg(
    dg_duplex: float,
    dg_self: float,
    dg_target: float,
) -> float:
    """Compute net binding ΔG from component energies.

    The thermodynamic cycle for ASO binding:
    1. ASO unfolds (cost: -ΔG_self, since ΔG_self is negative)
    2. Target opens at binding site (cost: -ΔG_target)
    3. ASO:target duplex forms (gain: ΔG_duplex)

    ΔG_net = ΔG_duplex - ΔG_self - ΔG_target

    More negative ΔG_net = more favorable binding.

    Parameters
    ----------
    dg_duplex : float
        Duplex formation ΔG (negative = favorable).
    dg_self : float
        Self-structure ΔG (negative = more self-structure).
    dg_target : float
        Target structure ΔG (negative = more structured target).

    Returns
    -------
    float
        Net binding ΔG in kcal/mol.
    """
    return dg_duplex - dg_self - dg_target


def classify_binding(
    dg_net: float,
    strong_threshold: float = -10.0,
    moderate_threshold: float = -5.0,
) -> str:
    """Classify binding strength by net binding ΔG.

    Parameters
    ----------
    dg_net : float
        Net binding ΔG in kcal/mol.
    strong_threshold : float
        ΔG at or below this is "strong" (default -10.0 kcal/mol).
    moderate_threshold : float
        ΔG at or below this (but above strong_threshold) is "moderate"
        (default -5.0 kcal/mol). Above this is "weak".

    Notes
    -----
    The default thresholds are heuristic defaults, not constants from a
    single publication. They are calibrated so that typical 18-22 nt
    gapmer designs — whose perfect-duplex ΔG37 values fall roughly in
    the -20 to -30 kcal/mol range under the SantaLucia (1998) DNA and
    Xia et al. (1998) RNA parameter sets — classify as "strong" after
    subtracting self-structure and target-opening costs. Tune both
    thresholds per chemistry and assay context; they are exposed as
    parameters precisely because no universal published cutoff exists.

    Returns
    -------
    str
        One of "strong", "moderate", or "weak".
    """
    if dg_net <= strong_threshold:
        return "strong"
    elif dg_net <= moderate_threshold:
        return "moderate"
    else:
        return "weak"


def compute_thermodynamic_summary(
    aso_sequence: str,
    target_sequence: str | None = None,
    chemistry: str = "gapmer",
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    tm_optimal: float = 52.0,
    tm_min: float = 40.0,
    tm_max: float = 65.0,
    self_structure_max: float = -3.0,
    method: str = "nearest_neighbor",
    temperature: float = DEFAULT_TEMPERATURE_C,
    nn_fraction: float = _DEFAULT_NN_FRACTION,
    modification_pattern: list[str] | None = None,
) -> ThermodynamicSummary:
    """Compute a complete thermodynamic summary for an ASO candidate.

    Combines Tm, duplex ΔG, self-structure ΔG, target structure ΔG,
    and net binding energy into a single summary.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str or None
        Target subsequence at binding site. Required for "rnaduplex"
        duplex method and "rnafold" target method.
    chemistry : str
        Chemistry name (default "gapmer").
    na_conc : float
        Sodium concentration in M (default 0.1).
    oligo_conc : float
        Oligo concentration in M (default 1 µM).
    mg_conc : float
        Magnesium concentration in M (0 disables Mg2+ correction).
    dntp_conc : float
        dNTP concentration in M (chelates Mg2+).
    tm_optimal : float
        Optimal Tm in °C (default 52.0).
    tm_min : float
        Minimum acceptable Tm (default 40.0).
    tm_max : float
        Maximum acceptable Tm (default 65.0).
    self_structure_max : float
        Self-structure threshold (default -3.0 kcal/mol).
    method : str
        Computation method: "nearest_neighbor" (fast) or
        "rnafold"/"rnaduplex" (ViennaRNA, accurate).
    temperature : float
        Temperature for ΔG calculations (default 37 °C).
    nn_fraction : float
        Fraction of the perfect-duplex ΔG attributed to intramolecular
        self-structure in the NN approximation (default 0.3, in (0, 1]).

    Returns
    -------
    ThermodynamicSummary
        Complete thermodynamic profile.
    """
    # Tm
    tm = calculate_tm_with_chemistry(
        aso_sequence,
        chemistry=chemistry,
        na_conc=na_conc,
        oligo_conc=oligo_conc,
        mg_conc=mg_conc,
        dntp_conc=dntp_conc,
        modification_pattern=modification_pattern,
    )
    tm_score = score_tm(tm, tm_optimal=tm_optimal, tm_min=tm_min, tm_max=tm_max)

    # Duplex ΔG
    if method == "rnaduplex" and target_sequence is not None:
        duplex = calculate_duplex_dg(
            aso_sequence, target_sequence, chemistry=chemistry, method="rnaduplex"
        )
    else:
        duplex = calculate_duplex_dg(
            aso_sequence,
            chemistry=chemistry,
            method="nearest_neighbor",
            temperature=temperature,
            modification_pattern=modification_pattern,
        )

    # Self-structure ΔG
    self_struct = calculate_self_structure_dg(
        aso_sequence,
        chemistry=chemistry,
        method=method if method != "rnaduplex" else "rnafold",
        self_structure_max=self_structure_max,
        temperature=temperature,
        nn_fraction=nn_fraction,
    )

    # Target structure ΔG
    if target_sequence is not None:
        target_struct = calculate_target_structure_dg(
            target_sequence,
            method=method if method != "rnaduplex" else "rnafold",
            temperature=temperature,
        )
        dg_target = target_struct.dg_target
    else:
        # Without target sequence, assume no target structure cost
        dg_target = 0.0

    # Net binding energy
    dg_net = calculate_net_binding_dg(
        dg_duplex=duplex.dg_duplex,
        dg_self=self_struct.dg_self,
        dg_target=dg_target,
    )

    # Classify binding
    binding_class = classify_binding(dg_net)

    return ThermodynamicSummary(
        tm=tm,
        tm_score=tm_score,
        dg_duplex=duplex.dg_duplex,
        dg_self=self_struct.dg_self,
        self_structure_penalty=self_struct.penalty,
        dg_target=dg_target,
        dg_net=dg_net,
        binding_class=binding_class,
    )
