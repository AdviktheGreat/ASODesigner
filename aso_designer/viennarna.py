"""ViennaRNA integration layer for ASODesigner.

Wraps ViennaRNA command-line tools (RNAfold, RNAplfold, RNAcofold,
RNAduplex) with robust subprocess management, output parsing, and
installation verification. Inspired by GuideX's RNAfold verification
but redesigned with multi-tool support, conda/homebrew/source install
guidance, and multi-strategy output parsing.
"""

import re
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from aso_designer.utils.exceptions import ViennaRNAError
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import dna_to_rna

logger = get_logger(__name__)

# Minimum required ViennaRNA version (major.minor.patch)
_MIN_VERSION = (2, 4, 0)

# Install instructions for different platforms
_INSTALL_INSTRUCTIONS = """
ViennaRNA installation instructions:

  macOS (Homebrew):
    brew install vienna-rna

  Conda:
    conda install -c bioconda viennarna

  Ubuntu/Debian:
    sudo apt-get install viennarna

  From source:
    https://www.tbi.univie.ac.at/RNA/#source
"""


@dataclass
class ViennaRNAVersion:
    """Parsed ViennaRNA version information."""

    major: int
    minor: int
    patch: int

    @property
    def version_string(self) -> str:
        """Formatted version string (e.g. '2.6.4')."""
        return f"{self.major}.{self.minor}.{self.patch}"

    def __ge__(self, other: tuple[int, int, int]) -> bool:
        return (self.major, self.minor, self.patch) >= other


def verify_viennarna() -> ViennaRNAVersion:
    """Verify ViennaRNA is installed and meets the minimum version.

    Checks for RNAfold in PATH and parses its version string.

    Returns
    -------
    ViennaRNAVersion
        Parsed version information.

    Raises
    ------
    ViennaRNAError
        If RNAfold is not found or the version is too old.
    """
    rnafold_path = shutil.which("RNAfold")
    if rnafold_path is None:
        raise ViennaRNAError("RNAfold not found in PATH.\n" + _INSTALL_INSTRUCTIONS)

    try:
        result = subprocess.run(
            [rnafold_path, "--version"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        )
    except subprocess.SubprocessError as e:
        raise ViennaRNAError(f"Failed to execute RNAfold --version: {e}") from e

    version = _parse_version(result.stdout)
    if version is None:
        raise ViennaRNAError(
            f"Could not parse ViennaRNA version from output: {result.stdout[:100]}"
        )

    if not version >= _MIN_VERSION:
        raise ViennaRNAError(
            f"ViennaRNA {_MIN_VERSION[0]}.{_MIN_VERSION[1]}.{_MIN_VERSION[2]}+ required, "
            f"found {version.version_string}.\n" + _INSTALL_INSTRUCTIONS
        )

    logger.info("ViennaRNA %s verified at %s", version.version_string, rnafold_path)
    return version


def _parse_version(output: str) -> ViennaRNAVersion | None:
    """Parse a ViennaRNA version string from command output.

    Handles formats like:
    - "RNAfold 2.6.4"
    - "ViennaRNA 2.5.1"
    - "2.4.6"
    """
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", output)
    if match:
        return ViennaRNAVersion(
            major=int(match.group(1)),
            minor=int(match.group(2)),
            patch=int(match.group(3)),
        )
    return None


def _get_tool_path(tool_name: str) -> str:
    """Get the path to a ViennaRNA tool, raising if not found.

    Parameters
    ----------
    tool_name : str
        Tool name: "RNAfold", "RNAplfold", "RNAcofold", "RNAduplex".

    Returns
    -------
    str
        Full path to the tool executable.

    Raises
    ------
    ViennaRNAError
        If the tool is not found in PATH.
    """
    path = shutil.which(tool_name)
    if path is None:
        raise ViennaRNAError(f"{tool_name} not found in PATH.\n" + _INSTALL_INSTRUCTIONS)
    return path


# ---------------------------------------------------------------------------
# RNAfold wrapper
# ---------------------------------------------------------------------------


@dataclass
class RNAfoldResult:
    """Result of an RNAfold computation.

    Attributes
    ----------
    sequence : str
        RNA sequence that was folded.
    structure : str
        Dot-bracket secondary structure notation.
    mfe : float
        Minimum free energy in kcal/mol.
    """

    sequence: str
    structure: str
    mfe: float


def run_rnafold(
    sequence: str,
    no_ps: bool = True,
    partition_function: bool = False,
) -> RNAfoldResult:
    """Run RNAfold on a sequence and parse the MFE structure.

    Converts DNA input (T → U) and parses the output using a
    multi-strategy approach: line-based parsing first, regex fallback
    second. Inspired by GuideX's _calculate_mfe but extended to
    return the full structure, not just MFE.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive).
    no_ps : bool
        Suppress PostScript output (default True).
    partition_function : bool
        Compute partition function in addition to MFE (default False).

    Returns
    -------
    RNAfoldResult
        Sequence, dot-bracket structure, and MFE.

    Raises
    ------
    ViennaRNAError
        If RNAfold fails or output cannot be parsed.
    """
    rna_seq = dna_to_rna(sequence)
    tool_path = _get_tool_path("RNAfold")

    args = [tool_path]
    if no_ps:
        args.append("--noPS")
    if partition_function:
        args.append("-p")

    input_data = f">{sequence}\n{rna_seq}\n"

    try:
        result = subprocess.run(
            args,
            input=input_data,
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except subprocess.SubprocessError as e:
        raise ViennaRNAError(f"RNAfold execution failed: {e}") from e

    return _parse_rnafold_output(rna_seq, result.stdout)


def _parse_rnafold_output(rna_seq: str, stdout: str) -> RNAfoldResult:
    """Parse RNAfold stdout with multi-strategy approach.

    Strategy 1: Find the line starting with the RNA sequence,
    extract the dot-bracket and MFE from the following line.
    Strategy 2: Regex fallback for MFE only.
    """
    lines = stdout.strip().split("\n")

    # Skip header lines (lines starting with '>')
    seq_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith(">"):
            continue
        if line.strip() and not line.startswith(">"):
            seq_line_idx = i
            break

    if seq_line_idx is not None and seq_line_idx + 1 < len(lines):
        structure_line = lines[seq_line_idx + 1].strip()
        # Structure line: "((((...)))) (-2.30)"
        # Split on last space to separate structure from energy
        parts = structure_line.rsplit(" ", 1)
        if len(parts) == 2:
            structure = parts[0]
            mfe_str = parts[1].strip("()")
            try:
                mfe = float(mfe_str)
                return RNAfoldResult(sequence=rna_seq, structure=structure, mfe=mfe)
            except ValueError:
                pass

    # Strategy 2: Regex fallback for MFE
    mfe_match = re.search(r"\(\s*(-?\d+\.\d+)\s*\)", stdout)
    if mfe_match:
        mfe = float(mfe_match.group(1))
        # Try to extract structure from the line before the MFE
        structure = ""
        for line in lines:
            if "(" in line and ")" in line and not line.startswith(">"):
                bracket_part = line.rsplit(" ", 1)[0]
                if all(c in "()." for c in bracket_part):
                    structure = bracket_part
                    break
        return RNAfoldResult(sequence=rna_seq, structure=structure, mfe=mfe)

    raise ViennaRNAError(f"Could not parse RNAfold output: {stdout[:200]}")


# ---------------------------------------------------------------------------
# RNAplfold wrapper (local accessibility)
# ---------------------------------------------------------------------------


@dataclass
class RNAplfoldResult:
    """Result of an RNAplfold computation.

    Attributes
    ----------
    positions : list[int]
        1-based positions in the sequence.
    unpaired_probs : list[float]
        Probability of being unpaired at each position (0.0–1.0) for the
        shortest requested unpaired-region length.
    unpaired_lengths : list[int]
        The requested unpaired-region lengths, in output-column order.
    unpaired_probs_by_length : dict[int, list[float]]
        Per-length unpaired probabilities: for each requested length L,
        the probability that a stretch of L nucleotides starting at each
        position is unpaired.
    """

    positions: list[int]
    unpaired_probs: list[float]
    unpaired_lengths: list[int] = field(default_factory=lambda: [1])
    unpaired_probs_by_length: dict[int, list[float]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.unpaired_probs_by_length:
            length = self.unpaired_lengths[0] if self.unpaired_lengths else 1
            self.unpaired_probs_by_length = {length: list(self.unpaired_probs)}


def run_rnaplfold(
    sequence: str,
    window_size: int = 150,
    max_span: int = 100,
    unpaired_region: int | Sequence[int] = 1,
) -> RNAplfoldResult:
    """Run RNAplfold to compute local unpaired probabilities.

    RNAplfold computes the probability that each nucleotide is
    unpaired within a sliding window. This is the key tool for
    ASO target accessibility analysis. Multiple unpaired-region
    lengths can be requested at once, producing per-length
    probabilities suitable for weighting multi-nucleotide unpaired
    stretches.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive, T → U conversion applied).
    window_size : int
        Sliding window size in nucleotides (default 150).
    max_span : int
        Maximum base pair distance (default 100).
    unpaired_region : int or sequence of int
        Unpaired-region length(s) to consider (default 1 =
        per-nucleotide accessibility). A sequence such as ``(1, 2, 3,
        4)`` is passed to RNAplfold as ``-u 1,2,3,4`` and yields one
        probability column per length.

    Returns
    -------
    RNAplfoldResult
        Per-position unpaired probabilities, per requested length.

    Raises
    ------
    ViennaRNAError
        If RNAplfold fails or output cannot be parsed.
    """
    import tempfile

    lengths = (
        [unpaired_region] if isinstance(unpaired_region, int) else sorted(set(unpaired_region))
    )
    if not lengths or any(x < 1 for x in lengths):
        raise ViennaRNAError(f"Invalid unpaired_region lengths: {lengths} (must be >= 1)")

    rna_seq = dna_to_rna(sequence)
    tool_path = _get_tool_path("RNAplfold")

    with tempfile.TemporaryDirectory() as tmpdir:
        args = [
            tool_path,
            "-W",
            str(window_size),
            "-L",
            str(max_span),
            "-u",
            ",".join(str(x) for x in lengths),
            "--id",
            "asodesigner",
            "-o",
            tmpdir + "/",
        ]

        input_data = f">asodesigner\n{rna_seq}\n"

        try:
            subprocess.run(
                args,
                input=input_data,
                capture_output=True,
                text=True,
                check=True,
                timeout=300,
                cwd=tmpdir,
            )
        except subprocess.SubprocessError as e:
            raise ViennaRNAError(f"RNAplfold execution failed: {e}") from e

        # Find the lunp (unpaired probability) output file
        lunp_files = list(Path(tmpdir).glob("*_lunp"))
        if not lunp_files:
            raise ViennaRNAError("RNAplfold did not produce an unpaired probability file")

        return _parse_rnaplfold_lunp(lunp_files[0], tuple(lengths))


def _parse_rnaplfold_lunp(
    filepath: Path, unpaired_lengths: tuple[int, ...] = (1,)
) -> RNAplfoldResult:
    """Parse an RNAplfold lunp (unpaired probability) output file.

    For a single requested length the file format is::

        #unpaired probabilities
        #i   p_u   u
        1    0.95  1
        2    0.87  1

    With several requested lengths (``-u 1,2,3,4``) the rows carry one
    probability column per requested length after the position::

        1    0.95  0.80  0.70  0.60

    The parser accepts both layouts: a row with exactly one float per
    requested length is taken as-is; a row with one extra trailing
    near-integer column equal to a requested length (the legacy ``u``
    column) is trimmed before assignment.
    """
    lengths = list(unpaired_lengths) or [1]
    k = len(lengths)
    positions: list[int] = []
    columns: list[list[float]] = [[] for _ in range(k)]

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                pos = int(parts[0])
            except ValueError:
                continue
            floats: list[float] = []
            for token in parts[1:]:
                try:
                    floats.append(float(token))
                except ValueError:
                    break
            if len(floats) == k:
                probs = floats
            elif len(floats) == k + 1 and any(
                abs(floats[-1] - length) < 1e-9 for length in lengths
            ):
                probs = floats[:-1]
            else:
                raise ViennaRNAError(
                    f"Could not parse lunp row with {len(floats)} value(s) for "
                    f"{k} requested length(s): {line[:100]}"
                )
            positions.append(pos)
            for idx, value in enumerate(probs):
                columns[idx].append(value)

    if not positions:
        raise ViennaRNAError(f"No unpaired probabilities found in {filepath}")

    probs_by_length = dict(zip(lengths, columns, strict=True))
    logger.debug("Parsed %d unpaired probabilities from %s", len(positions), filepath)
    return RNAplfoldResult(
        positions=positions,
        unpaired_probs=probs_by_length[min(lengths)],
        unpaired_lengths=lengths,
        unpaired_probs_by_length=probs_by_length,
    )


# ---------------------------------------------------------------------------
# RNAup wrapper (accessibility by RNA opening energy)
# ---------------------------------------------------------------------------


@dataclass
class RNAupResult:
    """Result of an RNAup computation.

    Attributes
    ----------
    sequence : str
        RNA sequence analyzed.
    positions : list[int]
        1-based positions in the sequence.
    unpaired_probs : list[float]
        Probability p(0) that each individual nucleotide is unpaired.
    region_probs : list[float]
        Probability that the region of ``region_length`` nucleotides
        starting at each position is unbound.
    opening_dg : list[float]
        Per-position opening energy ΔG_open in kcal/mol (the energy to
        open the region; larger = less accessible).
    region_length : int
        Length of the unpaired region the probabilities refer to.
    """

    sequence: str
    positions: list[int]
    unpaired_probs: list[float]
    region_probs: list[float]
    opening_dg: list[float]
    region_length: int


def _parse_rnaup_stdout(stdout: str, sequence: str) -> RNAupResult:
    """Parse RNAup stdout into an :class:`RNAupResult`.

    RNAup prints a header (``RNAup output in file: ...``), the input
    FASTA, then a region block::

        >t_1
        AUGCUAGC...
        [region 4]
           1    0.98    0.90    -0.1053
           2    ...

    Each data row holds the 1-based position ``i``, the per-nucleotide
    unpaired probability ``p(0)``, the probability that the region of
    length L starting at ``i`` is unbound, and the opening energy
    ``dG(i)`` in kcal/mol.
    """
    region_length: int | None = None
    positions: list[int] = []
    unpaired_probs: list[float] = []
    region_probs: list[float] = []
    opening_dg: list[float] = []

    in_block = False
    for line in stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        match = re.search(r"\[region\s+(\d+)\]", stripped)
        if match:
            region_length = int(match.group(1))
            in_block = True
            continue
        if not in_block:
            continue
        parts = stripped.split()
        try:
            pos = int(parts[0])
        except ValueError:
            continue  # header/comment lines inside the block
        floats: list[float] = []
        for token in parts[1:]:
            try:
                floats.append(float(token))
            except ValueError:
                break
        if len(floats) < 3:
            raise ViennaRNAError(f"Could not parse RNAup output row: {stripped[:100]}")
        positions.append(pos)
        unpaired_probs.append(floats[0])
        region_probs.append(floats[1])
        opening_dg.append(floats[2])

    if region_length is None or not positions:
        raise ViennaRNAError(f"Could not parse RNAup output: {stdout[:200]}")

    logger.debug("Parsed %d RNAup region-%d rows", len(positions), region_length)
    return RNAupResult(
        sequence=sequence,
        positions=positions,
        unpaired_probs=unpaired_probs,
        region_probs=region_probs,
        opening_dg=opening_dg,
        region_length=region_length,
    )


def run_rnaup(
    sequence: str,
    unpaired_region: int = 4,
) -> RNAupResult:
    """Run RNAup to compute per-site RNA opening energies.

    RNAup computes, for each nucleotide, the probability that the
    nucleotide and the following ``unpaired_region``-long stretch remain
    unbound, and converts that into an opening energy
    ΔG_open = −kT·ln(p_unbound). This complements the RNAplfold
    unpaired-probability screen with an explicit per-site opening
    energy.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive, T → U conversion applied).
    unpaired_region : int
        Maximum length of the unpaired region to consider (default 4,
        RNAup's own default).

    Returns
    -------
    RNAupResult
        Per-position unpaired probabilities, region probabilities, and
        opening energies.

    Raises
    ------
    ViennaRNAError
        If RNAup fails or its output cannot be parsed.
    """
    import tempfile

    rna_seq = dna_to_rna(sequence)
    tool_path = _get_tool_path("RNAup")

    with tempfile.TemporaryDirectory() as tmpdir:
        args = [
            tool_path,
            "--id",
            "asodesigner",
            "-u",
            str(unpaired_region),
        ]

        input_data = f">asodesigner\n{rna_seq}\n"

        try:
            proc = subprocess.run(
                args,
                input=input_data,
                capture_output=True,
                text=True,
                check=True,
                timeout=300,
                cwd=tmpdir,
            )
        except subprocess.SubprocessError as e:
            raise ViennaRNAError(f"RNAup execution failed: {e}") from e

        return _parse_rnaup_stdout(proc.stdout, rna_seq)


# ---------------------------------------------------------------------------
# RNAcofold wrapper (two-sequence hybridization with intramolecular structure)
# ---------------------------------------------------------------------------


@dataclass
class RNAcofoldResult:
    """Result of an RNAcofold computation.

    Attributes
    ----------
    sequence1 : str
        First RNA sequence.
    sequence2 : str
        Second RNA sequence.
    structure : str
        Dot-bracket structure with '&' separator.
    mfe : float
        Minimum free energy of the hybrid in kcal/mol.
    """

    sequence1: str
    sequence2: str
    structure: str
    mfe: float


def run_rnacofold(
    sequence1: str,
    sequence2: str,
    no_ps: bool = True,
) -> RNAcofoldResult:
    """Run RNAcofold on two sequences and parse the hybrid MFE.

    RNAcofold computes the minimum free energy of two sequences
    that can form both intra- and inter-molecular base pairs.
    This is useful for computing the binding energy of an ASO
    to its target when both can form self-structure.

    Parameters
    ----------
    sequence1 : str
        First sequence (DNA or RNA, T → U conversion applied).
    sequence2 : str
        Second sequence (DNA or RNA, T → U conversion applied).
    no_ps : bool
        Suppress PostScript output (default True).

    Returns
    -------
    RNAcofoldResult
        Sequences, dot-bracket structure, and MFE.

    Raises
    ------
    ViennaRNAError
        If RNAcofold fails or output cannot be parsed.
    """
    rna1 = dna_to_rna(sequence1)
    rna2 = dna_to_rna(sequence2)
    tool_path = _get_tool_path("RNAcofold")

    args = [tool_path]
    if no_ps:
        args.append("--noPS")

    input_data = f">seq1\n{rna1}&{rna2}\n"

    try:
        result = subprocess.run(
            args,
            input=input_data,
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except subprocess.SubprocessError as e:
        raise ViennaRNAError(f"RNAcofold execution failed: {e}") from e

    return _parse_rnacofold_output(rna1, rna2, result.stdout)


def _parse_rnacofold_output(rna1: str, rna2: str, stdout: str) -> RNAcofoldResult:
    """Parse RNAcofold stdout.

    Output format:
    >seq1
    AUGCAUGC&AUGCAUGC
    ((((...))))&((((...))))  (-4.60)
    """
    lines = stdout.strip().split("\n")

    for i, line in enumerate(lines):
        if line.startswith(">"):
            continue
        if "&" in line:
            # Next line has the structure and MFE
            if i + 1 < len(lines):
                structure_line = lines[i + 1].strip()
                parts = structure_line.rsplit(" ", 1)
                if len(parts) == 2:
                    structure = parts[0]
                    mfe_str = parts[1].strip("()")
                    try:
                        mfe = float(mfe_str)
                        return RNAcofoldResult(
                            sequence1=rna1,
                            sequence2=rna2,
                            structure=structure,
                            mfe=mfe,
                        )
                    except ValueError:
                        pass

    # Regex fallback
    mfe_match = re.search(r"\(\s*(-?\d+\.\d+)\s*\)", stdout)
    if mfe_match:
        mfe = float(mfe_match.group(1))
        structure = ""
        for line in lines:
            if "&" in line and not line.startswith(">"):
                bracket_part = line.rsplit(" ", 1)[0]
                if all(c in "().&." for c in bracket_part):
                    structure = bracket_part
                    break
        return RNAcofoldResult(sequence1=rna1, sequence2=rna2, structure=structure, mfe=mfe)

    raise ViennaRNAError(f"Could not parse RNAcofold output: {stdout[:200]}")


# ---------------------------------------------------------------------------
# RNAduplex wrapper (inter-molecular hybridization only)
# ---------------------------------------------------------------------------


@dataclass
class RNAduplexResult:
    """Result of an RNAduplex computation.

    Attributes
    ----------
    sequence1 : str
        First RNA sequence.
    sequence2 : str
        Second RNA sequence.
    structure : str
        Dot-bracket structure with '&' separator.
    mfe : float
        Minimum free energy of the duplex in kcal/mol.
    start1 : int
        Start position of hybrid in sequence 1 (1-based).
    end1 : int
        End position of hybrid in sequence 1 (1-based).
    start2 : int
        Start position of hybrid in sequence 2 (1-based).
    end2 : int
        End position of hybrid in sequence 2 (1-based).
    """

    sequence1: str
    sequence2: str
    structure: str
    mfe: float
    start1: int
    end1: int
    start2: int
    end2: int


def run_rnaduplex(
    sequence1: str,
    sequence2: str,
    no_ps: bool = True,
) -> RNAduplexResult:
    """Run RNAduplex on two sequences and parse the duplex MFE.

    RNAduplex computes the minimum free energy of the optimal
    inter-molecular hybridization between two sequences, ignoring
    intramolecular structure. This is the primary tool for ASO:target
    binding energy computation.

    Parameters
    ----------
    sequence1 : str
        First sequence (DNA or RNA, T → U conversion applied).
    sequence2 : str
        Second sequence (DNA or RNA, T → U conversion applied).
    no_ps : bool
        Suppress PostScript output (default True).

    Returns
    -------
    RNAduplexResult
        Sequences, structure, MFE, and hybridization positions.

    Raises
    ------
    ViennaRNAError
        If RNAduplex fails or output cannot be parsed.
    """
    rna1 = dna_to_rna(sequence1)
    rna2 = dna_to_rna(sequence2)
    tool_path = _get_tool_path("RNAduplex")

    args = [tool_path]
    if no_ps:
        args.append("--noPS")

    input_data = f"{rna1}\n{rna2}\n"

    try:
        result = subprocess.run(
            args,
            input=input_data,
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except subprocess.SubprocessError as e:
        raise ViennaRNAError(f"RNAduplex execution failed: {e}") from e

    return _parse_rnaduplex_output(rna1, rna2, result.stdout)


def _parse_rnaduplex_output(rna1: str, rna2: str, stdout: str) -> RNAduplexResult:
    """Parse RNAduplex stdout.

    Output format:
    .(((((((((...)))))))))...  &...((((((....))))))...  (-12.40)
    1,15  : 30,16
    """
    lines = stdout.strip().split("\n")

    structure = ""
    mfe = 0.0
    start1 = end1 = start2 = end2 = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith(">"):
            continue

        # Structure + MFE line: contains '&' and '(' ')'
        if "&" in line and ("(" in line or ")" in line):
            parts = line.rsplit(" ", 1)
            if len(parts) == 2:
                structure = parts[0]
                mfe_str = parts[1].strip("()")
                try:
                    mfe = float(mfe_str)
                except ValueError:
                    pass

        # Position line: "1,15  : 30,16"
        if ":" in line and "," in line:
            match = re.match(r"(\d+),(\d+)\s*:\s*(\d+),(\d+)", line)
            if match:
                start1 = int(match.group(1))
                end1 = int(match.group(2))
                start2 = int(match.group(3))
                end2 = int(match.group(4))

    if structure and mfe != 0.0:
        return RNAduplexResult(
            sequence1=rna1,
            sequence2=rna2,
            structure=structure,
            mfe=mfe,
            start1=start1,
            end1=end1,
            start2=start2,
            end2=end2,
        )

    raise ViennaRNAError(f"Could not parse RNAduplex output: {stdout[:200]}")
