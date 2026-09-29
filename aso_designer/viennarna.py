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
from dataclasses import dataclass

from aso_designer.utils.exceptions import ViennaRNAError
from aso_designer.utils.logger import get_logger

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
