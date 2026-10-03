"""Command-line interface for ASODesigner.

Provides a command-based CLI using argparse with subcommands for
the full ASO design pipeline, scoring pre-generated candidates,
generating reports, and displaying configuration info.

Designed from scratch — GuideX uses a simple argparse setup; we
add subcommands, shared argument groups, verbose logging control,
and structured output formatting.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from aso_designer import __version__
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Argument parser construction
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with subcommands.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser with 'design', 'score', 'report', and
        'info' subcommands.
    """
    parser = argparse.ArgumentParser(
        prog="aso-designer",
        description=(
            "ASODesigner — in-silico antisense oligonucleotide (ASO) "
            "candidate design and ranking pipeline."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  aso-designer design --target BRCA1 --chemistry gapmer\n"
            "  aso-designer score --candidates candidates.json\n"
            "  aso-designer report --results ranked.json --output-dir reports/\n"
            "  aso-designer info --chemistry gapmer\n"
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"ASODesigner v{__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        help="Available commands",
    )

    _add_design_subparser(subparsers)
    _add_score_subparser(subparsers)
    _add_report_subparser(subparsers)
    _add_info_subparser(subparsers)

    return parser


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments common to multiple subcommands."""
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to a YAML config file (overrides defaults).",
    )
    parser.add_argument(
        "--chemistry",
        type=str,
        default="gapmer",
        choices=["gapmer", "mixmer", "2ome", "moe", "morpholino", "full"],
        help="ASO chemistry (default: gapmer).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG-level) logging.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress all output except errors.",
    )


def _add_design_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Add the 'design' subcommand."""
    parser = subparsers.add_parser(
        "design",
        help="Run the full ASO design pipeline from target to ranked report.",
        description=(
            "Run the complete pipeline: fetch transcript, generate "
            "candidates, filter, compute features, score, rank, and "
            "generate reports."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(parser)

    parser.add_argument(
        "--target",
        type=str,
        required=True,
        help="Target gene symbol (e.g. BRCA1) or transcript accession (e.g. NM_007294).",
    )
    parser.add_argument(
        "--organism",
        type=str,
        default="Homo sapiens",
        help="Organism for NCBI lookup (default: Homo sapiens).",
    )
    parser.add_argument(
        "--local-fasta",
        type=str,
        default=None,
        help="Path to a local FASTA file (bypasses NCBI).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Output directory for reports (default: results).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Keep only the top N candidates.",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="Minimum score threshold (default: 0.0).",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip HTML report generation.",
    )
    parser.add_argument(
        "--no-offtarget",
        action="store_true",
        help="Skip off-target analysis (faster but less specific).",
    )
    parser.add_argument(
        "--no-variant",
        action="store_true",
        help="Skip variant/SNP annotation.",
    )


def _add_score_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Add the 'score' subcommand."""
    parser = subparsers.add_parser(
        "score",
        help="Score pre-generated candidates from a JSON file.",
        description=(
            "Load candidates from a JSON file (from candidate_serializer), "
            "score them using the scoring engine, and output ranked results."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(parser)

    parser.add_argument(
        "--candidates",
        type=str,
        required=True,
        help="Path to a JSON file of candidates (from save_candidates).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON file for scored candidates (default: stdout).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Keep only the top N candidates.",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="Minimum score threshold (default: 0.0).",
    )


def _add_report_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Add the 'report' subcommand."""
    parser = subparsers.add_parser(
        "report",
        help="Generate JSON and HTML reports from saved results.",
        description=("Load a JSON results file and generate formatted JSON and HTML reports."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(parser)

    parser.add_argument(
        "--results",
        type=str,
        required=True,
        help="Path to a JSON results file (from 'design' or 'score').",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="reports",
        help="Output directory for reports (default: reports).",
    )
    parser.add_argument(
        "--basename",
        type=str,
        default="aso_report",
        help="Base filename for reports (default: aso_report).",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip HTML report generation.",
    )


def _add_info_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Add the 'info' subcommand."""
    parser = subparsers.add_parser(
        "info",
        help="Display configuration, chemistry presets, and version info.",
        description=(
            "Show ASODesigner version, default configuration, and chemistry preset details."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(parser)

    parser.add_argument(
        "--show-weights",
        action="store_true",
        help="Show scoring weights for the selected chemistry.",
    )
    parser.add_argument(
        "--show-presets",
        action="store_true",
        help="Show all available chemistry presets.",
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    """ASODesigner CLI main entry point.

    Parameters
    ----------
    argv : Sequence[str] or None
        Command-line arguments. If None, uses sys.argv.

    Returns
    -------
    int
        Exit code (0 = success, 1 = error).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 0

    # Configure logging level
    if getattr(args, "verbose", False):
        _configure_logging("DEBUG")
    elif getattr(args, "quiet", False):
        _configure_logging("ERROR")
    else:
        _configure_logging("INFO")

    logger.debug("Parsed args: %s", args)

    # Dispatch to subcommand handler
    handlers = {
        "design": _handle_design,
        "score": _handle_score,
        "report": _handle_report,
        "info": _handle_info,
    }

    handler = handlers.get(args.command)
    if handler is None:
        parser.print_help()
        return 1

    try:
        return handler(args)
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        return 130
    except Exception as e:
        logger.error("Command '%s' failed: %s", args.command, e)
        if getattr(args, "verbose", False):
            raise
        return 1


def _configure_logging(level: str) -> None:
    """Configure the root logger level."""
    import logging

    logging.getLogger("aso_designer").setLevel(getattr(logging, level.upper(), logging.INFO))


# ---------------------------------------------------------------------------
# Subcommand handlers (stubs — implemented in Step 43)
# ---------------------------------------------------------------------------


def _handle_design(args: argparse.Namespace) -> int:
    """Handle the 'design' subcommand."""
    logger.info("Design command not yet implemented")
    print("Design command is not yet implemented. See Step 43.")
    return 0


def _handle_score(args: argparse.Namespace) -> int:
    """Handle the 'score' subcommand."""
    logger.info("Score command not yet implemented")
    print("Score command is not yet implemented. See Step 43.")
    return 0


def _handle_report(args: argparse.Namespace) -> int:
    """Handle the 'report' subcommand."""
    logger.info("Report command not yet implemented")
    print("Report command is not yet implemented. See Step 43.")
    return 0


def _handle_info(args: argparse.Namespace) -> int:
    """Handle the 'info' subcommand."""
    logger.info("Info command not yet implemented")
    print("Info command is not yet implemented. See Step 43.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
