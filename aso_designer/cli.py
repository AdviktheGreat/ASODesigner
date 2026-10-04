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
        "--transcriptome",
        type=str,
        default=None,
        help="Path to a transcriptome FASTA for off-target search "
        "(overrides off_target.transcriptome_path in the config).",
    )
    parser.add_argument(
        "--vcf",
        type=str,
        default=None,
        help="Path to a VCF file with known variants in transcript "
        "coordinates (overrides variant.vcf_path in the config).",
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
# CLI utility functions
# ---------------------------------------------------------------------------


def format_score_table(scored_candidates: list, top_n: int = 10) -> str:
    """Format scored candidates as a text table for CLI output.

    Parameters
    ----------
    scored_candidates : list[CandidateScore]
        Scored and ranked candidates.
    top_n : int
        Maximum rows to display (default 10).

    Returns
    -------
    str
        Formatted text table.
    """
    if not scored_candidates:
        return "No candidates to display."

    lines: list[str] = []
    header = (
        f"{'Rank':<6} {'Sequence':<22} {'Start':<8} {'End':<8} "
        f"{'Chem':<10} {'Score':<8} {'Binding':<10}"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for c in scored_candidates[:top_n]:
        rank = str(c.rank or "-")
        seq = c.sequence[:20] + ("..." if len(c.sequence) > 20 else "")
        binding_str = ""

        lines.append(
            f"{rank:<6} {seq:<22} {c.start:<8} {c.end:<8} "
            f"{c.chemistry:<10} {c.final_score:<8.4f} {binding_str:<10}"
        )

    return "\n".join(lines)


def format_summary(stats: dict, indent: int = 0) -> str:
    """Format a summary dict as aligned key-value pairs.

    Parameters
    ----------
    stats : dict
        Summary statistics (from ReportSummary.to_dict()).
    indent : int
        Indentation level (default 0).

    Returns
    -------
    str
        Formatted text.
    """
    prefix = " " * indent
    lines: list[str] = []
    max_key = max(len(k) for k in stats) if stats else 0
    for key, value in stats.items():
        if isinstance(value, float):
            val_str = f"{value:.4f}"
        else:
            val_str = str(value)
        lines.append(f"{prefix}{key:<{max_key + 2}} {val_str}")
    return "\n".join(lines)


def print_banner(title: str, width: int = 60) -> None:
    """Print a section banner to stdout."""
    print()
    print("=" * width)
    print(f"  {title}")
    print("=" * width)


def validate_output_dir(path: str) -> str:
    """Validate and normalize an output directory path.

    Parameters
    ----------
    path : str
        Output directory path.

    Returns
    -------
    str
        Normalized path.
    """
    from pathlib import Path

    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return str(p.resolve())


def validate_input_file(path: str, description: str = "input file") -> str:
    """Validate that an input file exists and is readable.

    Parameters
    ----------
    path : str
        Input file path.
    description : str
        Human-readable description for error messages.

    Returns
    -------
    str
        Normalized path.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    from pathlib import Path

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{description} not found: {path}")
    if not p.is_file():
        raise FileNotFoundError(f"{description} is not a file: {path}")
    return str(p.resolve())


def validate_target(target: str) -> str:
    """Validate a target gene symbol or transcript accession.

    Parameters
    ----------
    target : str
        Gene symbol (e.g. "BRCA1") or accession (e.g. "NM_007294").

    Returns
    -------
    str
        Normalized target string (stripped, uppercased).

    Raises
    ------
    ValueError
        If the target is empty or contains invalid characters.
    """
    if not target or not target.strip():
        raise ValueError("Target cannot be empty.")
    normalized = target.strip().upper()
    # Allow alphanumeric, dots, underscores, dashes (for accessions)
    import re

    if not re.match(r"^[A-Z0-9._-]+$", normalized):
        raise ValueError(
            f"Invalid target '{target}': use gene symbol (BRCA1) or accession (NM_007294)."
        )
    return normalized


def is_accession(target: str) -> bool:
    """Check if a target string looks like an NCBI accession.

    Accessions match patterns like NM_001123, NR_123456, XM_456789.

    Parameters
    ----------
    target : str
        Target string.

    Returns
    -------
    bool
        True if the target looks like an accession.
    """
    import re

    return bool(re.match(r"^[A-Z]{1,4}_\d+(\.\d+)?$", target.strip().upper()))


# ---------------------------------------------------------------------------
# CLI progress and result display
# ---------------------------------------------------------------------------


def print_progress(stage: str, message: str) -> None:
    """Print a pipeline stage progress message.

    Parameters
    ----------
    stage : str
        Stage name (e.g. "FETCH", "FILTER", "SCORE").
    message : str
        Human-readable progress message.
    """
    print(f"  [{stage:>8}] {message}")


def print_summary_box(title: str, stats: dict) -> None:
    """Print a summary box with statistics.

    Parameters
    ----------
    title : str
        Box title.
    stats : dict
        Statistics to display.
    """
    print()
    print(f"  ┌─ {title} " + "─" * max(0, 50 - len(title)) + "┐")
    for key, value in stats.items():
        if isinstance(value, float):
            val_str = f"{value:.4f}"
        else:
            val_str = str(value)
        key_str = key.replace("_", " ").title()
        print(f"  │ {key_str:<24} {val_str:>20} │")
    print("  └" + "─" * 53 + "┘")


def print_results_table(scored: list, top_n: int = 10) -> None:
    """Print scored candidates as a formatted table to stdout.

    Parameters
    ----------
    scored : list[CandidateScore]
        Scored and ranked candidates.
    top_n : int
        Maximum rows to display (default 10).
    """
    print()
    print(format_score_table(scored, top_n=top_n))
    if len(scored) > top_n:
        print(f"  ... and {len(scored) - top_n} more candidates")


def print_error(message: str) -> None:
    """Print an error message to stderr."""
    print(f"Error: {message}", file=sys.stderr)


def print_success(message: str) -> None:
    """Print a success message to stdout."""
    print(f"  ✓ {message}")


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def _handle_design(args: argparse.Namespace) -> int:
    """Handle the 'design' subcommand — full pipeline."""
    from aso_designer.pipeline import ASOPipeline, PipelineConfig

    # Validate target
    try:
        target = validate_target(args.target)
    except ValueError as e:
        print_error(str(e))
        return 1

    # Validate local FASTA if provided
    local_fasta = None
    if args.local_fasta:
        try:
            local_fasta = validate_input_file(args.local_fasta, "FASTA file")
        except FileNotFoundError as e:
            print_error(str(e))
            return 1

    # Build pipeline config
    config = PipelineConfig.from_args(args)
    config.target = target
    config.local_fasta = local_fasta

    print_banner("ASODesigner — Design Pipeline")
    print_progress("CONFIG", f"Target: {target}, Chemistry: {config.chemistry}")
    print_progress("CONFIG", f"Output: {config.output_dir}")

    # Run pipeline
    pipeline = ASOPipeline(config)
    result = pipeline.run()

    # Display results
    print_banner("Pipeline Results")
    print_summary_box(
        "Summary",
        {
            "target": result.target,
            "transcript_id": result.transcript_id,
            "transcript_length": result.transcript_length,
            "total_candidates": result.total_candidates,
            "passed_filters": result.passed_filters,
            "ranked_candidates": result.ranked_candidates,
            "top_score": result.top_score,
        },
    )

    if result.errors:
        print()
        print("  Warnings/Errors:")
        for err in result.errors:
            print(f"    ! {err}")

    if result.report_paths:
        print()
        print("  Reports generated:")
        for fmt, path in result.report_paths.items():
            print(f"    {fmt.upper():>6}: {path}")

    if result.success:
        print_success("Pipeline completed successfully.")
        return 0
    else:
        print_error("Pipeline did not produce ranked candidates.")
        return 1


def _handle_score(args: argparse.Namespace) -> int:
    """Handle the 'score' subcommand — score pre-generated candidates."""
    import json as json_module

    from aso_designer.candidate import ASOCandidate
    from aso_designer.reporting import build_report_data, generate_report
    from aso_designer.scoring import rank_and_filter, score_candidates_batch

    # Validate input file
    try:
        candidates_path = validate_input_file(args.candidates, "Candidates file")
    except FileNotFoundError as e:
        print_error(str(e))
        return 1

    print_banner("ASODesigner — Score Candidates")

    # Load candidates
    print_progress("LOAD", f"Reading candidates from {candidates_path}")
    with open(candidates_path) as f:
        data = json_module.loads(f.read())

    candidates_data = data.get("candidates", data if isinstance(data, list) else [])
    candidates = [ASOCandidate.from_dict(c) for c in candidates_data]
    print_progress("LOAD", f"Loaded {len(candidates)} candidates")

    if not candidates:
        print_error("No candidates found in file.")
        return 1

    # Score
    print_progress("SCORE", f"Scoring with chemistry: {args.chemistry}")
    features_list = [{} for _ in candidates]
    scored = score_candidates_batch(candidates, features_list)

    # Rank and filter
    ranking = rank_and_filter(
        scored,
        min_score=args.min_score,
        top_n=args.top_n,
    )
    print_progress(
        "RANK",
        f"Ranked {ranking.ranked_candidates} candidates (top score: {ranking.top_score:.4f})",
    )

    # Display results
    print_results_table(ranking.candidates, top_n=20)

    # Output
    if args.output:
        report = build_report_data(
            scored_candidates=ranking.candidates,
            total_candidates=len(candidates),
            chemistry=args.chemistry,
        )
        from pathlib import Path

        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        paths = generate_report(
            report,
            out_path.parent,
            basename=out_path.stem,
            generate_html=False,
        )
        print_success(f"Scored candidates saved to {paths['json']}")
    else:
        # Print to stdout as JSON
        print()
        print_banner("Scored Candidates (JSON)")
        output_data = {
            "total_candidates": len(candidates),
            "ranked_candidates": ranking.ranked_candidates,
            "top_score": round(ranking.top_score, 4),
            "candidates": [c.to_dict() for c in ranking.candidates],
        }
        print(json_module.dumps(output_data, indent=2))

    print_success("Scoring complete.")
    return 0


def _handle_report(args: argparse.Namespace) -> int:
    """Handle the 'report' subcommand — generate reports from saved results."""
    from aso_designer.reporting import generate_report, load_json_report

    # Validate input file
    try:
        results_path = validate_input_file(args.results, "Results file")
    except FileNotFoundError as e:
        print_error(str(e))
        return 1

    print_banner("ASODesigner — Generate Reports")

    # Load results
    print_progress("LOAD", f"Reading results from {results_path}")
    report_data = load_json_report(results_path)
    print_progress("LOAD", f"Loaded {len(report_data.candidates)} candidates")

    if not report_data.candidates:
        print_error("No candidates found in results file.")
        return 1

    # Generate reports
    output_dir = validate_output_dir(args.output_dir)
    print_progress("REPORT", f"Generating reports in {output_dir}")

    paths = generate_report(
        report_data,
        output_dir,
        basename=args.basename,
        generate_html=not args.no_html,
    )

    print()
    print("  Reports generated:")
    for fmt, path in paths.items():
        print(f"    {fmt.upper():>6}: {path}")

    print_success("Report generation complete.")
    return 0


def _handle_info(args: argparse.Namespace) -> int:
    """Handle the 'info' subcommand — display config, presets, and version."""
    from aso_designer import __author__, __email__, __version__
    from aso_designer.config_loader import load_config

    print_banner("ASODesigner — Info")

    # Version info
    print(f"  Version:    {__version__}")
    print(f"  Author:     {__author__}")
    print(f"  Email:      {__email__}")
    print()

    # Load config
    config = load_config(args.config, chemistry=args.chemistry)

    # Show candidate config
    print("  Candidate Configuration:")
    cand = config.candidate
    print(f"    Length:         {cand.length} nt")
    print(f"    Step:           {cand.step}")
    print(f"    GC min:         {cand.gc_min}")
    print(f"    GC max:         {cand.gc_max}")
    print(f"    Homopolymer max: {cand.homopolymer_max}")
    print(f"    Chemistry:      {cand.chemistry}")
    print()

    # Show scoring weights
    if args.show_weights:
        print("  Scoring Weights:")
        weights = config.scoring.weights
        for name, weight in weights.model_dump().items():
            print(f"    {name:<20} {weight:.2f}")
        print()

    # Show chemistry presets
    if args.show_presets:
        print("  Chemistry Presets:")
        from aso_designer.config_loader import _CHEMISTRY_PRESETS_PATH, load_yaml

        if _CHEMISTRY_PRESETS_PATH.exists():
            presets = load_yaml(_CHEMISTRY_PRESETS_PATH)
            for chem_name, preset in presets.items():
                print(f"    {chem_name}:")
                for key, value in preset.items():
                    print(f"      {key:<24} {value}")
                print()
        else:
            print("    (presets file not found)")
        print()

    # Show reporting config
    print("  Reporting Configuration:")
    rpt = config.reporting
    print(f"    Output dir:           {rpt.output_dir}")
    print(f"    Generate HTML:        {rpt.generate_html}")
    print(f"    Off-target details:   {rpt.include_offtarget_details}")
    print()

    print_success("Info displayed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
