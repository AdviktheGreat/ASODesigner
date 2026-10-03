#!/usr/bin/env python3
"""Individual module usage example.

Demonstrates using ASODesigner modules individually rather than
through the full pipeline. This is useful for custom workflows
or integration into existing pipelines.

Usage:
    python examples/module_usage.py
"""

from aso_designer.candidate_generator import generate_candidates
from aso_designer.config_loader import load_config
from aso_designer.filters import run_all_filters
from aso_designer.reporting import assemble_report
from aso_designer.scoring import score_rank_and_filter
from aso_designer.utils.sequence import gc_content, validate_sequence


def main() -> None:
    """Use individual modules to build a custom pipeline."""
    import os
    import tempfile

    # 1. Load configuration
    config = load_config(chemistry="gapmer")
    print("Configuration loaded:")
    print(f"  ASO length: {config.candidate.length}")
    print(f"  GC range: {config.candidate.gc_min}-{config.candidate.gc_max}")
    print()

    # 2. Define a target sequence
    target_sequence = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 10
    print(f"Target sequence: {len(target_sequence)} bp")

    # Validate the sequence
    validate_sequence(target_sequence, allow_rna=False)
    print(f"GC content: {gc_content(target_sequence):.2%}")
    print()

    # 3. Generate candidates
    candidates = generate_candidates(target_sequence, config.candidate)
    print(f"Generated {len(candidates)} candidates")

    # 4. Apply pre-filters
    run_all_filters(
        candidates,
        gc_min=config.candidate.gc_min,
        gc_max=config.candidate.gc_max,
        homopolymer_max=config.candidate.homopolymer_max,
    )
    passed = [c for c in candidates if c.passes_filters]
    print(f"Passed filters: {len(passed)} / {len(candidates)}")
    print()

    # 5. Score and rank (without computed features — uses neutral defaults)
    features_list = [{} for _ in passed]
    ranking = score_rank_and_filter(
        passed,
        features_list,
        min_score=0.0,
        top_n=5,
    )

    # 6. Display top candidates
    print("Top 5 Candidates:")
    print("-" * 70)
    for c in ranking.candidates:
        print(f"  Rank {c.rank}: {c.sequence} (pos {c.start}-{c.end}, score={c.final_score:.4f})")
    print()

    # 7. Generate reports
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = assemble_report(
            scored_candidates=ranking.candidates,
            output_dir=tmpdir,
            target_name="EXAMPLE",
            target_id="NM_EXAMPLE",
            chemistry="gapmer",
            total_candidates=len(candidates),
            passed_filters=len(passed),
            basename="example_report",
        )
        print(f"Reports generated: {list(paths.keys())}")
        for fmt, path in paths.items():
            print(f"  {fmt}: {os.path.basename(path)} ({os.path.getsize(path)} bytes)")


if __name__ == "__main__":
    main()
