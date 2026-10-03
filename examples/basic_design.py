#!/usr/bin/env python3
"""Basic ASO design example.

Demonstrates the simplest way to run the ASODesigner pipeline:
generate candidates from a local FASTA file, score them, and
generate reports.

Usage:
    python examples/basic_design.py
"""

from aso_designer.pipeline import run_pipeline


def main() -> None:
    """Run a basic ASO design pipeline."""
    # Create a test FASTA file (in practice, use a real transcript)
    import os
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write a test transcript
        fasta_path = os.path.join(tmpdir, "test_transcript.fasta")
        with open(fasta_path, "w") as f:
            f.write(">NM_TEST_001 Example transcript\n")
            f.write("ATGCATGCATGCATGCATGC" * 20)  # 400 bp transcript
            f.write("\n")

        output_dir = os.path.join(tmpdir, "results")

        # Run the pipeline
        result = run_pipeline(
            target="NM_TEST_001",
            chemistry="gapmer",
            local_fasta=fasta_path,
            output_dir=output_dir,
            top_n=10,
            min_score=0.0,
            # Skip optional stages for this example (no ViennaRNA installed)
            run_accessibility=False,
            run_thermodynamics=False,
            run_off_target=False,
            run_variant=False,
        )

        # Display results
        print("\nPipeline complete!")
        print(f"  Target:          {result.target}")
        print(f"  Transcript:      {result.transcript_id} ({result.transcript_length} bp)")
        print(f"  Total candidates: {result.total_candidates}")
        print(f"  Passed filters:  {result.passed_filters}")
        print(f"  Ranked:          {result.ranked_candidates}")
        print(f"  Top score:       {result.top_score:.4f}")
        print(f"  Stages run:      {result.stages_run}")
        if result.errors:
            print(f"  Errors:          {result.errors}")
        print(f"  Reports:         {result.report_paths}")


if __name__ == "__main__":
    main()
