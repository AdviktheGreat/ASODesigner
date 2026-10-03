#!/usr/bin/env python3
"""Batch design example.

Demonstrates running the ASO design pipeline for multiple targets
and generating per-target reports.

Usage:
    python examples/batch_design.py
"""

import os
import tempfile

from aso_designer.pipeline import run_pipeline


def main() -> None:
    """Run batch design for multiple targets."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test FASTA files for multiple targets
        targets = [
            ("NM_GENE_A", "ATGCATGCATGCATGCATGC" * 20),
            ("NM_GENE_B", "GCATGCATGCATGCATGCAT" * 20),
            ("NM_GENE_C", "TGCATGCATGCATGCATGCAT" * 20),
        ]

        results = []

        for accession, sequence in targets:
            fasta_path = os.path.join(tmpdir, f"{accession}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{accession} Test gene\n{sequence}\n")

            output_dir = os.path.join(tmpdir, "results", accession)

            result = run_pipeline(
                target=accession,
                chemistry="gapmer",
                local_fasta=fasta_path,
                output_dir=output_dir,
                top_n=5,
                run_accessibility=False,
                run_thermodynamics=False,
                run_off_target=False,
                run_variant=False,
            )
            results.append(result)

        # Display batch summary
        print("Batch Design Results")
        print("=" * 70)
        print(f"{'Target':<15} {'Length':<8} {'Candidates':<12} {'Ranked':<8} {'Top Score':<10}")
        print("-" * 70)
        for r in results:
            print(
                f"{r.target:<15} "
                f"{r.transcript_length:<8} "
                f"{r.total_candidates:<12} "
                f"{r.ranked_candidates:<8} "
                f"{r.top_score:<10.4f}"
            )

        print(f"\nTotal targets: {len(results)}")
        print(f"Successful:    {sum(1 for r in results if r.success)}")
        print(f"Total candidates: {sum(r.total_candidates for r in results)}")
        print(f"Total ranked:     {sum(r.ranked_candidates for r in results)}")


if __name__ == "__main__":
    main()
