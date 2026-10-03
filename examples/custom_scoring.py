#!/usr/bin/env python3
"""Custom scoring example.

Demonstrates how to create a custom scoring engine with non-default
weights and use it to score candidates.

Usage:
    python examples/custom_scoring.py
"""

from aso_designer.candidate import ASOCandidate
from aso_designer.scoring import (
    ScoreContext,
    create_score_engine,
    rank_and_filter,
    score_candidates_batch,
)


def main() -> None:
    """Score candidates with a custom scoring engine."""
    # Create some test candidates
    candidates = [
        ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, gc_content=0.5, region="CDS"),
        ASOCandidate("GCGCGCGCGCGCGCGCGCGC", 20, 40, 20, gc_content=1.0, region="3UTR"),
        ASOCandidate("AATTAATTAATTAATTAATT", 40, 60, 20, gc_content=0.1, region="5UTR"),
        ASOCandidate("ATGCATGCATGCATGCATGC", 60, 80, 20, gc_content=0.5, region="CDS"),
    ]

    # Create a custom engine with modified weights
    # Emphasize accessibility and off-target specificity
    engine = create_score_engine(
        weights={
            "accessibility": 0.30,
            "thermodynamics": 0.15,
            "gc_content": 0.10,
            "off_target": 0.25,
            "variant": 0.10,
            "self_structure": 0.05,
            "positional": 0.05,
        },
        chemistry_bonus=0.02,  # Small bonus for gapmer chemistry
    )

    # Build features for each candidate
    features_list = [
        {
            "accessibility_score": 0.8,
            "tm_score": 0.9,
            "binding_class": "strong",
            "off_target_score": 0.95,
            "variant_penalty": 0.1,
            "variant_risk": "low_risk",
            "self_structure_penalty": 0.3,
            "dg_self": -2.0,
            "transcript_length": 2000,
        },
        {
            "accessibility_score": 0.3,
            "tm_score": 0.5,
            "binding_class": "moderate",
            "off_target_score": 0.7,
            "variant_penalty": 0.0,
            "variant_risk": "",
            "self_structure_penalty": 0.6,
            "dg_self": -5.0,
            "transcript_length": 2000,
        },
        {
            "accessibility_score": 0.6,
            "tm_score": 0.7,
            "binding_class": "moderate",
            "off_target_score": 0.9,
            "variant_penalty": 0.0,
            "variant_risk": "",
            "self_structure_penalty": 0.1,
            "dg_self": -0.5,
            "transcript_length": 2000,
        },
        {
            "accessibility_score": 0.9,
            "tm_score": 0.85,
            "binding_class": "strong",
            "off_target_score": 0.98,
            "variant_penalty": 0.0,
            "variant_risk": "",
            "self_structure_penalty": 0.2,
            "dg_self": -1.5,
            "transcript_length": 2000,
        },
    ]

    # Score candidates
    scored = score_candidates_batch(candidates, features_list, engine)

    # Rank and filter
    ranking = rank_and_filter(scored, min_score=0.0, top_n=10)

    # Display results
    print("Custom Scoring Results")
    print("=" * 70)
    print(f"{'Rank':<6} {'Sequence':<22} {'Score':<8} {'Region':<8}")
    print("-" * 70)
    for c in ranking.candidates:
        print(f"{c.rank:<6} {c.sequence:<22} {c.final_score:<8.4f} {c.sequence[:0]:<8}")

    print(f"\nTop score:    {ranking.top_score:.4f}")
    print(f"Median score: {ranking.median_score:.4f}")
    print(f"Total:        {ranking.total_candidates}")
    print(f"Ranked:       {ranking.ranked_candidates}")

    # Show component breakdown for the top candidate
    if ranking.candidates:
        top = ranking.candidates[0]
        print(f"\nTop candidate component breakdown:")
        for name, score in top.component_scores.items():
            print(f"  {name:<20} {score:.4f}")


if __name__ == "__main__":
    main()
