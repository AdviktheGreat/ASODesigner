# ASODesigner Methodology

This document describes the scientific approach behind ASODesigner's
candidate design, scoring, and ranking pipeline.

## Overview

ASODesigner implements a multi-stage computational pipeline for antisense
oligonucleotide (ASO) candidate design. The pipeline generates all
possible ASO candidates of a specified length across a target transcript,
filters out candidates with poor sequence properties, computes
thermodynamic and accessibility features, screens for off-target
binding, annotates variant overlaps, and ranks candidates by a
transparent multi-factor composite score.

### Design principles

1. **Transparency** — Every score component is documented with its
   weight, formula, and scientific basis. No black-box scoring.
2. **Configurability** — All parameters (length, GC bounds, weights,
   chemistry presets) are user-configurable via YAML.
3. **Graceful degradation** — Optional stages (accessibility,
   off-target, variant) can be skipped without breaking the pipeline.
4. **Reproducibility** — All external calls (NCBI, ViennaRNA) are
   cached and mockable for testing.

## Pipeline stages

### 1. Transcript fetching

The pipeline retrieves the target transcript sequence from one of:

- **NCBI RefSeq** via the Entrez E-utilities API (ESearch + EFetch)
- **Local FASTA file** (user-provided)
- **Disk cache** (previous fetches, with configurable TTL)

The `FetchManager` tries sources in order: cache → NCBI → local file.
Transcripts are validated for ambiguous bases, minimum length, and ORF
detection (CDS annotation).

### 2. Candidate generation

A sliding window of configurable length (default 20 nt, range 10–30)
is moved across the transcript with a configurable step size (default 1).
Each window becomes an `ASOCandidate` with:

- Sequence (DNA alphabet, uppercase)
- Start/end positions (0-based transcript coordinates)
- GC content
- Region annotation (5UTR, CDS, 3UTR) if transcript metadata is available
- Chemistry type (gapmer, mixmer, 2'-OMe, MOE, morpholino)

### 3. Pre-filters

Candidates are screened by four filters before expensive feature
computation:

- **GC content** — Reject candidates outside [0.30, 0.60] (configurable)
- **Homopolymer** — Reject candidates with >4 consecutive identical bases
- **End stability** — Check GC content at 5' and 3' ends (3 nt windows)
- **Ambiguous bases** — Reject any IUPAC ambiguity codes (N, R, Y, etc.)

Failed candidates are retained with `passes_filters=False` and a
`filter_reason` in their features dict, allowing post-hoc inspection.

### 4. RNA accessibility

Target RNA secondary structure is predicted using ViennaRNA:

- **RNAfold** — Full-transcript MFE structure
- **RNAplfold** — Local accessibility profile (per-position unpaired
  probabilities with configurable window size)

Accessibility scores are computed by aggregating unpaired probabilities
across the ASO binding window. A score of 1.0 means the target site is
fully accessible (unpaired); 0.0 means fully structured.

Results are cached to disk to avoid expensive re-computation.

### 5. Thermodynamics

See [Thermodynamics Methodology](#thermodynamics) below for details.

### 6. Off-target analysis

See [Off-target Methodology](#off-target-analysis) below for details.

### 7. Variant/SNP awareness

See [Variant Methodology](#variant-snp-awareness) below for details.

### 8. Scoring and ranking

See [Scoring Methodology](#scoring) below for details.

### 9. Reporting

The pipeline generates two report formats:

- **JSON** — Structured, machine-readable report with summary
  statistics, candidate details, and optional detail sections
- **HTML** — Self-contained (no external dependencies) report with
  embedded CSS, summary cards, candidate table with score bars,
  SVG score distribution chart, and detail tables

Reports are written to a configurable output directory.

## Thermodynamics

### Nearest-neighbor (NN) model

ASODesigner uses the unified nearest-neighbor (NN) model for melting
temperature (Tm) and free energy (ΔG) predictions:

- **DNA/DNA** — SantaLucia 1998 unified parameters (PNAS 95:1460-1465)
- **RNA/RNA** — Turner 2004 parameters (Nucleic Acids Res 38:D280-D282)

The NN model computes thermodynamic properties by summing contributions
from each dinucleotide step in the duplex, plus initiation and symmetry
corrections.

### Melting temperature (Tm)

Tm is calculated as:

$$T_m = \frac{\Delta H \times 1000}{\Delta S + R \times \ln(C_T / x)} - 273.15$$

Where:
- $\Delta H$ = total enthalpy (kcal/mol), summed over NN steps
- $\Delta S$ = total entropy (cal/mol·K), summed over NN steps
- $R$ = gas constant (1.987 cal/mol·K)
- $C_T$ = total strand concentration
- $x$ = 4 for non-self-complementary, 1 for self-complementary

#### Salt correction

The SantaLucia salt correction adjusts ΔS rather than Tm directly:

$$\Delta S_{corrected} = \Delta S + 0.368 \times (N-1) \times \ln([Na^+])$$

Where $N$ is the number of nucleotides and $[Na^+]$ is the sodium
concentration in molar.

#### Chemistry-aware Tm

Each chemistry has preset Tm adjustments (e.g., 2'-OMe increases Tm
due to increased binding affinity). These are applied as additive
corrections to the computed Tm.

### Duplex ΔG

Duplex formation free energy is computed at a reference temperature
(37°C = 310.15 K):

$$\Delta G = \Delta H - T_{ref} \times \Delta S / 1000$$

Two methods are available:
- **Nearest-neighbor** — Fast, uses the NN parameter tables
- **RNAduplex** — ViennaRNA tool, more accurate for structured targets

### Self-structure ΔG

Self-structure (intramolecular folding) of the ASO is computed via:
- **RNAfold** (ViennaRNA) — MFE structure prediction
- **Nearest-neighbor** — Approximate using NN parameters

A chemistry-aware self-structure penalty is applied based on the
self-structure ΔG and the chemistry preset's self-structure weight.

### Target structure ΔG

The target RNA structure at the binding site is computed to estimate
the energetic cost of unfolding the target for ASO binding.

### Net binding energy

The net binding energy combines all thermodynamic terms:

$$\Delta G_{net} = \Delta G_{duplex} - \Delta G_{self} - \Delta G_{target}$$

A more negative $\Delta G_{net}$ indicates stronger net binding.

### Binding classification

Candidates are classified by their thermodynamic profile:

| Class | Criteria |
|---|---|
| Strong | $\Delta G_{net} < -10$ kcal/mol and Tm in optimal range |
| Moderate | $-10 \leq \Delta G_{net} < -5$ kcal/mol |
| Weak | $\Delta G_{net} \geq -5$ kcal/mol |

## Off-target analysis

### Transcriptome reference

ASODesigner builds a `TranscriptomeReference` from a set of transcript
sequences. A k-mer index (default k=12) is constructed for fast
seed-based searching. Each transcript is stored with its accession,
gene name, and sequence.

### Seed-and-extend search

For each ASO candidate, the search engine:

1. **Seed extraction** — Extracts k-mers from the ASO's seed region
   (default: first 8 nt, the most critical for RNase H recruitment)
2. **Index lookup** — Finds transcriptome positions matching any seed
   k-mer
3. **Extension** — Extends each seed match to full-length alignment,
   counting mismatches in seed and non-seed regions
4. **Classification** — Classifies hits as perfect, seed-match, or
   seed-mismatch

### Hit filtering

Hits are filtered by:
- Maximum total mismatches (default 3)
- Maximum seed mismatches (default 0 for strict seed filtering)
- On-target exclusion (the intended target is excluded from
  off-target counts)

### Off-target scoring

The off-target specificity score is computed as:

$$\text{score} = \frac{1}{1 + \text{weighted\_hits}}$$

Where `weighted_hits` sums penalties for each off-target hit:
- Perfect match: weight 1.0
- Seed match (0 seed mismatches, >0 total): weight 0.5
- Seed mismatch: weight 0.2

A score of 1.0 means no off-target hits (highly specific); 0.0 means
many perfect-match off-targets.

### Specificity classification

| Class | Score Range |
|---|---|
| Specific | ≥ 0.8 |
| Moderate | 0.5–0.8 |
| Non-specific | < 0.5 |

## Variant/SNP awareness

### Variant data sources

ASODesigner fetches variant data from:

1. **dbSNP** — NCBI ESearch/ESummary API for known SNPs
2. **ClinVar** — NCBI API for clinical significance
3. **VCF files** — User-provided variant call format files

The `VariantFetchManager` tries API sources first, then falls back to
VCF. All variants are deduplicated by position and allele.

### Variant annotation

For each ASO candidate, the `VariantAnnotator` checks for overlapping
variants:

- **Position overlap** — Does any variant fall within the ASO binding
  site (start to end)?
- **Seed region** — Does any variant fall within the seed region
  (first 8 nt)? Seed-region variants are penalized more heavily.
- **Pathogenic flag** — Is any overlapping variant marked as
  pathogenic or likely pathogenic in ClinVar?

### Variant penalty

The variant penalty combines:

- **MAF-weighted penalty** — Common variants (high MAF) receive a
  smaller penalty than rare variants
- **Seed-region multiplier** — Variants in the seed region receive a
  2× penalty multiplier
- **Pathogenic bonus** — Pathogenic variants add a fixed penalty
  regardless of MAF

### Risk classification

| Risk Level | Criteria |
|---|---|
| High risk | Pathogenic variant in seed region |
| Moderate risk | Variant overlap in seed region, or pathogenic in non-seed |
| Low risk | Variant overlap in non-seed region only |
| No risk | No overlapping variants |

## Scoring

### Score components

ASODesigner uses a pluggable component architecture for scoring. Each
component implements the `ScoreComponent` protocol and returns a
sub-score in [0.0, 1.0]. The `ScoreEngine` combines components into a
weighted composite score.

### Component details

#### 1. Accessibility (weight: 0.25)

Uses the RNAplfold-based accessibility score directly. If accessibility
was not computed, returns a neutral 0.5.

$$\text{score}_{accessibility} = \text{accessibility\_score}$$

#### 2. Thermodynamics (weight: 0.20)

Combines Tm score and binding class:

$$\text{score}_{thermo} = 0.6 \times \text{Tm\_score} + 0.4 \times \text{binding\_class\_score}$$

Binding class scores: strong=1.0, moderate=0.5, weak=0.0.

#### 3. GC content (weight: 0.10)

Optimal range is 30–60%. Within range: 1.0. Outside: linear decrease
to 0.0 at GC=0.0 or GC=1.0.

#### 4. Off-target (weight: 0.20)

Uses the off-target specificity score directly. If not computed,
returns neutral 0.5.

#### 5. Variant (weight: 0.10)

Converts variant penalty to a score:

$$\text{score}_{variant} = 1.0 - \text{variant\_penalty}$$

Clamped to [0, 1]. If no variant data, returns 1.0 (no penalty).

#### 6. Self-structure (weight: 0.10)

Converts self-structure penalty to a score:

$$\text{score}_{self} = 1.0 - \text{self\_structure\_penalty}$$

Clamped to [0, 1]. If no self-structure data, returns 1.0.

#### 7. Positional (weight: 0.05)

Region-based base score with center-of-transcript position bonus:

| Region | Base Score |
|---|---|
| CDS | 1.0 |
| 5UTR | 0.7 |
| 3UTR | 0.5 |
| noncoding | 0.6 |
| unknown | 0.3 |

A position bonus (up to 0.1) is added for candidates near the center
of the transcript, decreasing linearly to 0 at the ends.

### Composite score

$$\text{composite} = \sum_{i} w_i \times \text{score}_i$$

$$\text{final} = \text{clamp}(\text{composite} + \text{chemistry\_bonus}, 0, 1)$$

Where $w_i$ are the component weights (summing to 1.0) and
`chemistry_bonus` is an optional per-chemistry additive bonus.

### Ranking and filtering

Candidates are ranked by final score (descending). The pipeline
supports:

- **Minimum score threshold** — Candidates below `min_score` are removed
- **Top-N selection** — Keep only the top N candidates
- **Rank assignment** — Each candidate gets a rank (1 = best)

### Chemistry-aware adjustments

Each chemistry preset can adjust:
- Tm (additive correction)
- Self-structure weight (multiplier on self-structure penalty)
- Scoring bonus (additive to composite score)

These are applied automatically based on the `--chemistry` flag.
