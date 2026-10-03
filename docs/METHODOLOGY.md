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
