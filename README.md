# ASODesigner

**In-silico antisense oligonucleotide (ASO) candidate design and ranking.**

ASODesigner is a computational pipeline for designing, scoring, and ranking
ASO candidates against a target transcript. It integrates RNA target
accessibility analysis, thermodynamic feature computation, transcriptome-wide
off-target screening, variant/SNP awareness, and interpretable multi-factor
candidate scoring.

## Features

- **Transcript fetching** — retrieve target mRNA sequences from NCBI RefSeq
  or load from local FASTA files, with a cache layer for repeated runs.
- **Candidate generation** — sliding-window generation of ASO candidates
  (15–25 nt) across the target transcript with support for gapmer, mixmer,
  2'-OMe, MOE, and morpholino chemistries.
- **Pre-filters** — GC content, homopolymer, end stability, and ambiguous
  base filters to reject poor candidates before expensive computation.
- **RNA accessibility** — predict target RNA secondary structure and
  identify accessible binding sites using ViennaRNA (RNAfold, RNAplfold).
- **Thermodynamic scoring** — compute melting temperature (Tm) via the
  SantaLucia 1998 nearest-neighbor model, duplex ΔG, self-structure ΔG,
  target structure ΔG, and net binding energy.
- **Off-target analysis** — transcriptome-wide seed-and-extend search for
  partial-complement matches with seed-region weighting and specificity
  classification.
- **Variant awareness** — flag candidates overlapping known SNPs or
  pathogenic variants from dbSNP/ClinVar APIs or user-provided VCF files.
- **Interpretable scoring** — transparent, multi-factor candidate ranking
  with 7 score components, documented weights, and chemistry-aware bonuses.
- **Reporting** — structured JSON and self-contained HTML reports with
  summary statistics, score distribution charts, and per-candidate detail
  tables.

## Installation

```bash
pip install -e .
```

Requires **ViennaRNA** (≥ 2.4.0) for RNA structure prediction:

```bash
# macOS
brew install viennarna

# Conda
conda install -c bioconda viennarna
```

ViennaRNA is optional — the pipeline runs without it but skips accessibility
and structure-based thermodynamic stages.

## Quick Start

### Design ASO candidates for a target gene

```bash
# From NCBI (requires network)
asodesigner design --target BRCA1 --chemistry gapmer --output-dir results/

# From a local FASTA file
asodesigner design --target NM_007294 --local-fasta transcript.fasta --output-dir results/

# Keep only top 10 candidates with score ≥ 0.5
asodesigner design --target BRCA1 --top-n 10 --min-score 0.5 --output-dir results/

# Skip off-target and variant analysis for faster runs
asodesigner design --target BRCA1 --no-offtarget --no-variant --output-dir results/
```

### Score pre-generated candidates

```bash
# Score candidates from a JSON file (from save_candidates)
asodesigner score --candidates candidates.json --top-n 20

# Save scored results to a file
asodesigner score --candidates candidates.json --output scored.json
```

### Generate reports from saved results

```bash
# Generate JSON and HTML reports
asodesigner report --results scored.json --output-dir reports/

# JSON only
asodesigner report --results scored.json --output-dir reports/ --no-html
```

### Display configuration info

```bash
# Show default config
asodesigner info

# Show scoring weights and chemistry presets
asodesigner info --show-weights --show-presets
```

## Architecture

ASODesigner is organized into the following modules:

```
aso_designer/
├── __init__.py              # Package init, version info
├── __main__.py              # python -m aso_designer entry point
├── cli.py                   # CLI parser, commands, and utilities
├── pipeline.py              # ASOPipeline orchestration
├── config_loader.py         # YAML config loading and merging
├── config_schema.py         # Pydantic config validation
├── transcript_fetcher.py    # NCBI Entrez fetcher
├── local_loader.py          # Local FASTA loader
├── cache_manager.py         # Transcript cache with TTL
├── fetch_manager.py         # Unified fetch fallback chain
├── quality.py               # Transcript quality and CDS detection
├── candidate.py             # ASOCandidate dataclass
├── candidate_generator.py   # Sliding-window generation
├── chemistry.py             # Gapmer/Mixmer/Full chemistry modeling
├── filters.py               # Pre-filters (GC, homopolymer, end, ambiguous)
├── candidate_serializer.py  # JSON save/load, batch I/O, dedup
├── viennarna.py             # ViennaRNA wrappers (RNAfold, RNAplfold, etc.)
├── accessibility.py         # Target structure and accessibility scoring
├── thermodynamics.py        # NN parameters, Tm, ΔG, binding classification
├── off_target.py            # Transcriptome reference, seed-and-extend search
├── variant.py               # Variant/SNP fetch, VCF, annotation, scoring
├── scoring.py               # Multi-factor scoring engine and ranking
├── reporting.py             # JSON/HTML report generation
├── config/
│   ├── aso_config.yaml      # Default configuration
│   └── chemistry_presets.yaml  # Per-chemistry adjustments
└── utils/
    ├── logger.py            # Structured logging
    ├── exceptions.py        # ASO-specific exception hierarchy
    └── sequence.py          # Sequence validation, GC, complement, conversion
```

### Pipeline stages

1. **Fetch** — Retrieve target transcript from NCBI or local FASTA
2. **Generate** — Sliding-window candidate generation (15–25 nt)
3. **Filter** — Pre-filters: GC content, homopolymer, end stability, ambiguous bases
4. **Accessibility** — RNAplfold-based per-position accessibility scoring
5. **Thermodynamics** — Tm, duplex ΔG, self-structure ΔG, net binding energy
6. **Off-target** — Transcriptome-wide seed-and-extend search
7. **Variant** — dbSNP/ClinVar overlap annotation and penalty
8. **Score** — Multi-factor weighted composite scoring and ranking
9. **Report** — JSON and HTML report generation

Each stage can be toggled via CLI flags or configuration.

## Disclaimer

ASODesigner is an **in-silico research and candidate-ranking tool**. It does
not claim experimental efficacy, safety, or clinical utility for any
candidate produced. All candidates require experimental validation before
any downstream application.

## License

MIT
