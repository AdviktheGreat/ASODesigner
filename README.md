# ASODesigner

**In-silico antisense oligonucleotide (ASO) candidate design and ranking.**

ASODesigner is a computational pipeline for designing, scoring, and ranking
ASO candidates against a target transcript. It integrates RNA target
accessibility analysis, thermodynamic feature computation, transcriptome-wide
off-target screening, variant/SNP awareness, and interpretable multi-factor
candidate scoring.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Architecture](#architecture)
- [Configuration](#configuration)
- [Scoring](#scoring)
- [Testing](#testing)
- [Programmatic Usage](#programmatic-usage)
- [Contributing](#contributing)
- [Disclaimer](#disclaimer)
- [License](#license)

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
- **Motif screening** — flags immunostimulatory motifs (CpG content,
  canonical TLR9 hexamers) and G-quadruplex-prone runs per candidate.
- **Reporting** — structured JSON and self-contained HTML reports with
  summary statistics, score distribution charts, and per-candidate detail
  tables; optional CSV and FASTA exports for candidates.

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

ViennaRNA is optional. The thermodynamics stage defaults to the fast
nearest-neighbor method, which needs no external tools. When ViennaRNA is
missing, the accessibility stage and the ViennaRNA-based thermodynamic
methods degrade gracefully: the pipeline logs a warning and continues with
neutral scores for the affected features instead of failing.

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
├── chemistry.py             # Gapmer/Mixmer/2'-OMe/MOE/morpholino modeling
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

## Configuration

ASODesigner uses a layered YAML configuration system:

1. **Bundled defaults** (`aso_designer/config/aso_config.yaml`)
2. **Chemistry presets** (`aso_designer/config/chemistry_presets.yaml`)
3. **User config** (passed via `--config`)

Later layers override earlier ones. All configuration is validated by a
Pydantic schema with cross-field validation (e.g., scoring weights must
sum to 1.0).

### Key configuration sections

```yaml
candidate:
  length: 20          # ASO length (10–30 nt)
  step: 1             # Sliding window step
  gc_min: 0.30        # Minimum GC content
  gc_max: 0.60        # Maximum GC content
  homopolymer_max: 4  # Max consecutive identical bases
  chemistry: gapmer

scoring:
  weights:
    accessibility: 0.25
    thermodynamics: 0.20
    gc_content: 0.10
    off_target: 0.20
    variant: 0.10
    self_structure: 0.10
    positional: 0.05

reporting:
  output_dir: results
  generate_html: true
  include_offtarget_details: true
```

### Chemistry presets

Each chemistry has preset adjustments for Tm, self-structure weight, and
bonuses. Available chemistries: `gapmer`, `mixmer`, `2ome`, `moe`,
`morpholino`.

```bash
# View presets for a specific chemistry
asodesigner info --chemistry mixmer --show-presets
```

### Custom configuration

```bash
# Use a custom config file
asodesigner design --target BRCA1 --config my_config.yaml --output-dir results/
```

## Scoring

Candidates are scored by a weighted composite of 7 components:

| Component | Weight | Description |
|---|---|---|
| Accessibility | 0.25 | RNAplfold-based target site accessibility |
| Thermodynamics | 0.20 | Tm score + binding class classification |
| GC Content | 0.10 | Optimal range 30–60% |
| Off-target | 0.20 | Transcriptome-wide specificity |
| Variant | 0.10 | SNP/pathogenic overlap penalty |
| Self-structure | 0.10 | Self-folding penalty |
| Positional | 0.05 | Region preference (CDS > 5UTR > 3UTR) |

Each component returns a score in [0, 1]. The composite score is the
weighted sum, plus an optional chemistry bonus, clamped to [0, 1].
Candidates are ranked by composite score and filtered by configurable
thresholds.

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=aso_designer

# Run specific test module
pytest tests/test_scoring.py -v
```

The test suite includes 810+ tests covering all modules, with mocked
ViennaRNA and NCBI calls for reproducibility.

## Programmatic Usage

ASODesigner can be used as a Python library:

```python
from aso_designer.pipeline import run_pipeline

# Run the full pipeline
result = run_pipeline(
    target="BRCA1",
    chemistry="gapmer",
    local_fasta="transcript.fasta",
    output_dir="results/",
    top_n=10,
    min_score=0.5,
)

print(f"Generated {result.total_candidates} candidates")
print(f"Ranked {result.ranked_candidates} candidates")
print(f"Top score: {result.top_score:.4f}")
print(f"Reports: {result.report_paths}")
```

### Individual module usage

```python
from aso_designer.candidate_generator import generate_candidates
from aso_designer.config_loader import load_config
from aso_designer.filters import run_all_filters
from aso_designer.scoring import score_rank_and_filter

# Load config
config = load_config(chemistry="gapmer")

# Generate candidates
candidates = generate_candidates("ATGCATGC...", config.candidate)

# Apply pre-filters
run_all_filters(candidates, gc_min=0.30, gc_max=0.60)

# Score and rank
features = [{} for _ in candidates if c.passes_filters]
ranking = score_rank_and_filter(
    [c for c in candidates if c.passes_filters],
    features,
    min_score=0.5,
    top_n=10,
)

for c in ranking.candidates:
    print(f"Rank {c.rank}: {c.sequence} score={c.final_score:.4f}")
```

### Custom scoring

```python
from aso_designer.scoring import create_score_engine, ScoreContext

# Create engine with custom weights
engine = create_score_engine(
    weights={"accessibility": 0.4, "thermodynamics": 0.3, "gc_content": 0.3},
)

# Score a single candidate
ctx = ScoreContext(
    sequence="ATGCATGCATGCATGCATGC",
    gc_content=0.5,
    accessibility_score=0.8,
    tm_score=0.9,
    binding_class="strong",
)
result = engine.score(ctx)
print(f"Final score: {result.final_score:.4f}")
```

## Contributing

### Development setup

```bash
git clone https://github.com/AdviktheGreat/ASODesigner.git
cd ASODesigner
pip install -e ".[dev]"
```

### Code style

- Linting: `ruff check aso_designer/ tests/`
- Formatting: `ruff format aso_designer/ tests/`
- Tests: `pytest`

All code must pass `ruff check`, `ruff format --check`, and `pytest`
before committing.

### Project structure

The codebase follows a 52-step implementation plan across 12 phases.
Each module is self-contained with its own test file. See
`docs/METHODOLOGY.md` for the scientific approach and `docs/API.md`
for the module API reference.

## Disclaimer

ASODesigner is an **in-silico research and candidate-ranking tool**. It does
not claim experimental efficacy, safety, or clinical utility for any
candidate produced. All candidates require experimental validation before
any downstream application.

## License

MIT
