# Usage Guide

## Commands

ASODesigner provides five CLI commands:

| Command | Description |
|---|---|
| `design` | Full pipeline: fetch → generate → filter → score → report |
| `score` | Score pre-generated candidates from a JSON file |
| `report` | Generate JSON/HTML reports from saved results |
| `info` | Display configuration, weights, and presets |
| `calibrate` | Recommend scoring weights from measured activity data (CSV) |

## design

Run the complete ASO design pipeline.

### Basic usage

```bash
asodesigner design --target BRCA1 --output-dir results/
```

### From a local FASTA file

```bash
asodesigner design --target NM_007294 --local-fasta transcript.fasta --output-dir results/
```

### All options

```
asodesigner design --target <TARGET> [options]

Required:
  --target TARGET          Gene symbol (BRCA1) or accession (NM_007294)

Common:
  --config PATH            YAML config file (overrides defaults)
  --chemistry CHEM         ASO chemistry: gapmer, mixmer, 2ome, moe,
                           morpholino, lna
  --profile PROFILE        Named scoring weight profile: balanced,
                           efficacy, specificity, safety
  --organism ORG           Organism for NCBI lookup (default: Homo sapiens)
  --local-fasta PATH       Local FASTA file (bypasses NCBI)
  --transcriptome PATH     Transcriptome FASTA for off-target search
                           (enables off-target analysis; overrides
                           off_target.transcriptome_path)
  --vcf PATH               VCF file with known variants in transcript
                           coordinates (enables variant analysis;
                           overrides variant.vcf_path)
  --output-dir DIR         Output directory for reports (default: results)
  --top-n N                Keep only the top N candidates
  --min-score SCORE        Minimum score threshold (default: 0.0)
  --csv                    Also export ranked candidates as CSV
                           ({basename}_candidates.csv)
  --export-fasta           Also export candidate sequences as FASTA
                           ({basename}_candidates.fasta)
  --no-html                Skip HTML report generation
  --no-offtarget           Skip off-target analysis
  --no-variant             Skip variant/SNP annotation
  --verbose                Enable DEBUG-level logging
  --quiet                  Suppress all output except errors
```

### Examples

```bash
# Gapmer design for BRCA1, top 20 candidates
asodesigner design --target BRCA1 --chemistry gapmer --top-n 20 --output-dir results/

# Mixmer design from local FASTA, min score 0.5
asodesigner design --target NM_007294 --chemistry mixmer --local-fasta transcripts.fasta --min-score 0.5

# Full analysis with off-target search and variant annotation
asodesigner design --target BRCA1 \
  --transcriptome refseq_mrna_human.fna \
  --vcf brca1_variants.vcf \
  --output-dir results/

# Fast run without off-target or variant analysis
asodesigner design --target BRCA1 --no-offtarget --no-variant --output-dir results/

# Custom config file
asodesigner design --target BRCA1 --config my_config.yaml --output-dir results/

# 2'-OMe chemistry with verbose logging
asodesigner design --target DMD --chemistry 2ome --verbose --output-dir results/
```

## score

Score pre-generated candidates from a JSON file (created by
`save_candidates` or the `candidate_serializer` module).

### Basic usage

```bash
# Score and print to stdout
asodesigner score --candidates candidates.json --top-n 20

# Save scored results to a file
asodesigner score --candidates candidates.json --output scored.json
```

### All options

```
asodesigner score --candidates PATH [options]

Required:
  --candidates PATH        JSON file of candidates

Common:
  --output PATH            Output JSON file (default: stdout)
  --top-n N                Keep only the top N candidates
  --min-score SCORE        Minimum score threshold (default: 0.0)
  --config PATH            YAML config file
  --chemistry CHEM         ASO chemistry
  --verbose                Enable DEBUG-level logging
  --quiet                  Suppress all output except errors
```

## report

Generate JSON and HTML reports from a saved results file.

### Basic usage

```bash
asodesigner report --results scored.json --output-dir reports/
```

### All options

```
asodesigner report --results PATH [options]

Required:
  --results PATH           JSON results file

Common:
  --output-dir DIR         Output directory (default: reports)
  --basename NAME          Base filename (default: aso_report)
  --csv                    Also export ranked candidates as CSV
                           ({basename}_candidates.csv)
  --export-fasta           Also export candidate sequences as FASTA
                           ({basename}_candidates.fasta)
  --no-html                Skip HTML report generation
  --config PATH            YAML config file
  --chemistry CHEM         ASO chemistry
  --verbose                Enable DEBUG-level logging
  --quiet                  Suppress all output except errors
```

## info

Display ASODesigner configuration, scoring weights, and chemistry presets.

### Basic usage

```bash
# Show default config
asodesigner info

# Show scoring weights
asodesigner info --show-weights

# Show chemistry presets
asodesigner info --show-presets

# Show everything for a specific chemistry
asodesigner info --chemistry mixmer --show-weights --show-presets
```

## calibrate

Recommend scoring weights from measured ASO activity data.

### Basic usage

```bash
asodesigner calibrate --csv results.csv
```

The CSV needs a header row with the seven component-score columns
(`accessibility`, `thermodynamics`, `gc_content`, `off_target`,
`variant`, `self_structure`, `positional`), a numeric `activity` column
(higher = better, e.g. percent knockdown), and at least nine data rows.
Extra columns are ignored.

### Options

```
asodesigner calibrate --csv PATH [--ridge F] [--out PATH]

  --csv PATH        CSV with component scores and measured activity
  --ridge F         Ridge penalty (default 1.0; increase when component
                    scores are collinear)
  --out PATH        Optionally write a scoring weights YAML snippet
```

The command prints the recommended non-negative weights (summing to
1.0), the fit quality (R²), and a warning when R² < 0.2. Pass `--out`
to write a `scoring.weights` snippet you can paste into a config YAML.

### Examples

```bash
# Fit weights and save them for future design runs
asodesigner calibrate --csv screen_2026.csv --ridge 2.0 --out weights.yaml
asodesigner design --target BRCA1 --config config_with_weights.yaml
```

## Output files

The `design` and `report` commands always generate two files in the
output directory:

- `{basename}.json` — Machine-readable JSON report with summary,
  candidates (including motif-screen features), and optional detail
  sections
- `{basename}.html` — Self-contained HTML report with embedded CSS,
  summary cards, score distribution chart, and candidate table

Optional exports (enabled with `--csv` / `--export-fasta`):

- `{basename}_candidates.csv` — Ranked candidates as a table (core
  columns, per-component scores, and motif-screen columns)
- `{basename}_candidates.fasta` — Candidate sequences as FASTA records
  (for oligo ordering / downstream tools)

The `basename` defaults to `{target}_report` for the `design` command
and `aso_report` for the `report` command.

## Environment variables

| Variable | Description |
|---|---|
| `NCBI_API_KEY` | NCBI API key for higher rate limits (10/s vs 3/s) |
| `NCBI_EMAIL` | Email for NCBI identification |

## Programmatic usage

See the [README](../README.md#programmatic-usage) for Python API examples.
