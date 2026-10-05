# Changelog

All notable changes to ASODesigner are documented in this file.

## [Unreleased]

### Fixed

- Removed the invalid `full` chemistry from the CLI choices and docs;
  it was never supported by the chemistry factory.
- Corrected the RNA nearest-neighbor parameters to the NNDB/Xia 1998
  values and rebuilt RNA duplex initiation (intermolecular initiation,
  terminal AU penalty, symmetry correction).
- Fixed an off-by-one in ORF detection that skipped one codon per
  recorded ORF.
- Multi-allelic dbSNP MAF values (comma/pipe separated lists) now
  parse correctly, taking the maximum; extra colon segments and
  unparseable values are handled.

### Changed

- Pipeline stage errors now catch only expected failure classes
  (ASODesignError, RuntimeError, OSError, ValueError); programming
  errors propagate. Thermodynamics failures degrade non-fatally like
  accessibility and off-target.
- Narrowed broad `except Exception` blocks in candidate generation,
  accessibility caching, config validation, chemistry presets, and
  NCBI API clients to the specific exception types they can raise.
- Scoring components are registered via direct class references
  instead of importlib lookups.
- README and API docs updated: real test counts, accurate
  ViennaRNA-optional description, chemistry list.

### Added

- Real pipeline wiring for accessibility, thermodynamics, off-target,
  and variant stages; `--transcriptome` and `--vcf` CLI flags;
  configurable thermodynamics method and NN self-structure fraction.
- ViennaRNA result disk cache (opt-in via `thermodynamics.cache_dir`),
  chemistry preset caching, and `lru_cache` memoization of the
  nearest-neighbor thermodynamic core.
- K-mer set-based shared-region detection across transcript isoforms
  (O(L·n) vs the previous O(L²·n·m) scan).
- Memory-managed batched candidate generation in the pipeline stage.
- Off-target verification now pre-filters candidates by seed-region
  mismatches and uses early-terminating Hamming distance.

### Tests

- Grew from 602 to 687 tests, including pipeline error-contract,
  cache-invalidation, batched-generation equivalence, and NCBI
  failure-path coverage.

## [1.0.0] — 2026-10-04

Initial release of ASODesigner. All 52 implementation steps complete
across 12 phases. 602 tests passing, 85% code coverage.

### Highlights

- Full ASO design pipeline: fetch → generate → filter → score → report
- 7 score components with transparent weighted composite scoring
- 5 chemistry types: gapmer, mixmer, 2'-OMe, MOE, morpholino
- JSON and HTML report generation with SVG score distribution charts
- CLI with 4 subcommands: design, score, report, info
- Configurable via YAML with chemistry presets and Pydantic validation
- 602 tests with 85% code coverage
- Comprehensive documentation: README, methodology, API reference,
  installation guide, usage guide, and examples

### Statistics

- 28 Python modules in `aso_designer/`
- 18 test files with 602 tests
- 54 files passing ruff check and format
- 85% test coverage (excluding network-dependent and CLI subprocess code)
- 6 example scripts and configurations

## [0.0.1] — 2026-10-03

Initial development release. All core pipeline stages implemented and
tested (550+ tests passing).

### Phase 1: Package skeleton (Steps 1-4)

- Package directory structure with `aso_designer/` and `tests/`
- Utility infrastructure: structured logger, ASO-specific exception
  hierarchy, sequence utilities (validation, complement, GC content,
  homopolymer, DNA↔RNA conversion)
- Configuration system: YAML configs with inline documentation,
  chemistry presets, Pydantic schema with cross-field validation,
  config loader with layered merge
- Development tooling: pyproject.toml, ruff.toml, CI workflow,
  requirements files

### Phase 2: Transcript fetching (Steps 5-8)

- NCBI Entrez fetcher (ESearch + EFetch)
- Local FASTA loader (Biopython)
- Cache manager with TTL
- FetchManager with NCBI → cache → local file fallback chain
- Transcript quality validation, CDS detection, multi-isoform
  TranscriptSet, batch quality report

### Phase 3: Candidate generation (Steps 9-13)

- ASOCandidate dataclass with full metadata
- Sliding-window candidate generation with exon-aware coordinate
  mapping and memory-managed batching
- Chemistry modeling: GapmerDesign, MixmerDesign, FullChemistryDesign,
  factory function, mixmer scheme parsing
- Pre-filters: GC content, homopolymer, end stability, ambiguous bases
- JSON serialization with round-trip validation, batch directory I/O,
  deduplication

### Phase 4: ViennaRNA & accessibility (Steps 14-17)

- ViennaRNA wrappers: RNAfold, RNAplfold, RNAcofold, RNAduplex with
  multi-strategy output parsing
- TargetStructurePredictor with disk caching
- Per-position and windowed accessibility scoring
- Accessibility classification and profile export (JSON/CSV)

### Phase 5: Thermodynamics (Steps 18-22)

- SantaLucia 1998 DNA NN parameters (all 16 dinucleotides)
- Turner 2004 RNA NN parameters
- Tm calculation with SantaLucia salt correction (ΔS adjustment)
- Chemistry-aware Tm with preset adjustments
- Duplex ΔG (RNAduplex + NN model)
- Self-structure ΔG (RNAfold + NN) with chemistry-aware penalty
- Target structure ΔG, net binding energy, binding classification
- ThermodynamicSummary dataclass

### Phase 6: Off-target analysis (Steps 23-27)

- TranscriptomeReference with k-mer index
- OffTargetSearchEngine with seed-and-extend strategy
- Hit classification and filtering
- OffTargetResult summary
- Weighted penalty scoring, specificity classification, batch scoring

### Phase 7: Variant/SNP awareness (Steps 28-31)

- Variant and VariantSet data structures with position-based indexing
- DbSNP API client (ESearch/ESummary)
- ClinVar API client
- VariantFetchManager with API → VCF fallback
- VCF parser with INFO field parsing
- VariantAnnotator with overlap detection and seed/non-seed classification
- MAF-weighted penalty, pathogenic flagging, risk classification

### Phase 8: Scoring & ranking (Steps 32-36)

- ScoreContext, ScoreResult, ScoreComponent protocol
- ScoreEngine with weighted composite scoring
- 7 score components: Accessibility, Thermodynamics, GCContent,
  OffTarget, Variant, SelfStructure, Positional
- Config-driven engine creation
- Convenience functions: score_candidate, score_candidates_batch,
  build_features_dict, rank_candidates, filter_by_score,
  rank_and_filter, score_rank_and_filter
- RankingResult dataclass

### Phase 9: Reporting (Steps 37-40)

- ReportData and ReportSummary dataclasses with round-trip serialization
- JSON report generation and loading
- HTML report with embedded CSS, summary cards, candidate table,
  SVG score distribution chart, and detail sections
- ReportAssembler with fluent interface
- Config-driven and batch report generation

### Phase 10: CLI & integration (Steps 41-44)

- CLI module with argparse, four subcommands (design, score, report, info)
- `__main__.py` entry point for `python -m aso_designer`
- Argument validation helpers (validate_target, validate_input_file,
  is_accession)
- CLI display helpers (progress, summary box, results table)
- PipelineConfig and PipelineResult dataclasses
- ASOPipeline class with stage dispatch and graceful degradation
- run_pipeline convenience function
- Full CLI command implementation for all four subcommands
- 47 CLI and integration tests

### Phase 11: Documentation (Steps 45-48)

- README with features, installation, quick start, architecture,
  configuration, scoring, testing, programmatic usage, and contributing
- Methodology documentation (NN model, thermodynamics, off-target,
  variant, scoring)
- API reference for all modules
- Installation guide with ViennaRNA setup and troubleshooting
- CLI usage guide with all commands and options
- Example scripts and custom config
