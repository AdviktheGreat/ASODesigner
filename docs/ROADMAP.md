# ASODesigner Improvement Backlog v2 (Items 51–100)

Derived from a full codebase audit (October 2026) after improvements 1–50
were completed. Each item names the gap observed in the current code, so the
list stays verifiable if the original context is lost. Categories mirror the
audit split: core flow, analysis modules, orchestration/UX, infrastructure.

Current state at time of writing: ~22.7k lines of Python across 26 modules,
814 passing tests, ~95% coverage, CI on Python 3.10–3.12, JSON/HTML/CSV/FASTA
reporting, motif screening, disk-cached ViennaRNA results and accessibility
profiles, batched generation, k-mer off-target index.

## A. Thermodynamics & biophysics

51. **DNA/RNA hybrid NN parameters.** `thermodynamics.py:397` routes every
    chemistry to either SantaLucia DNA/DNA or Turner RNA/RNA parameters; a
    DNA-gap gapmer binding target RNA is physically an RNA/DNA hybrid. Add
    RNA/DNA hybrid parameter tables (Sugimoto 1995) and route gap-containing
    chemistries to them.
52. **Mg²⁺/divalent salt corrections.** Only the monovalent SantaLucia ΔS
    correction exists (`calculate_tm`, thermodynamics.py:490-530). Add the
    Owczarzy 2004/2008 divalent correction with `mg_conc`/`dntp_conc`
    config keys and validation.
53. **Mismatch and terminal NN terms.** The NN tables contain only matched
    pairs; add terminal-mismatch, internal-mismatch, and dangling-end terms
    so duplex ΔG near SNP sites or with overhangs is not silently over- or
    under-estimated.
54. **RNAup accessibility method.** `viennarna.py` wraps RNAfold/RNAplfold/
    RNAcofold/RNAduplex but not RNAup. Add `run_rnaup` and an
    `accessibility.method = {rnaplfold, rnaup}` option so per-site opening
    energy can complement the unpaired-probability screen.
55. **Multi-length unpaired probabilities.** RNAplfold runs with `-u 1` only
    (viennarna.py:299); run `-u 1,2,3,4` and expose per-length probabilities
    so the accessibility score can weight multi-nt unpaired stretches.
56. **Partition-function parsing.** `run_rnafold` accepts
    `partition_function=True` but never parses the `-p` output
    (viennarna.py:191-196). Parse base-pair probabilities and expose ensemble
    pairing probability as a candidate feature.
57. **ASO homodimer / self-complementarity screen.** Self-structure is
    estimated via an NN heuristic or RNAfold MFE only; add an RNAduplex-based
    ASO-vs-ASO dimer screen (3' complementarity is the classic ASO failure
    mode) as a filter or report flag.
58. **RNase H1 cleavage-site model.** RNase H eligibility is a boolean
    (`chemistry.py`). Add a positional cleavage-likelihood score for the DNA
    gap (gap length, position within gap, 5'/3' context) as a report feature
    or scoring component.

## B. Chemistry & candidate QC

59. **Strict gapmer geometry.** `GapmerDesign.build_pattern` silently mutates
    `gap_size` when length ≠ 2·wing+gap (chemistry.py:91-99). Validate the
    geometry at config load, reject impossible combinations with a clear
    error, and keep the design immutable.
60. **LNA/BNA chemistry.** Add an `lna` chemistry with per-position
    affinity-adjusted NN parameters and its own preset
    (`chemistry_presets.yaml`); extend `VALID_CHEMISTRIES`, schema enum, CLI
    choices.
61. **Phosphorothioate backbone modeling.** Chemistry is label-annotation
    only. Add a PS/PO backbone pattern (e.g. full-PS gapmer default),
    represent it in `modification_pattern`, and include it in FASTA/CSV
    output headers for ordering tools.
62. **Mixmer RNase H correctness.** `MixmerDesign.is_rnase_h_eligible`
    hardcodes False although `has_continuous_dna_gap` exists to answer
    exactly that question (chemistry.py:167-171 vs 345-371). Use it, and
    validate `mixmer.modification` against the modification enum in
    `config_schema.py`.
63. **End-stability filter defaults.** `filter_end_stability` defaults
    (`min_end_gc=0.0`, `max_end_gc=1.0`) disable the filter entirely
    (filters.py:118), while its docstring promises chemistry-specific
    thresholds. Give it real, per-chemistry defaults and wire chemistry
    through `run_all_filters`.

## C. Scoring & ranking

64. **Deterministic tie-breaking.** `rank_candidates` sorts stably by
    final_score only (scoring.py:996-1025): ties keep input order and get
    distinct ranks with no indication. Add secondary sort keys and a tie
    group in the output.
65. **Pareto/multi-objective ranking.** Only a single weighted sum exists.
    Add an optional nondominated-sort mode over efficacy (thermo/access),
    specificity (off-target), and safety (variant/motif) axes.
66. **Named weight profiles.** Add `scoring.profile = {balanced, efficacy,
    specificity, safety}` weight presets selected via config/CLI, so users
    don't hand-edit seven floats.
67. **chemistry_bonus as a schema field.** `ScoringConfig` lacks the field,
    so `create_engine_from_config` reads it through a `model_dump()` hack
    (scoring.py:898-903). Declare it (with bounds) and validate.
68. **Consume `self_structure_weight_adjust`.** The key is present in all
    five chemistry presets but never read anywhere (config_loader.py:157-180).
    Apply it to the self_structure weight, and honor the morpholino
    `rnase_h_eligible: false` preset flag.
69. **Weight-calibration utility.** Add a small tool that fits component
    weights (ridge regression) against a user-supplied CSV of
    candidates with measured activity, emitting recommended weights; keeps
    the interpretability story honest.
70. **Data-completeness transparency.** Neutral fallbacks (0.5 accessibility,
    1.0 variant, etc.) are indistinguishable from computed values in reports.
    Emit per-candidate flags for which components were computed vs neutral,
    and optionally re-rank with neutral components excluded.
71. **Consolidate duplicate scoring representations.** `ASOCandidate.scores`
    is populated nowhere while `CandidateScore` carries the real data.
    Remove the dead field or populate it consistently; one representation.

## D. Off-target & specificity

72. **Persistent transcriptome index.** The k-mer index (k=12) is rebuilt in
    memory from the FASTA on every run (off_target.py:169-195,
    pipeline.py:540-542). Serialize the index (with a cache key covering the
    FASTA hash and k) so repeat runs skip parsing/building.
73. **Parallel off-target search.** `search_batch` is a serial Python loop
    (off_target.py:438). Add multiprocessing over transcript shards with a
    worker-count config knob.
74. **Isoform-aware scoring.** `TranscriptSet.shared_regions` (quality.py:449)
    is computed but consumed by nothing. Add an isoform-conservation bonus/
    penalty component so multi-isoform data actually affects ranking.
75. **Configurable/spaced seeds.** Seed extraction is fixed at ASO[2:8) with
    a k=12 index (off_target.py). Make seed span and k configurable and add
    a spaced-seed option to recover 1-mismatch seed sites currently missed
    by construction.
76. **Top-N off-target hit details in reports.** Only aggregate
    OffTargetResult data is stored. Persist the top hits per candidate
    (transcript, position, mismatch class) and render them in JSON/HTML.

## E. Variant awareness

77. **VCF parser overhaul.** `parse_vcf` ignores CSQ/ANN consequence fields
    and FORMAT/genotypes, and multi-allelic ALT rows are misclassified as
    "mnv" (variant.py:860, :1150). Parse consequences, explode multi-allelic
    entries, and read AF from INFO/FORMAT.
78. **ClinVar MAF enrichment.** The ClinVar path never provides MAF
    (`maf=-1.0`, variant.py:770-900). Join rsIDs against dbSNP (already
    fetched) or ClinVar frequency fields so common-vs-rare weighting works
    for ClinVar-sourced variants too.
79. **Variant-fetch caching.** dbSNP/ClinVar results are refetched every run
    (variant.py:709-783). Add a TTL disk cache keyed by gene/accession like
    the transcript cache.
80. **Seed-SNP tolerance analysis.** For candidates with known common-variant
    overlaps, report which seed positions are affected and whether a
    single-position mismatch would break the seed, surfacing design
    alternatives instead of just a penalty.

## F. Transcript acquisition & annotation

81. **Annotation-driven regions.** Regions are inferred purely from ORF
    scanning (quality.py:279-345). Accept GenBank flatfiles or GFF3 to set
    real 5'UTR/CDS/3'UTR boundaries, falling back to ORF inference.
82. **Reverse-strand ORF detection.** `find_orfs` scans only the forward
    three frames (quality.py:61-105); scan both strands (or honor record
    strand metadata) so antisense-orientation transcripts aren't mislabeled.
83. **FetchManager FASTA fallback fix.** The local-FASTA fallback caches all
    records of a multi-gene FASTA under one gene id (fetch_manager.py:116-122),
    so later gene lookups can return the wrong transcript. Key per record.
84. **Unify transcript caches.** EntrezFetcher's JSON cache and
    TranscriptCacheManager overlap for the same data; also
    `resolve_gene_to_transcripts` and `fetch_gene` collide on the
    `gene:{symbol}:{organism}` cache id (transcript_fetcher.py). One cache,
    namespaced ids.
85. **Coordinate-consistency validation.** VCF CHROM is trusted blindly;
    validate chromosome/sequence consistency against the fetched transcript
    (and warn on transcript- vs genomic-coordinate confusion) before
    annotating.

## G. Pipeline, CLI, UX

86. **Stage timing instrumentation.** No per-stage durations exist anywhere
    (pipeline.py). Time each stage into `PipelineResult` and the report
    metadata; print a breakdown with a flag.
87. **Unreachable stage toggles.** `run_accessibility`/`run_thermodynamics`
    exist in config but `PipelineConfig.from_args` never maps them and the
    design CLI has no flags (pipeline.py:98-124, cli.py:102-190). Add
    `--no-accessibility`/`--no-thermodynamics`.
88. **Progress reporting and --dry-run.** Stages run silently; add a progress
    callback (tqdm is a declared dependency) and a `--dry-run` that reports
    candidate counts/estimated memory without running analysis stages.
89. **Exit-code mapping and CLI consistency.** Everything failure-shaped
    exits 1 (cli.py:306-357). Map error classes to documented codes
    (config/network/ViennaRNA). Also fix the entry-point vs prog-name
    mismatch (`asodesigner` in setup.py:33 vs `aso-designer` in cli.py:37)
    and remove the dead `binding_str` column and unused `format_summary`
    (cli.py:394-412).
90. **Results round-trip fidelity.** `score` writes a report-shaped dict, so
    `score → score` doesn't parse, and it discards any features present in
    the input (cli.py:706, cli.py:735-741). Emit `save_candidates` schema,
    carry features through, extend `validate_roundtrip` to compare features,
    and stop `default=str` silently stringifying values
    (candidate_serializer.py:198-217).
91. **`validate` subcommand.** Pre-flight checks without running the
    pipeline: config file validation, transcriptome FASTA readability/size,
    VCF parse, ViennaRNA presence/version, local FASTA validation.
92. **Persistent run logs.** The rotating file handler exists but no call
    site enables it (utils/logger.py:60-71). Write `results/{basename}.log`
    per run by default and offer `--log-format json` structured logging.

## H. Reporting & I/O

93. **Accessibility profile plots.** Embed an SVG sparkline of the unpaired
    probability across each top candidate's binding window (± flanks) in the
    HTML report; the data is already cached by the accessibility stage.
94. **Report diffing.** Add a `diff` subcommand that compares two report
    JSONs: added/removed/dropped candidates and score deltas — useful for
    config A/B testing.
95. **BED export and browser links.** Emit candidate coordinates as BED for
    genome-browser inspection and add UCSC/NCBI deep links to the HTML
    candidate table.

## I. Packaging, CI, code quality

96. **Migrate packaging to pyproject.toml.** Metadata lives in setup.py while
    pyproject.toml holds only pytest config, and requirements*.txt duplicate
    the dependency list with visible drift. Single-source everything in
    `[project]`.
97. **Type checking.** mypy is a dev dependency but unconfigured. Add
    `[tool.mypy]` (py310, relaxed start) and a CI job; annotate public APIs.
98. **Coverage gate and marker hygiene.** CI uploads coverage.xml but
    enforces nothing. Add `--cov-fail-under` (e.g. 90) with a coverage
    config, and apply the declared-but-never-used `integration` marker to
    test_integration.py so its 4 tests don't run unmarked in CI.
99. **Pre-commit hooks and Makefile.** No pre-commit config and no Makefile;
    lint/format happen only in CI. Add `.pre-commit-config.yaml` (ruff +
    ruff-format) and `make test/lint/format/coverage` targets.
100. **Dependency pruning and container story.** plotly, click, tqdm, pandas,
     scipy are declared but unused in any module (tqdm could back item 88,
     plotly could back item 93). Prune or justify, and add a Dockerfile with
     a ViennaRNA layer for reproducible environments.
