# ASODesigner API Reference

Module-by-module API reference for ASODesigner.

## Table of Contents

- [Core Modules](#core-modules)
  - [aso_designer.candidate](#aso_designercandidate)
  - [aso_designer.chemistry](#aso_designerchemistry)
  - [aso_designer.filters](#aso_designerfilters)
  - [aso_designer.candidate_serializer](#aso_designercandidate_serializer)
- [Analysis Modules](#analysis-modules)
  - [aso_designer.viennarna](#aso_designerviennarna)
  - [aso_designer.accessibility](#aso_designeraccessibility)
  - [aso_designer.thermodynamics](#aso_designerthermodynamics)
- [Off-target & Variant](#off-target--variant)
  - [aso_designer.off_target](#aso_designeroff_target)
  - [aso_designer.variant](#aso_designervariant)
- [Scoring & Reporting](#scoring--reporting)
  - [aso_designer.scoring](#aso_designerscoring)
  - [aso_designer.reporting](#aso_designerreporting)
- [Pipeline & CLI](#pipeline--cli)
  - [aso_designer.pipeline](#aso_designerpipeline)
  - [aso_designer.cli](#aso_designercli)
- [Configuration](#configuration)
  - [aso_designer.config_loader](#aso_designerconfig_loader)
  - [aso_designer.config_schema](#aso_designerconfig_schema)
- [Utilities](#utilities)
  - [aso_designer.utils.sequence](#aso_designerutilssequence)
  - [aso_designer.utils.exceptions](#aso_designerutilsexceptions)

---

## Core Modules

### aso_designer.candidate

The `ASOCandidate` dataclass carries a candidate through the pipeline.

```python
from aso_designer.candidate import ASOCandidate

cand = ASOCandidate(
    sequence="ATGCATGCATGCATGCATGC",
    start=100, end=120, length=20,
    chemistry="gapmer",
    gc_content=0.5,
    region="CDS",
)
```

**Key attributes:** `sequence`, `start`, `end`, `length`, `chemistry`,
`gc_content`, `region`, `passes_filters`, `features`, `scores`, `rank`.

**Methods:** `to_dict()`, `from_dict(data)`, `__repr__()`.

### aso_designer.chemistry

Chemistry modeling for gapmer, mixmer, 2'-OMe, MOE, and morpholino
modification designs.

```python
from aso_designer.chemistry import create_chemistry_design

design = create_chemistry_design("gapmer", length=20, dna_window=10)
print(design.modification_pattern)  # ["MOE", "MOE", ..., "DNA", ..., "MOE", "MOE"]
```

**Classes:** `GapmerDesign`, `MixmerDesign`, `FullChemistryDesign`.

**Factory:** `create_chemistry_design(chemistry, length, **kwargs)`.

### aso_designer.filters

Pre-filters for candidate screening before feature computation.

```python
from aso_designer.filters import run_all_filters

result = run_all_filters(candidates, gc_min=0.30, gc_max=0.60)
print(result["total_passed"], "/", len(candidates))
```

**Functions:** `filter_gc_content()`, `filter_homopolymer()`,
`filter_end_stability()`, `filter_ambiguous_bases()`,
`run_all_filters()`.

### aso_designer.candidate_serializer

JSON serialization for candidate lists.

```python
from aso_designer.candidate_serializer import save_candidates, load_candidates

save_candidates(candidates, "candidates.json", metadata={"gene": "BRCA1"})
loaded = load_candidates("candidates.json")
```

**Functions:** `save_candidates()`, `load_candidates()`,
`save_candidates_batch()`, `load_candidates_batch()`,
`deduplicate_candidates()`, `validate_round_trip()`.

---

## Analysis Modules

### aso_designer.viennarna

ViennaRNA command-line tool wrappers.

```python
from aso_designer.viennarna import run_rnafold, run_rnaplfold, run_rnaduplex

# Predict MFE structure
result = run_rnafold("AUGCAUGC...")
print(result.mfe, result.structure)

# Local accessibility
result = run_rnaplfold("AUGCAUGC...", window_size=70)

# Duplex energy
result = run_rnaduplex("ASO_sequence", "target_sequence")
print(result.dg, result.start_5, result.end_5)
```

**Functions:** `verify_viennarna()`, `run_rnafold()`, `run_rnaplfold()`,
`run_rnacofold()`, `run_rnaduplex()`.

**Result classes:** `RNAfoldResult`, `RNAplfoldResult`, `RNAcofoldResult`,
`RNAduplexResult`.

### aso_designer.accessibility

Target RNA structure prediction and accessibility scoring.

```python
from aso_designer.accessibility import TargetStructurePredictor

predictor = TargetStructurePredictor()
structure = predictor.predict_structure("AUGCAUGC...", use_cache=True)
accessibility = predictor.compute_accessibility_profile("AUGCAUGC...")
score = predictor.score_window(accessibility, start=10, end=30)
```

**Classes:** `TargetStructurePredictor`, `AccessibilityProfile`,
`TargetStructure`.

### aso_designer.thermodynamics

Thermodynamic calculations using NN model and ViennaRNA.

```python
from aso_designer.thermodynamics import calculate_tm, calculate_duplex_dg

tm = calculate_tm("ATGCATGCATGCATGCATGC", chemistry="gapmer")
print(f"Tm = {tm:.1f} °C")

dg = calculate_duplex_dg("ATGCATGC", "AUGCAUGC", method="nearest_neighbor")
print(f"ΔG = {dg.dg:.2f} kcal/mol")
```

**Functions:** `calculate_tm()`, `calculate_duplex_dg()`,
`calculate_self_structure_dg()`, `calculate_target_structure_dg()`,
`compute_thermodynamic_summary()`.

**Classes:** `ThermodynamicSummary`, `DuplexResult`, `SelfStructureResult`.

---

## Off-target & Variant

### aso_designer.off_target

Transcriptome-wide off-target search.

```python
from aso_designer.off_target import TranscriptomeReference, OffTargetSearchEngine

ref = TranscriptomeReference(kmer_size=12)
ref.add_transcript("NM_001", "GENE1", "ATGCATGC...")

engine = OffTargetSearchEngine(ref)
hits = engine.search("ATGCATGCATGCATGCATGC", max_mismatches=3)
result = classify_and_filter_hits(hits, exclude_on_target=True)
```

**Classes:** `TranscriptomeReference`, `OffTargetSearchEngine`,
`OffTargetHit`, `OffTargetResult`.

**Functions:** `calculate_off_target_score()`, `classify_specificity()`,
`score_off_target()`, `score_off_target_batch()`.

### aso_designer.variant

Variant/SNP awareness with dbSNP, ClinVar, and VCF support.

```python
from aso_designer.variant import VariantFetchManager, VariantAnnotator

fm = VariantFetchManager()
variants = fm.fetch("BRCA1", "1")

annotator = VariantAnnotator(seed_size=8)
annotation = annotator.annotate(candidate, variants)
print(annotation.penalty, annotation.has_pathogenic)
```

**Classes:** `Variant`, `VariantSet`, `DbSNPClient`, `ClinVarClient`,
`VariantFetchManager`, `VariantAnnotator`.

**Functions:** `load_vcf()`, `load_vcf_filtered()`,
`score_variant()`, `classify_variant_risk()`,
`filter_candidates_by_variant_risk()`.

---

## Scoring & Reporting

### aso_designer.scoring

Multi-factor scoring engine and ranking.

```python
from aso_designer.scoring import create_score_engine, score_rank_and_filter

engine = create_score_engine()
# Score with 7 components: accessibility, thermodynamics, gc_content,
# off_target, variant, self_structure, positional

result = score_rank_and_filter(
    candidates, features_list, min_score=0.5, top_n=10
)
for c in result.candidates:
    print(f"Rank {c.rank}: {c.sequence} score={c.final_score:.4f}")
```

**Classes:** `ScoreContext`, `ScoreResult`, `CompositeScore`,
`CandidateScore`, `RankingResult`, `ScoreEngine`, and 7 score
component classes.

**Functions:** `create_score_engine()`, `create_engine_from_config()`,
`score_candidate()`, `score_candidates_batch()`, `rank_candidates()`,
`filter_by_score()`, `rank_and_filter()`, `score_rank_and_filter()`,
`build_features_dict()`.

### aso_designer.reporting

JSON and HTML report generation.

```python
from aso_designer.reporting import assemble_report

paths = assemble_report(
    scored_candidates=ranking.candidates,
    output_dir="reports/",
    target_name="BRCA1",
    target_id="NM_007294",
    chemistry="gapmer",
    basename="brca1_report",
)
print(paths)  # {"json": Path, "html": Path}
```

**Classes:** `ReportData`, `ReportSummary`, `ReportAssembler`.

**Functions:** `generate_json_report()`, `load_json_report()`,
`generate_html_report()`, `generate_report()`, `build_report_data()`,
`assemble_report()`, `assemble_report_from_config()`,
`generate_report_from_ranking()`, `generate_reports_batch()`.

---

## Pipeline & CLI

### aso_designer.pipeline

Pipeline orchestration.

```python
from aso_designer.pipeline import run_pipeline, ASOPipeline, PipelineConfig

# Convenience function
result = run_pipeline(
    target="BRCA1",
    chemistry="gapmer",
    local_fasta="transcript.fasta",
    output_dir="results/",
    top_n=10,
)

# Or use the class directly
config = PipelineConfig(target="BRCA1", chemistry="gapmer")
pipeline = ASOPipeline(config)
result = pipeline.run()
```

**Classes:** `PipelineConfig`, `PipelineResult`, `ASOPipeline`.

**Functions:** `run_pipeline()`.

### aso_designer.cli

Command-line interface.

```bash
# Full pipeline
asodesigner design --target BRCA1 --chemistry gapmer --output-dir results/

# Score pre-generated candidates
asodesigner score --candidates candidates.json --top-n 20

# Generate reports from saved results
asodesigner report --results scored.json --output-dir reports/

# Show config info
asodesigner info --show-weights --show-presets
```

**Entry point:** `main()` (also via `python -m aso_designer`).

---

## Configuration

### aso_designer.config_loader

YAML configuration loading with layered merge.

```python
from aso_designer.config_loader import load_config

config = load_config("my_config.yaml", chemistry="gapmer")
print(config.candidate.length)  # 20
print(config.scoring.weights.accessibility)  # 0.25
```

**Functions:** `load_config()`, `load_yaml()`.

### aso_designer.config_schema

Pydantic models for configuration validation.

**Classes:** `ASOConfig`, `CandidateConfig`, `GapmerConfig`,
`MixmerConfig`, `ScoringConfig`, `ScoringWeightsConfig`,
`ReportingConfig`, `AccessibilityConfig`, `ThermodynamicsConfig`,
`OffTargetConfig`, `VariantConfig`.

---

## Utilities

### aso_designer.utils.sequence

Sequence utilities.

```python
from aso_designer.utils.sequence import (
    validate_sequence, gc_content, reverse_complement,
    has_homopolymer, dna_to_rna, rna_to_dna,
)

gc = gc_content("ATGCATGC")  # 0.5
rna = dna_to_rna("ATGC")  # "AUGC"
```

**Functions:** `validate_sequence()`, `gc_content()`,
`reverse_complement()`, `has_homopolymer()`, `dna_to_rna()`,
`rna_to_dna()`, `complement()`.

### aso_designer.utils.exceptions

ASO-specific exception hierarchy.

```
ASOError (base)
├── ConfigurationError
├── TranscriptFetchError
├── CandidateGenerationError
├── ChemistryError
├── FilterError
├── ViennaRNAError
├── AccessibilityError
├── ThermodynamicsError
├── OffTargetError
├── VariantError
└── ScoringError
```
