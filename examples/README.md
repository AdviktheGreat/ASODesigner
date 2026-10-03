# ASODesigner Examples

This directory contains example scripts and configurations for
ASODesigner.

## Files

| File | Description |
|---|---|
| `basic_design.py` | Simplest pipeline invocation — generate, filter, score, report |
| `custom_config.yaml` | Custom YAML configuration with modified weights and parameters |
| `custom_scoring.py` | Custom scoring engine with non-default weights and chemistry bonus |
| `batch_design.py` | Multi-target batch processing with per-target reports |
| `module_usage.py` | Individual module usage for custom workflows |

## Running examples

All examples are self-contained and can be run directly:

```bash
# Basic pipeline
python examples/basic_design.py

# Custom scoring
python examples/custom_scoring.py

# Batch design
python examples/batch_design.py

# Individual module usage
python examples/module_usage.py
```

The examples create temporary files and directories, so they don't
leave any artifacts on your system.

## Using the custom config

```bash
# Use the custom config with the CLI
asodesigner design --target BRCA1 --config examples/custom_config.yaml --output-dir results/
```

## Prerequisites

- ASODesigner installed (`pip install -e .`)
- The basic, batch, and module examples work without ViennaRNA
  (they skip accessibility and thermodynamics stages)
- For full pipeline runs with all stages, install ViennaRNA
  (see `docs/INSTALL.md`)
