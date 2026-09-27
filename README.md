# ASODesigner

**In-silico antisense oligonucleotide (ASO) candidate design and ranking.**

ASODesigner is a computational pipeline for designing, scoring, and ranking
ASO candidates against a target transcript. It integrates RNA target
accessibility analysis, thermodynamic feature computation, transcriptome-wide
off-target screening, variant/SNP awareness, and interpretable multi-factor
candidate scoring.

## Features

- **Transcript fetching** — retrieve target mRNA sequences from NCBI RefSeq
  or load from local FASTA files.
- **Candidate generation** — generate ASO candidates across the target
  transcript with support for gapmer, mixmer, 2'-OMe, MOE, and morpholino
  chemistries.
- **RNA accessibility** — predict target RNA secondary structure and identify
  accessible binding sites using ViennaRNA.
- **Thermodynamic scoring** — compute melting temperature (Tm), duplex free
  energy, self-structure energy, and net binding energy.
- **Off-target analysis** — transcriptome-wide search for partial-complement
  matches with seed-region weighting.
- **Variant awareness** — flag candidates overlapping known SNPs or
  pathogenic variants from dbSNP/ClinVar or user-provided VCF files.
- **Interpretable scoring** — transparent, multi-factor candidate ranking
  with documented score components and weights.

## Installation

```bash
pip install -e .
```

Requires **ViennaRNA** (≥ 2.4.0) for RNA structure prediction:

```bash
# macOS
brew install vienna-rna

# Conda
conda install -c bioconda viennarna
```

See `docs/INSTALL.md` for detailed instructions.

## Quick Start

```bash
# Design ASO candidates for a target gene
asodesigner design --gene DMD --organism "Homo sapiens" --output results/

# Design from a local FASTA file
asodesigner design --input transcript.fasta --output results/
```

See `docs/USAGE.md` for full command reference.

## Disclaimer

ASODesigner is an **in-silico research and candidate-ranking tool**. It does
not claim experimental efficacy, safety, or clinical utility for any
candidate produced. All candidates require experimental validation before
any downstream application.

## License

MIT
