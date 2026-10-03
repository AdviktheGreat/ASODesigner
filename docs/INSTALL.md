# Installation Guide

## Prerequisites

- Python 3.10 or higher
- pip (comes with Python)
- Git (for cloning the repository)

## Step 1: Clone the repository

```bash
git clone https://github.com/AdviktheGreat/ASODesigner.git
cd ASODesigner
```

## Step 2: Install ASODesigner

### Option A: Editable install (recommended for development)

```bash
pip install -e .
```

This installs ASODesigner in "editable" mode — changes to the source
code are immediately reflected without reinstalling.

### Option B: Regular install

```bash
pip install .
```

### Option C: Development install with test dependencies

```bash
pip install -e ".[dev]"
```

This installs ASODesigner plus development tools (pytest, ruff, mypy,
responses).

### Option D: With VCF support

```bash
pip install -e ".[vcf]"
```

This installs `cyvcf2` for fast VCF file parsing.

## Step 3: Install ViennaRNA (optional but recommended)

ViennaRNA is required for RNA structure prediction (accessibility and
structure-based thermodynamic calculations). Without it, the pipeline
runs but skips these stages.

### macOS (Homebrew)

```bash
brew install viennarna
```

### Conda

```bash
conda install -c bioconda viennarna
```

### Linux (from source)

```bash
# Download from https://www.tbi.univie.ac.at/RNA/
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/viennarna-2.6.4.tar.gz
tar xzf viennarna-2.6.4.tar.gz
cd viennarna-2.6.4
./configure --prefix=/usr/local
make
sudo make install
```

### Verify installation

```bash
# Check ViennaRNA is on PATH
which RNAfold
RNAfold --version
```

Or from Python:

```python
from aso_designer.viennarna import verify_viennarna
verify_viennarna()  # Raises ViennaRNAError if not found
```

## Step 4: Verify the installation

```bash
# Check version
asodesigner --version

# Run info command
asodesigner info

# Run the test suite
pytest
```

## Troubleshooting

### "command not found: asodesigner"

Make sure the installation directory is on your PATH. With `pip install -e .`,
the `asodesigner` command is installed in your Python's bin directory.

```bash
# Find where it was installed
pip show asodesigner | grep Location
```

### ViennaRNA tools not found

Ensure ViennaRNA binaries (RNAfold, RNAplfold, RNAduplex, RNAcofold)
are on your PATH:

```bash
echo $PATH | tr ':' '\n' | grep -i vienna
```

If not, add the ViennaRNA bin directory to your PATH:

```bash
export PATH="/usr/local/bin:$PATH"  # or wherever ViennaRNA was installed
```

### Import errors

Make sure all dependencies are installed:

```bash
pip install -r requirements.txt
```

### NCBI API rate limits

NCBI limits API requests to 3/second without an API key. To increase
this to 10/second, set the `NCBI_API_KEY` environment variable:

```bash
export NCBI_API_KEY="your_api_key_here"
```

You can obtain an API key at https://www.ncbi.nlm.nih.gov/account/.
