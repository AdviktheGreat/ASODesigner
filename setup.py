from setuptools import find_packages, setup

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="asodesigner",
    version="1.0.0",
    packages=find_packages(),
    package_data={
        "aso_designer": ["config/*.yaml"],
    },
    include_package_data=True,
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.81",
        "numpy>=1.24",
        "scipy>=1.10",
        "pandas>=2.0",
        "plotly>=5.14",
        "pydantic>=2.0",
        "pyyaml>=6.0",
        "requests>=2.28",
        "tqdm>=4.65",
        "click>=8.1",
    ],
    extras_require={
        "dev": [
            "pytest>=7.3",
            "pytest-cov>=4.1",
            "responses>=0.23",
            "ruff>=0.1",
            "mypy>=1.5",
        ],
        "vcf": [
            "cyvcf2>=0.30",
        ],
    },
    entry_points={
        "console_scripts": [
            "asodesigner=aso_designer.cli:main",
        ],
    },
    author="AdviktheGreat",
    author_email="advikanand7@gmail.com",
    description=(
        "In-silico antisense oligonucleotide candidate design and ranking "
        "with RNA accessibility, thermodynamics, off-target analysis, and "
        "variant awareness."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AdviktheGreat/ASODesigner",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
