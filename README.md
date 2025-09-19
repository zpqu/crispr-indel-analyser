# CRISPR Indel Analyser

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python: 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)

CRISPR Indel Aanalyser is a lightweight, modular Python pipeline for **demultiplexing and indel frequency analysis** of CRISPR sequencing data. Supports barcode-based sample separation and accurate indel detection.

Ideal for analyzing base editing, prime editing, or knockout efficiency from amplicon-seq data.

---

## Installation

```bash
# Clone the repo
git clone https://github.com/zpqu/crispr-indel-analyser.git
cd crispr-indel-analyser

# (Optional) Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# Install from source
pip install -e .

# Verigy the installation
crispr-indel-analyser
```

## Quick Start

### 1. Prepare Input files

**meta.csv** -- Sample meta file

```csv
sample,target,up,down,fp,rp
sample1,AAACCCGGG,TT,AA,AGGTCA,CTAGCT
sample2,TGCAGCTA,CA,TG,GATCCA,CGAAGT
```

| Field | Description |
|-------|-------------|
| `sample` | Sample name |
| `target` | Target editing sequence |
| `up` / `down` | Flanking sequences (used to locate window) |
| `fp` / `rp` | Forward and reverse PCR primers (used as barcodes) |

**mixed.fq.gz** -- FASTQ file with barcoded samples

### 2. Run the Pipeline

```bash
crispr-indel-analyser \
  --fastq data/mixed.fq.gz \
  --meta-csv data/meta.csv \
  --demux-dir demultiplexed \
  --result-dir results \
  --mismatch 1
```

**Options**

| Option | Description |
|--------|-------------|
| `--fastq` | Path to input FASTQ (required) |
| `--meta-csv` | Path to metadata CSV (required) |
| `--demux-dir` | Output directory for split FASTQs (default: `demultiplexed`) |
| `--result-dir` | Output directory for analysis results (default: `results`) |
| `--mismatch` | Allowed barcode mismatches: 0, 1, or 2 (default: `0`) |

### Output

| File Pattern | Description |
|--------------|-------------|
| `demultiplexed/*.fq.gz` | One FASTQ per sample |
| `results/*.summary.txt` | Indel counts and editing efficiency |
| `results/summary.csv` | Combined results for all samples |

Each result includes:

| Field | Description |
|-------|-------------|
| `notINDEL` | Unedited reads |
| `NumIns` / `NumDels` | Insertion and deletion counts |
| `Ins_per` / `Dels_per` | Editing efficiency (%) |

## License
MIT © 2025 Your Name
See [LICENSE](LICENSE) for details.
