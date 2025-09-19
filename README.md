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
python -m venv .venv
source .venv/bin/activate  # Linux/Mac

# Install from source
pip install -e .

# Verify the installation
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

**Get help/options**
```bash
crispr-indel-analyser -h
```

**Run the pipeline with example data**

```bash
crispr-indel-analyser \
  --fastq data/test_40k_1.merged.fastq.gz \
  --meta-csv data/sample_meta_data.csv
```

**Required options**

| Option | Description |
|--------|-------------|
| `--fastq` | Path to input FASTQ (required) |
| `--meta-csv` | Path to metadata CSV (required) |

### Output

| File Pattern | Description |
|--------------|-------------|
| `demultiplexed/*.fq.gz` | One FASTQ per sample |
| `results/*.summary.txt` | Indel counts and editing efficiency |
| `results/summary.csv` | Combined results for all samples |

Each result includes:

| Field | Description |
|-------|-------------|
| `num_ins` / `num_del` | Insertion and deletion counts |
| `num_other` | Unedited reads |
| `num_skip` |  Reads without matched barcodes |
| `per_ins` / `per_del` | Editing efficiency (%), (num_ins or num_del)/(num_ins + num_del + num_other)|
| `pos_ins` / `pos_del` | Insertion and deletion positions counts (based on target window) |

## License
MIT © 2025 Zhipeng Qu

See [LICENSE](LICENSE) for details.
