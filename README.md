# PlasAnnotatoR

Ensemble-based plasmid classification and annotation pipeline for metagenomes and genomes.

## Overview

PlasAnnotatoR combines four plasmid classifiers using weighted voting (AUC-based weights) and annotates detected plasmids against multiple reference databases. Results are presented in an interactive HTML report with a taxonomic network visualization.

## Requirements

### System dependencies

- micromamba (required)
- conda or Miniconda (required for PLASMe)

All other dependencies (BLAST+, DIAMOND, Prodigal, Python packages) 
are installed automatically by the installation script.

## Installation

### 1. Clone the repository

    git clone https://github.com/HaydeSaracho/PlasAnnotatoR.git
    cd PlasAnnotatoR

### 2. Create environments (~30-60 min, requires ~5 GB disk)
Creates four conda/micromamba environments:
- **plasclass**: Python 3.7 + scikit-learn 0.21.3
- **plasme**: Python 3.9 + PyTorch 1.11 (requires Miniconda)
- **plasmidhunter**: Python 3.10 + BLAST + DIAMOND + Prodigal
- **plasannotator**: Python 3.10 + all annotation tools

```bash
bash install.sh
```

### 3. Configure paths (~1 min)
Auto-detects environment paths and generates config.yaml.

    bash configure.sh

### 4. Download databases (~2-4 hours, requires ~25 GB disk)
Downloads all required databases:
- RF model (139 MB, Zenodo DOI: 10.5281/zenodo.20348779)
- PlasAnn functional databases (Zenodo DOI: 10.5281/zenodo.20501577)
- PLSDB 2025 sequences (7 GB)
- PLSDB 2025 metadata (3.4 GB)
- CARD (4.4 MB)
- MIBiG 4.0 (28 MB, Zenodo DOI: 10.5281/zenodo.20350302)
- CAZy (1.2 GB, Zenodo DOI: 10.5281/zenodo.20350742)
- PLASMe DB (12.4 GB, Zenodo)

```bash
bash download_databases.sh
```

## Usage

Activate the main environment before running:

```bash
micromamba activate plasannotator
```

Run with the included test datasets:

```bash
# Full-length sequences (genomes)
python main.py -i test_data/test_full_length.fasta -o results/test_full_length -t 8

# Fragmented contigs (metagenomes)
python main.py -i test_data/test_fragmented.fasta -o results/test_fragmented -t 8
```

Full pipeline on your data:

```bash
python main.py -i input.fasta -o results/ -t 8
```

Skip annotation (classification only, recommended for computers with less than 8 GB RAM):

```bash
python main.py -i input.fasta -o results/ --skip-annotation -t 8
```

Skip network visualization:

```bash
python main.py -i input.fasta -o results/ --skip-network -t 8
```

Adjust ensemble threshold (default: 0.6):

```bash
python main.py -i input.fasta -o results/ --threshold 0.5 -t 8
```
## Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-i` / `--input` | Input FASTA file (contigs) | required |
| `-o` / `--output` | Output directory | `results/` |
| `-t` / `--threads` | Number of threads | `8` |
| `--threshold` | Ensemble score threshold for plasmid classification | `0.6` |
| `--skip-annotation` | Skip functional annotation. Recommended for computers with less than 8 GB RAM | `False` |
| `--skip-network` | Skip taxonomic network construction | `False` |

## Output

    results/
    ├── plasannotator_report.html     Interactive HTML report
    ├── ensemble_results.tsv          Per-contig classification scores
    ├── ensemble/                     Individual tool outputs
    ├── annotation/
    │   ├── annotation_results.tsv    Full annotation table
    │   ├── card_blast.tsv            AMR hits
    │   ├── mibig_blast.tsv           BGC hits
    │   ├── cazy_diamond.tsv          CAZyme hits
    │   └── plsdb_blast.tsv           PLSDB hits + taxonomy
    └── network/
        └── plasmid_network.html      Standalone network visualization

## Author

Hayde Saracho

## Co-author

Daniel G. Kurth
