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

### 2. Create environments

    bash install.sh

### 3. Configure paths

    bash configure.sh

### 4. Download databases

    bash download_databases.sh

## Usage

Activate the main environment before running:

    micromamba activate plasannotator

Full pipeline:

    python main.py -i input.fasta -o results/ -t 8

Skip annotation (classification only, much faster, recommended for computers with less than 8 GB RAM):

    python main.py -i input.fasta -o results/ --skip-annotation -t 8

Skip network visualization:

    python main.py -i input.fasta -o results/ --skip-network -t 8

## Arguments

| Argument | Description | Default |
|---|---|---|
| -i / --input | Input FASTA file (contigs) | required |
| -o / --output | Output directory | results/ |
| -t / --threads | Number of threads | 8 |
| -c / --config | Config file | config.yaml |
| --skip-annotation | Skip Layer 2 | False |
| --skip-network | Skip Layer 3 | False |

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

## Pipeline layers

    Layer 1 - Ensemble classifier
        PlasClass (AUC 0.9017) + PLASMe (AUC 0.9748) +
        PlasmidHunter (AUC 0.9178) + custom RF model (AUC 0.9872)
        weighted vote -> plasmid/chromosome label + ensemble score

    Layer 2 - Annotation
        BLAST vs CARD (AMR genes)
        BLAST vs MIBiG 4.0 (biosynthetic gene clusters)
        DIAMOND blastp vs CAZy (carbohydrate-active enzymes)
        BLAST vs PLSDB 2025 (identity + taxonomy)

    Layer 3 - Contextual network
        Interactive pyvis network: Genus -> Plasmid
        Edge length proportional to PLSDB identity

    Layer 4 - HTML report
        Statistics, annotated plasmid table, embedded network

## Custom RF model

Trained on PLSDB 2025 (72,556 plasmids) and chromosomal fragments from 3,362 RefSeq Complete genomes. Features: 5-mer frequencies (1,024 features). Fragments: 500 bp to 500 kbp.

| Tool | AUC |
|---|---|
| RF model (custom) | 0.9872 |
| PLASMe | 0.9748 |
| PlasmidHunter | 0.9178 |
| PlasClass | 0.9017 |

## Author

Hayde Saracho
