# PlasAnnotatoR

Ensemble-based plasmid classification and annotation pipeline for metagenomes and bacterial genomes.

## Overview

PlasAnnotatoR combines four plasmid classifiers using weighted voting (AUC-based weights) and annotates detected plasmids against multiple reference databases. Results are presented in an interactive HTML report with a taxonomic network visualization.

## Requirements

### System dependencies

- micromamba or conda/mamba
- NCBI BLAST+ (blastn, makeblastdb)
- DIAMOND
- Prodigal

### Tool environments

PlasAnnotatoR uses four separate environments due to dependency conflicts between tools.

| Tool | Environment | Key dependency |
|---|---|---|
| PlasClass | plasclass | Python 3.7, scikit-learn 0.21.3 |
| PLASMe | plasme | Python 3.9, PyTorch 1.11, CUDA 10.2 |
| PlasmidHunter | plasmidhunter | Python 3.10, scikit-learn 1.3.2 |
| PlasAnnotatoR | plasannotator | Python 3.10 |

## Installation

### 1. Clone the repository

    git clone https://github.com/HaydeSaracho/PlasAnnotatoR.git
    cd PlasAnnotatoR

### 2. Create environments

Run the installation script:

    bash install.sh

Or create environments manually:

    micromamba env create -f envs/plasclass_env.yml
    micromamba env create -f envs/plasmidhunter_env.yml
    micromamba env create -f envs/plasannotator_env.yml
    conda env create -f envs/plasme_env.yml

PLASMe also requires a separate installation of the tool and its database:

    conda activate plasme
    pip install git+https://github.com/HubertTang/PLASMe.git
    # Download PLASMe database following instructions at:
    # https://github.com/HubertTang/PLASMe

### 3. Configure paths

    cp config.yaml.example config.yaml

Edit config.yaml and set the paths to each environment and tool on your system:

    environments:
      plasclass:
        conda_base: /path/to/envs/plasclass
        python: /path/to/envs/plasclass/bin/python
        script: /path/to/envs/plasclass/bin/classify_fasta.py
      plasme:
        conda_base: /path/to/envs/plasme
        script: /path/to/PLASMe/PLASMe.py
        database: /path/to/PLASMe/DB
      plasmidhunter:
        conda_base: /path/to/envs/plasmidhunter
        script: /path/to/envs/plasmidhunter/bin/plasmidhunter

### 4. Download databases

Create the required directories:

    mkdir -p data/plsdb/meta data/indexes data/models

PLSDB 2025 (plasmid sequences and metadata):

    # Download from https://plsdb.github.io/plsdb/
    # Place sequences.fasta in data/plsdb/
    # Place nuccore.csv, taxonomy.csv, amr.tsv, typing.csv in data/plsdb/meta/

CARD (antimicrobial resistance genes):

    # Download from https://card.mcmaster.ca/download
    # Place nucleotide_fasta_protein_homolog_model.fasta in data/indexes/

MIBiG 4.0 (biosynthetic gene clusters):

    # Download from https://mibig.secondarymetabolites.org/download
    # Place mibig_4.0.fasta in data/indexes/

CAZy (carbohydrate-active enzymes):

    # Download from http://www.cazy.org/ or https://bcb.unl.edu/dbCAN2/download/
    # Place cazy.fasta in data/indexes/

RF model:

    # Download from the PlasAnnotatoR releases page on GitHub
    # Place rf_model.pkl in data/models/

## Usage

Activate the main environment before running:

    micromamba activate plasannotator

Full pipeline:

    python main.py -i input.fasta -o results/ -t 8

Skip annotation (classification only, much faster):

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
