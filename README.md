# PlasAnnotatoR

Ensemble-based plasmid classification and annotation pipeline for metagenomes and bacterial genomes.

## Overview

PlasAnnotatoR combines four plasmid classifiers using weighted voting (AUC-based weights) and annotates detected plasmids against multiple reference databases. Results are presented in an interactive HTML report with a taxonomic network visualization.

## Requirements

### Conda environments

| Tool | Environment |
|---|---|
| PlasClass | micromamba/envs/plasclass (Python 3.7, scikit-learn 0.21.3) |
| PLASMe | anaconda3/envs/plasme |
| PlasmidHunter | micromamba/envs/plasmidhunter |
| PlasAnnotatoR | micromamba/envs/plasannotator (Python 3.10) |

### External tools

- blastn / makeblastdb (NCBI BLAST+)
- diamond
- prodigal

## Databases

| Database | Path |
|---|---|
| PLSDB 2025 sequences | data/plsdb/sequences.fasta |
| PLSDB metadata | data/plsdb/meta/ |
| CARD | data/indexes/nucleotide_fasta_protein_homolog_model.fasta |
| MIBiG 4.0 | data/indexes/mibig_4.0.fasta |
| CAZy 2025 | data/indexes/cazy.fasta |
| RF model | data/models/rf_model.pkl |

## Installation

    git clone https://github.com/HaydeSaracho/PlasAnnotatoR.git
    cd PlasAnnotatoR
    micromamba activate plasannotator

## Usage

    # Full pipeline
    python main.py -i input.fasta -o results/ -t 8

    # Skip annotation
    python main.py -i input.fasta -o results/ --skip-annotation -t 8

    # Skip network
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