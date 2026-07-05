"""
Wrapper for custom Random Forest model
Model: data/models/rf_model.pkl (AUC 0.9913 on internal test set)
Trained on plasmids from PLSDB 2025, chromosomal fragments from RefSeq
and bacteriophage sequences from NCBI RefSeq viral
"""

import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import product
from Bio import SeqIO


def get_all_kmers(k=5):
    """Generates all possible k-mers."""
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=k)]


ALL_KMERS = get_all_kmers(5)
KMER_INDEX = {kmer: i for i, kmer in enumerate(ALL_KMERS)}


def sequence_to_vector(sequence):
    """Converts a sequence to a k-mer frequency vector."""
    sequence = sequence.upper()
    vector = np.zeros(len(ALL_KMERS))
    total = 0
    for i in range(len(sequence) - 5 + 1):
        kmer = sequence[i:i+5]
        if 'N' not in kmer and kmer in KMER_INDEX:
            vector[KMER_INDEX[kmer]] += 1
            total += 1
    if total > 0:
        vector = vector / total
    return vector


def run_rf_model(input_fasta, model_path="data/models/rf_model.pkl", threshold=0.5):
    """
    Classifies contigs using the custom Random Forest model.

    Args:
        input_fasta: path to input FASTA file
        model_path: path to .pkl model file
        threshold: probability threshold (default=0.5)

    Returns:
        DataFrame with columns: contig_id, rf_model_score, rf_model_label
    """
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError("[RF] Model not found: {}".format(model_path))

    print("[RF] Loading model from {}".format(model_path))
    with open(model_path, "rb") as f:
        model = pickle.load(f)

    print("[RF] Reading sequences from {}".format(input_fasta))
    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    contig_ids = [rec.id for rec in records]
    sequences = [str(rec.seq) for rec in records]

    print("[RF] Computing k-mer frequencies ({} sequences)...".format(len(sequences)))
    X = np.array([sequence_to_vector(seq) for seq in sequences])

    print("[RF] Predicting...")
    probas = model.predict_proba(X)
    rf_scores = probas[:, 1]

    df = pd.DataFrame({
        "contig_id": contig_ids,
        "rf_model_score": rf_scores
    })
    df["rf_model_label"] = df["rf_model_score"].apply(
        lambda x: "plasmid" if x >= threshold else "chromosome"
    )

    n_plasmid = (df["rf_model_label"] == "plasmid").sum()
    print("[RF] Done. {} contigs processed. {} plasmids.".format(len(df), n_plasmid))
    return df


if __name__ == "__main__":
    df = run_rf_model(
        input_fasta="/home/bionfo/PLASMe/test.fasta",
        model_path="data/models/rf_model.pkl"
    )
    print(df.head(10))