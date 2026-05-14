"""
Wrapper for PlasmidHunter
Environment: /home/bionfo/micromamba/envs/plasmidhunter
Script: plasmidhunter
Output: directory with predictions.tsv
"""

import subprocess
import pandas as pd
from pathlib import Path


def run_plasmidhunter(input_fasta, output_dir, threads=8, threshold=0.5):
    """
    Runs PlasmidHunter on a FASTA file.

    Args:
        input_fasta: path to input FASTA file
        output_dir: output directory
        threads: number of threads
        threshold: probability threshold (default=0.5)

    Returns:
        DataFrame with columns: contig_id, plasmidhunter_score, plasmidhunter_label
    """
    output_dir = Path(output_dir)

    if output_dir.exists():
        import shutil
        shutil.rmtree(output_dir)

    plasmidhunter = "/home/bionfo/micromamba/envs/plasmidhunter/bin/plasmidhunter"

    cmd = [
        plasmidhunter,
        "-i", str(input_fasta),
        "-o", str(output_dir),
        "-c", str(threads)
    ]

    print("[PlasmidHunter] Running: {}".format(' '.join(cmd)))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        raise RuntimeError("[PlasmidHunter] Error:\n{}".format(result.stderr))

    result_file = output_dir / "predictions.tsv"
    if not result_file.exists():
        raise FileNotFoundError("[PlasmidHunter] File not found: {}".format(result_file))

    df = pd.read_csv(result_file, sep="\t", index_col=0)
    df.index.name = "contig_id"
    df = df.reset_index()
    df = df.rename(columns={
        "Prediction (0: chromosome, 1: plasmid)": "prediction",
        "Probability of 0": "prob_chromosome",
        "Probability of 1": "plasmidhunter_score"
    })
    df["plasmidhunter_score"] = pd.to_numeric(df["plasmidhunter_score"], errors="coerce").fillna(0.0)
    df = df.drop_duplicates(subset="contig_id", keep="first")
    df["plasmidhunter_label"] = df["plasmidhunter_score"].apply(
        lambda x: "plasmid" if x >= threshold else "chromosome"
    )

    n_plasmid = (df["plasmidhunter_label"] == "plasmid").sum()
    print("[PlasmidHunter] Done. {} contigs processed. {} plasmids.".format(len(df), n_plasmid))
    return df[["contig_id", "plasmidhunter_score", "plasmidhunter_label"]]


if __name__ == "__main__":
    df = run_plasmidhunter(
        input_fasta="/home/bionfo/PLASMe/test.fasta",
        output_dir="/tmp/plasmidhunter_test"
    )
    print(df.head(10))