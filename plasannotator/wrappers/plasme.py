"""
Wrapper for PLASMe
Environment: /home/bionfo/anaconda3/envs/plasme
Script: /home/bionfo/PLASMe/PLASMe.py
Output: <output>_report.csv with contigs classified as plasmids
"""

import subprocess
import pandas as pd
from pathlib import Path


def run_plasme(input_fasta, output_dir, threads=8, threshold=0.5):
    """
    Runs PLASMe on a FASTA file.

    Args:
        input_fasta: path to input FASTA file
        output_dir: output directory
        threads: number of threads
        threshold: probability threshold (default=0.5)

    Returns:
        DataFrame with columns: contig_id, plasme_score, plasme_label
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "plasme_output"
    report_file = Path(str(output_file) + "_report.csv")

    # Remove if exists (PLASMe fails if output already exists)
    if output_file.exists():
        output_file.unlink()
    if report_file.exists():
        report_file.unlink()

    python = "/home/bionfo/anaconda3/envs/plasme/bin/python"
    script = "/home/bionfo/PLASMe/PLASMe.py"
    database = "/home/bionfo/PLASMe/DB"

    cmd = [
        python, script,
        str(input_fasta),
        str(output_file),
        "-d", database,
        "-t", str(threads),
        "-p", str(threshold)
    ]

    print("[PLASMe] Running: {}".format(' '.join(cmd)))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        raise RuntimeError("[PLASMe] Error:\n{}".format(result.stderr))

    if not report_file.exists():
        raise FileNotFoundError("[PLASMe] Report not found: {}".format(report_file))

    # Read report - only contains contigs classified as plasmids
    df_plasmid = pd.read_csv(report_file, sep="\t")
    df_plasmid = df_plasmid.rename(columns={"contig": "contig_id", "score": "plasme_score"})
    df_plasmid["plasme_label"] = "plasmid"

    # Read all contigs from FASTA to include chromosomal ones
    from Bio import SeqIO
    all_contigs = [rec.id for rec in SeqIO.parse(str(input_fasta), "fasta")]
    df_all = pd.DataFrame({"contig_id": all_contigs})

    # Merge to assign chromosome label to undetected contigs
    df = df_all.merge(df_plasmid[["contig_id", "plasme_score", "plasme_label"]], on="contig_id", how="left")
    df["plasme_label"] = df["plasme_label"].fillna("chromosome")
    df["plasme_score"] = df["plasme_score"].fillna(0.0)

    n_plasmid = df["plasme_label"].eq("plasmid").sum()
    print("[PLASMe] Done. {} contigs processed. {} plasmids.".format(len(df), n_plasmid))
    return df[["contig_id", "plasme_score", "plasme_label"]]


if __name__ == "__main__":
    df = run_plasme(
        input_fasta="/home/bionfo/PLASMe/test.fasta",
        output_dir="/tmp/plasme_test3"
    )
    print(df.head(10))