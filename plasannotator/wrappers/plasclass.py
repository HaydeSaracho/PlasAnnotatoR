"""
Wrapper for PlasClass
Environment: plasclass (Python 3.7 + scikit-learn 0.21.3)
Script: classify_fasta.py
Output: tab-separated file, no header
"""

import subprocess
import pandas as pd
import yaml
from pathlib import Path


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def run_plasclass(input_fasta: str, output_dir: str, threads: int = 8, config_path="config.yaml") -> pd.DataFrame:
    """
    Runs PlasClass on a FASTA file.

    Args:
        input_fasta: path to input FASTA file
        output_dir: output directory
        threads: number of processes
        config_path: path to config.yaml

    Returns:
        DataFrame with columns: contig_id, plasclass_score, plasclass_label
    """
    config = load_config(config_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "plasclass_output"

    python = config['environments']['plasclass']['python']
    script = config['environments']['plasclass']['script']

    cmd = [
        python, script,
        "-f", str(input_fasta),
        "-o", str(output_file),
        "-p", str(threads)
    ]

    print("[PlasClass] Running: {}".format(' '.join(cmd)))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError("[PlasClass] Error:\n{}".format(result.stderr))

    df = pd.read_csv(output_file, sep="\t", header=None, names=["contig_id", "plasclass_score"])
    df["plasclass_score"] = pd.to_numeric(df["plasclass_score"], errors="coerce").fillna(0.0)
    df = df.drop_duplicates(subset="contig_id", keep="first")
    df["plasclass_label"] = df["plasclass_score"].apply(
        lambda x: "plasmid" if x >= 0.5 else "chromosome"
    )

    n_plasmid = (df["plasclass_label"] == "plasmid").sum()
    print("[PlasClass] Done. {} contigs processed. {} plasmids.".format(len(df), n_plasmid))
    return df


if __name__ == "__main__":
    df = run_plasclass(
        input_fasta="/tmp/test.fa",
        output_dir="/tmp/plasclass_test"
    )
    print(df.head())
