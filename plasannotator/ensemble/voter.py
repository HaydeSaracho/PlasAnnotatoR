"""
Weighted voting module for PlasAnnotatoR v2.0
Combines predictions from PlasClass, PLASMe, PlasmidHunter and custom RF model
Weights based on AUC of each tool on independent validation datasets
(full-length sequences and fragmented contigs)
"""

import pandas as pd
import numpy as np
import yaml
from pathlib import Path


def load_config(config_path="config.yaml"):
    """Loads project configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def weighted_vote(dfs, weights, threshold=0.6):
    """
    Combines predictions from multiple tools using weighted voting.

    Args:
        dfs: dict with {tool_name: DataFrame}
        weights: dict with {tool_name: weight}
        threshold: ensemble score threshold (default: 0.6)

    Returns:
        DataFrame with ensemble results
    """
    total_weight = sum(weights[tool] for tool in dfs.keys())
    norm_weights = {tool: w / total_weight for tool, w in weights.items() if tool in dfs}

    merged = None
    for tool, df in dfs.items():
        score_col = "{}_score".format(tool)
        label_col = "{}_label".format(tool)
        df_clean = df[["contig_id", score_col, label_col]].copy()
        if merged is None:
            merged = df_clean
        else:
            merged = merged.merge(df_clean, on="contig_id", how="outer")

    ensemble_score = np.zeros(len(merged))
    for tool, weight in norm_weights.items():
        score_col = "{}_score".format(tool)
        if score_col in merged.columns:
            scores = merged[score_col].fillna(0.0).values
            ensemble_score += weight * scores

    merged["ensemble_score"] = ensemble_score
    merged["ensemble_label"] = merged["ensemble_score"].apply(
        lambda x: "plasmid" if x >= threshold else "chromosome"
    )

    n_plasmid = (merged["ensemble_label"] == "plasmid").sum()
    n_chrom = (merged["ensemble_label"] == "chromosome").sum()
    print("[Ensemble] Total contigs: {}".format(len(merged)))
    print("[Ensemble] Plasmids: {} | Chromosomes: {} | Threshold: {}".format(
        n_plasmid, n_chrom, threshold))

    return merged


def run_ensemble(input_fasta, output_dir, config_path="config.yaml",
                 threads=8, threshold=0.6):
    """
    Runs all tools and combines results.

    Args:
        input_fasta: path to input FASTA file
        output_dir: output directory
        config_path: path to config.yaml
        threads: number of threads
        threshold: ensemble score threshold (default: 0.6)

    Returns:
        DataFrame with ensemble results
    """
    from plasannotator.wrappers.plasclass import run_plasclass
    from plasannotator.wrappers.plasme import run_plasme
    from plasannotator.wrappers.plasmidhunter import run_plasmidhunter
    from plasannotator.wrappers.rf_model import run_rf_model

    config = load_config(config_path)
    weights = config['ensemble']['weights']
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    print("\n--- Running PlasClass ---")
    try:
        results['plasclass'] = run_plasclass(
            input_fasta=input_fasta,
            output_dir=str(output_dir / "plasclass"),
            threads=threads
        )
    except Exception as e:
        print("[PlasClass] Error: {}".format(e))

    print("\n--- Running PLASMe ---")
    try:
        results['plasme'] = run_plasme(
            input_fasta=input_fasta,
            output_dir=str(output_dir / "plasme"),
            threads=threads
        )
    except Exception as e:
        print("[PLASMe] Error: {}".format(e))

    print("\n--- Running PlasmidHunter ---")
    try:
        results['plasmidhunter'] = run_plasmidhunter(
            input_fasta=input_fasta,
            output_dir=str(output_dir / "plasmidhunter"),
            threads=threads
        )
    except Exception as e:
        print("[PlasmidHunter] Error: {}".format(e))

    print("\n--- Running RF model ---")
    try:
        results['rf_model'] = run_rf_model(
            input_fasta=input_fasta,
            model_path=config['environments']['rf_model']['model_path']
        )
    except Exception as e:
        print("[RF] Error: {}".format(e))

    print("\n--- Weighted voting ---")
    ensemble_df = weighted_vote(results, weights, threshold=threshold)

    output_file = output_dir / "ensemble_results.tsv"
    ensemble_df.to_csv(output_file, sep='\t', index=False)
    print("\n[Ensemble] Results saved to {}".format(output_file))

    return ensemble_df