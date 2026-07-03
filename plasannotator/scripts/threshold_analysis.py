"""
Threshold analysis for PlasAnnotatoR ensemble classifier.
Evaluates precision, recall and F1 across thresholds 0.1-0.9
for both full-length and fragmented validation sets.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import precision_score, recall_score, f1_score

BASE = Path("/home/bionfo/Escritorio/Hayde/PlasAnnotatoR")

# Pesos del ensemble (promedio full-length + fragmented)
WEIGHTS = {
    'plasclass': 0.8421,
    'plasme': 0.8895,
    'plasmidhunter': 0.7898,
    'rf_model': 0.8809
}

SCENARIOS = {
    'full_length': {
        'fasta': BASE / "data/final/validation_full_length.fasta",
        'ground_truth': BASE / "data/final/validation_full_length_ground_truth.tsv",
        'eval_dir': BASE / "data/models/tool_evaluations_full_length",
        'scores_output': BASE / "data/models/rf_scores_full_length.tsv",
    },
    'fragmented': {
        'fasta': BASE / "data/final/validation_fragmented.fasta",
        'ground_truth': BASE / "data/final/validation_fragmented_ground_truth.tsv",
        'eval_dir': BASE / "data/models/tool_evaluations_fragmented",
        'scores_output': BASE / "data/models/rf_scores_fragmented.tsv",
    }
}


def load_ground_truth(path):
    gt = pd.read_csv(path, sep='\t')
    gt['true_binary'] = (gt['true_label'] == 'plasmid').astype(int)
    return gt.set_index('contig_id')


def load_plasclass_scores(eval_dir):
    df = pd.read_csv(eval_dir / "plasclass_output", sep="\t", header=None,
                     names=["contig_id", "score"])
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.0)
    return df.drop_duplicates("contig_id").set_index("contig_id")["score"]


def load_plasme_scores(eval_dir):
    df = pd.read_csv(eval_dir / "plasme_output_report.csv", sep="\t")
    df = df.rename(columns={"contig": "contig_id", "score": "plasme_score"})
    df["plasme_score"] = pd.to_numeric(df["plasme_score"], errors="coerce").fillna(0.0)
    return df.drop_duplicates("contig_id").set_index("contig_id")["plasme_score"]


def load_plasmidhunter_scores(eval_dir):
    df = pd.read_csv(eval_dir / "plasmidhunter_output/predictions.tsv", sep="\t", index_col=0)
    df.index.name = "contig_id"
    df = df.reset_index().rename(columns={"Probability of 1": "score"})
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.0)
    return df.drop_duplicates("contig_id").set_index("contig_id")["score"]


def get_rf_scores(fasta_path, scores_output):
    if scores_output.exists():
        print("  RF scores already computed, loading from {}".format(scores_output))
        df = pd.read_csv(scores_output, sep="\t")
        return df.set_index("contig_id")["rf_score"]
    print("  Running RF model...")
    sys.path.insert(0, str(BASE))
    from plasannotator.wrappers.rf_model import run_rf_model
    rf_df = run_rf_model(
        input_fasta=str(fasta_path),
        model_path=str(BASE / "data/models/rf_model.pkl")
    )
    rf_df = rf_df.rename(columns={"rf_model_score": "rf_score"})
    rf_df[["contig_id", "rf_score"]].to_csv(scores_output, sep="\t", index=False)
    print("  RF scores saved to {}".format(scores_output))
    return rf_df.set_index("contig_id")["rf_score"]


def compute_ensemble_scores(contig_ids, scores_dict, weights):
    total_weight = sum(weights.values())
    norm_weights = {k: v / total_weight for k, v in weights.items()}
    ensemble = pd.Series(0.0, index=contig_ids)
    for tool, weight in norm_weights.items():
        tool_scores = scores_dict[tool].reindex(contig_ids).fillna(0.0)
        ensemble += weight * tool_scores
    return ensemble


def threshold_analysis(y_true, ensemble_scores, thresholds):
    results = []
    for t in thresholds:
        y_pred = (ensemble_scores >= t).astype(int)
        results.append({
            'threshold': t,
            'precision': precision_score(y_true, y_pred, zero_division=0),
            'recall': recall_score(y_true, y_pred, zero_division=0),
            'f1': f1_score(y_true, y_pred, zero_division=0)
        })
    return pd.DataFrame(results)


def main():
    thresholds = np.round(np.arange(0.1, 1.0, 0.1), 2)
    all_results = {}

    for scenario_name, scenario in SCENARIOS.items():
        print("\n=== {} ===".format(scenario_name.upper()))

        gt = load_ground_truth(scenario['ground_truth'])
        contig_ids = gt.index.tolist()
        y_true = gt['true_binary'].values

        print("  Loading tool scores...")
        scores_dict = {
            'plasclass': load_plasclass_scores(scenario['eval_dir']),
            'plasme': load_plasme_scores(scenario['eval_dir']),
            'plasmidhunter': load_plasmidhunter_scores(scenario['eval_dir']),
            'rf_model': get_rf_scores(scenario['fasta'], scenario['scores_output'])
        }

        ensemble_scores = compute_ensemble_scores(contig_ids, scores_dict, WEIGHTS)

        df = threshold_analysis(y_true, ensemble_scores, thresholds)
        all_results[scenario_name] = df

        print("\n  Threshold analysis:")
        print(df.to_string(index=False))

        best_f1_idx = df['f1'].idxmax()
        best = df.loc[best_f1_idx]
        print("\n  Best F1: {:.4f} at threshold {:.1f} (precision={:.4f}, recall={:.4f})".format(
            best['f1'], best['threshold'], best['precision'], best['recall']))
        print("  At threshold 0.5: F1={:.4f}, precision={:.4f}, recall={:.4f}".format(
            df[df['threshold']==0.5]['f1'].values[0],
            df[df['threshold']==0.5]['precision'].values[0],
            df[df['threshold']==0.5]['recall'].values[0]))

        df.to_csv(BASE / "data/models/threshold_analysis_{}.tsv".format(scenario_name),
                  sep="\t", index=False)

    # Figura comparativa
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {'precision': '#2196F3', 'recall': '#F44336', 'f1': '#4CAF50'}
    titles = {'full_length': 'Full-length sequences', 'fragmented': 'Fragmented contigs'}

    for ax, (scenario_name, df) in zip(axes, all_results.items()):
        for metric, color in colors.items():
            ax.plot(df['threshold'], df[metric], marker='o', label=metric.capitalize(),
                    color=color, linewidth=2)
        ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.7, label='Default (0.5)')
        ax.set_xlabel('Threshold', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title(titles[scenario_name], fontsize=13)
        ax.legend(fontsize=10)
        ax.set_xlim(0.05, 0.95)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)

    plt.suptitle('PlasAnnotatoR ensemble threshold analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(BASE / "data/models/threshold_analysis_comparison.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(BASE / "data/models/threshold_analysis_comparison.png", dpi=300, bbox_inches='tight')
    print("\nFigure saved to data/models/threshold_analysis_comparison.pdf/png")


if __name__ == "__main__":
    main()