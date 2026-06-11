"""
Evaluación de PlasClass, PLASMe, PlasmidHunter y RF model en el benchmark dataset.
Dataset independiente: 487 plásmidos de RefSeq/PlasmidScope + 487 cromosomas de Chromosome level.
Ground truth conocido en data/benchmark_plasmids/ground_truth.tsv
"""

import pickle
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score
from Bio import SeqIO

# Rutas
BASE = Path("/home/bionfo/Escritorio/Hayde/PlasAnnotatoR")
TEST_FASTA = BASE / "data/benchmark_plasmids/benchmark_974.fasta"
GROUND_TRUTH = BASE / "data/benchmark_plasmids/ground_truth.tsv"
TEST_PKL = BASE / "data/models/test_set.pkl"
RESULTS_DIR = BASE / "data/models/tool_evaluations"
REPORT_OUTPUT = BASE / "data/models/evaluation_report.txt"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_test_set():
    print("Cargando ground truth...")
    gt = pd.read_csv(GROUND_TRUTH, sep='\t')
    gt['true_binary'] = (gt['true_label'] == 'plasmid').astype(int)
    y_test = gt['true_binary'].values
    contig_ids = gt['contig_id'].tolist()
    print("Benchmark dataset: {} secuencias".format(len(y_test)))
    print("Plasmids: {} | Chromosomes: {}".format(
        sum(y_test), len(y_test) - sum(y_test)))

    # RF model - correr sobre el benchmark dataset
    print("Evaluando RF model...")
    from plasannotator.wrappers.rf_model import run_rf_model
    import sys
    sys.path.insert(0, str(BASE))
    rf_df = run_rf_model(
        input_fasta=str(TEST_FASTA),
        model_path=str(BASE / "data/models/rf_model.pkl")
    )
    rf_df = rf_df.set_index('contig_id')
    rf_scores = [float(rf_df.loc[cid, 'rf_model_score']) if cid in rf_df.index else 0.0 for cid in contig_ids]
    rf_auc = roc_auc_score(y_test, rf_scores)
    print("RF AUC: {:.4f}".format(rf_auc))

    return y_test, rf_scores, rf_auc, contig_ids


def get_contig_order():
    print("Leyendo orden de contigs...")
    ids = [rec.id for rec in SeqIO.parse(str(TEST_FASTA), "fasta")]
    print("Contigs en FASTA: {}".format(len(ids)))
    return ids


def evaluate_plasclass(y_test, contig_ids):
    """Corre PlasClass sobre el benchmark dataset y calcula AUC."""
    print("\n--- Evaluando PlasClass ---")
    output = RESULTS_DIR / "plasclass_output"

    python = "/home/bionfo/micromamba/envs/plasclass/bin/python"
    script = "/home/bionfo/micromamba/envs/plasclass/bin/classify_fasta.py"
    cmd = [python, script, "-f", str(TEST_FASTA), "-o", str(output), "-p", "8"]

    print("Ejecutando PlasClass...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        print("Error PlasClass: {}".format(result.stderr[:300]))
        return None

    df = pd.read_csv(output, sep="\t", header=None, names=["contig_id", "score"])
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.0)
    df = df.drop_duplicates(subset="contig_id", keep="first")
    df = df.set_index("contig_id")
    scores = [float(df.loc[cid, "score"]) if cid in df.index else 0.0 for cid in contig_ids]
    auc = roc_auc_score(y_test, scores)
    print("PlasClass AUC: {:.4f}".format(auc))
    return auc


def evaluate_plasmidhunter(y_test, contig_ids):
    print("\n--- Evaluando PlasmidHunter ---")
    output_dir = RESULTS_DIR / "plasmidhunter_output"

    if output_dir.exists():
        import shutil
        shutil.rmtree(output_dir)

    plasmidhunter = "/home/bionfo/micromamba/envs/plasmidhunter/bin/plasmidhunter"
    cmd = [plasmidhunter, "-i", str(TEST_FASTA), "-o", str(output_dir), "-c", "8"]

    print("Ejecutando PlasmidHunter...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        print("Error PlasmidHunter: {}".format(result.stderr[:300]))
        return None

    result_file = output_dir / "predictions.tsv"
    df = pd.read_csv(result_file, sep="\t", index_col=0)
    df.index.name = "contig_id"
    df = df.reset_index()
    df = df.rename(columns={"Probability of 1": "score"})
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.0)
    df = df.drop_duplicates(subset="contig_id", keep="first")
    df = df.set_index("contig_id")
    scores = [float(df.loc[cid, "score"]) if cid in df.index else 0.0 for cid in contig_ids]
    auc = roc_auc_score(y_test, scores)
    print("PlasmidHunter AUC: {:.4f}".format(auc))
    return auc


def evaluate_plasme(y_test, contig_ids):
    print("\n--- Evaluando PLASMe ---")
    output_file = RESULTS_DIR / "plasme_output"
    report_file = RESULTS_DIR / "plasme_output_report.csv"

    if output_file.exists():
        output_file.unlink()
    if report_file.exists():
        report_file.unlink()

    python = "/home/bionfo/miniforge3/envs/plasme/bin/python"
    script = "/home/bionfo/PLASMe/PLASMe.py"
    database = "/home/bionfo/PLASMe/DB"
    cmd = [python, script, str(TEST_FASTA), str(output_file),
           "-d", database, "-t", "8", "-p", "0.0"]

    print("Ejecutando PLASMe...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        print("Error PLASMe: {}".format(result.stderr[:300]))
        return None

    df = pd.read_csv(report_file, sep="\t")
    df = df.rename(columns={"contig": "contig_id", "score": "plasme_score"})
    df["plasme_score"] = pd.to_numeric(df["plasme_score"], errors="coerce").fillna(0.0)
    df = df.drop_duplicates(subset="contig_id", keep="first")
    df = df.set_index("contig_id")
    scores = [float(df.loc[cid, "plasme_score"]) if cid in df.index else 0.0 for cid in contig_ids]
    auc = roc_auc_score(y_test, scores)
    print("PLASMe AUC: {:.4f}".format(auc))
    return auc


def main():
    y_test, rf_scores, rf_auc, contig_ids = load_test_set()

    aucs = {'rf_model': rf_auc}

    auc_plasclass = evaluate_plasclass(y_test, contig_ids)
    if auc_plasclass:
        aucs['plasclass'] = auc_plasclass

    auc_plasmidhunter = evaluate_plasmidhunter(y_test, contig_ids)
    if auc_plasmidhunter:
        aucs['plasmidhunter'] = auc_plasmidhunter

    auc_plasme = evaluate_plasme(y_test, contig_ids)
    if auc_plasme:
        aucs['plasme'] = auc_plasme

    print("\n=== RESUMEN AUC ===")
    with open(REPORT_OUTPUT, 'w') as f:
        f.write("=== AUC en benchmark dataset independiente (n=974) ===\n\n")
        for tool, auc in sorted(aucs.items(), key=lambda x: x[1], reverse=True):
            line = "{}: {:.4f}".format(tool, auc)
            print(line)
            f.write(line + "\n")
        f.write("\nActualizar config.yaml con estos valores.\n")

    print("\nReporte guardado en {}".format(REPORT_OUTPUT))


if __name__ == "__main__":
    main()