import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score
from Bio import SeqIO

# Rutas
BASE = Path("/home/bionfo/Escritorio/Hayde/PlasAnnotatoR")
TEST_FASTA = BASE / "data/final/validation_fragmented.fasta"
GROUND_TRUTH = BASE / "data/final/validation_fragmented_ground_truth.tsv"
RESULTS_DIR = BASE / "data/models/tool_evaluations_fragmented"
REPORT_OUTPUT = BASE / "data/models/evaluation_report_fragmented.txt"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_ground_truth():
    print("Cargando ground truth...")
    gt = pd.read_csv(GROUND_TRUTH, sep='\t')
    gt['true_binary'] = (gt['true_label'] == 'plasmid').astype(int)
    y_test = gt['true_binary'].values
    contig_ids = gt['contig_id'].tolist()
    print("Validation dataset: {} secuencias".format(len(y_test)))
    print("Plasmids: {} | Chromosomes: {} | Phages: {}".format(
        (gt['true_label']=='plasmid').sum(),
        (gt['true_label']=='chromosome').sum(),
        (gt['true_label']=='phage').sum()))
    return y_test, contig_ids


def evaluate_rf(y_test, contig_ids):
    print("\n--- Evaluando RF model ---")
    import sys
    sys.path.insert(0, str(BASE))
    from plasannotator.wrappers.rf_model import run_rf_model
    rf_df = run_rf_model(
        input_fasta=str(TEST_FASTA),
        model_path=str(BASE / "data/models/rf_model.pkl")
    )
    rf_df = rf_df.set_index('contig_id')
    scores = [float(rf_df.loc[cid, 'rf_model_score']) if cid in rf_df.index else 0.0 for cid in contig_ids]
    auc = roc_auc_score(y_test, scores)
    print("RF AUC: {:.4f}".format(auc))
    return auc


def evaluate_plasclass(y_test, contig_ids):
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
    cmd = [plasmidhunter, "-i", str(TEST_FASTA), "-o", str(output_dir), "-c", "1"]
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
    y_test, contig_ids = load_ground_truth()
    aucs = {}

    auc_rf = evaluate_rf(y_test, contig_ids)
    if auc_rf: aucs['rf_model'] = auc_rf

    auc_plasclass = evaluate_plasclass(y_test, contig_ids)
    if auc_plasclass: aucs['plasclass'] = auc_plasclass

    auc_plasmidhunter = evaluate_plasmidhunter(y_test, contig_ids)
    if auc_plasmidhunter: aucs['plasmidhunter'] = auc_plasmidhunter

    auc_plasme = evaluate_plasme(y_test, contig_ids)
    if auc_plasme: aucs['plasme'] = auc_plasme

    print("\n=== RESUMEN AUC ===")
    with open(REPORT_OUTPUT, 'w') as f:
        f.write("=== AUC en validation dataset fragmented (n=5365) ===\n\n")
        for tool, auc in sorted(aucs.items(), key=lambda x: x[1], reverse=True):
            line = "{}: {:.4f}".format(tool, auc)
            print(line)
            f.write(line + "\n")

    print("\nReporte guardado en {}".format(REPORT_OUTPUT))


if __name__ == "__main__":
    main()