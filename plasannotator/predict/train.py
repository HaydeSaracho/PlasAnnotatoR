"""
Entrenamiento del modelo Random Forest para predicción de plásmidos.
"""

import pickle
import numpy as np
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import (
    classification_report, roc_auc_score,
    precision_score, recall_score, f1_score
)
from plasannotator.predict.kmer import extract_features
from plasannotator.config import logger

MODEL_DIR = Path(__file__).parent.parent.parent / "data" / "models"
MODEL_DIR.mkdir(parents=True, exist_ok=True)

PLASMID_FASTA  = Path("/mnt/d/PlasAnnotatoR/data/indexes/all_plasmids_nr_full.fasta")
CHROM_FASTA    = Path("/mnt/d/PlasAnnotatoR/data/chromosomes_fragments.fasta")


def load_dataset():
    """Carga y etiqueta positivos (1) y negativos (0)."""
    print("Cargando positivos (plasmidos)...")
    _, _, X_pos = extract_features(PLASMID_FASTA)
    y_pos = np.ones(len(X_pos), dtype=np.int8)
    print(f"  {len(X_pos):,} secuencias plasmídicas")

    print("Cargando negativos (cromosomas)...")
    _, _, X_neg = extract_features(CHROM_FASTA)
    y_neg = np.zeros(len(X_neg), dtype=np.int8)
    print(f"  {len(X_neg):,} fragmentos cromosómicos")

    X = np.vstack([X_pos, X_neg])
    y = np.concatenate([y_pos, y_neg])
    print(f"  Dataset total: {len(y):,} secuencias")
    return X, y


def train():
    X, y = load_dataset()

    # Split 80/20 — el 20% nunca toca el entrenamiento
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    print(f"\nTrain: {len(y_train):,} | Test: {len(y_test):,}")

    # k-fold estratificado sobre train
    print("\nValidación k-fold (k=5)...")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_metrics = []

    for fold, (tr_idx, val_idx) in enumerate(skf.split(X_train, y_train), 1):
        clf = RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            n_jobs=-1,
            random_state=42
        )
        clf.fit(X_train[tr_idx], y_train[tr_idx])
        y_pred = clf.predict(X_train[val_idx])
        y_prob = clf.predict_proba(X_train[val_idx])[:, 1]

        metrics = {
            "precision": precision_score(y_train[val_idx], y_pred),
            "recall":    recall_score(y_train[val_idx], y_pred),
            "f1":        f1_score(y_train[val_idx], y_pred),
            "auc":       roc_auc_score(y_train[val_idx], y_prob),
        }
        fold_metrics.append(metrics)
        print(f"  Fold {fold} — P:{metrics['precision']:.3f} R:{metrics['recall']:.3f} "
              f"F1:{metrics['f1']:.3f} AUC:{metrics['auc']:.3f}")

    # Promedio k-fold
    print("\nPromedio k-fold:")
    for m in ["precision", "recall", "f1", "auc"]:
        vals = [f[m] for f in fold_metrics]
        print(f"  {m}: {np.mean(vals):.3f} ± {np.std(vals):.3f}")

    # Modelo final sobre todo el train
    print("\nEntrenando modelo final...")
    clf_final = RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",
        n_jobs=-1,
        random_state=42
    )
    clf_final.fit(X_train, y_train)

    # Evaluación sobre test independiente
    print("\nEvaluación sobre dataset de prueba independiente:")
    y_pred_test = clf_final.predict(X_test)
    y_prob_test = clf_final.predict_proba(X_test)[:, 1]
    print(classification_report(y_test, y_pred_test,
                                target_names=["chromosome", "plasmid"]))
    print(f"AUC-ROC: {roc_auc_score(y_test, y_prob_test):.4f}")

    # Guardar modelo
    model_path = MODEL_DIR / "rf_model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump(clf_final, f)
    print(f"\nModelo guardado: {model_path}")


if __name__ == "__main__":
    train()