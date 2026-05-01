"""
Motor de predicción ensemble: k-mer + alineamiento Minimap2.
"""

import csv
import json
import pickle
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from Bio import SeqIO

from plasannotator.config import (
    CONFIDENCE_THRESHOLD,
    MIN_CONTIG_LENGTH,
    MINIMAP2_THREADS,
    logger,
)
from plasannotator.db.manager import get_index_path
from plasannotator.predict.kmer import extract_features, KMER_SIZE


# ---------------------------------------------------------------------------
# Clasificador ML (modelo entrenado)
# ---------------------------------------------------------------------------

_MODEL_PATH = Path(__file__).parent.parent.parent / "data" / "models" / "rf_model.pkl"
_model = None


def _load_model():
    """Carga el modelo RF entrenado (singleton)."""
    global _model
    if _model is None:
        if not _MODEL_PATH.exists():
            logger.warning("Modelo RF no encontrado — usando heurística provisional.")
            return None
        with open(_MODEL_PATH, "rb") as f:
            _model = pickle.load(f)
        logger.info(f"Modelo RF cargado: {_MODEL_PATH}")
    return _model


def _kmer_score(features: np.ndarray) -> np.ndarray:
    """
    Score de plasmidio basado en el modelo RF entrenado.
    Si el modelo no está disponible, usa heurística de varianza.
    """
    model = _load_model()
    if model is not None:
        return model.predict_proba(features)[:, 1].astype(np.float32)
    variances = np.var(features, axis=1)
    scores = (variances - variances.min()) / (variances.max() - variances.min() + 1e-9)
    return scores.astype(np.float32)


# ---------------------------------------------------------------------------
# Alineamiento Minimap2
# ---------------------------------------------------------------------------

def _run_minimap2(fasta_path: Path, index_path: Path, threads: int) -> dict[str, float]:
    """
    Alinea los contigs contra el índice de la DB con Minimap2.
    Retorna un dict {contig_id: alignment_score} normalizado (0-1).
    """
    cmd = [
        "minimap2",
        "-c",
        "--secondary=no",
        "--split-prefix=/tmp/plasannotator_split",
        "-t", str(threads),
        str(index_path),
        str(fasta_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.warning(f"Minimap2 warning: {result.stderr[:200]}")

    scores: dict[str, float] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        contig_id = parts[0]
        query_len = int(parts[1])
        matches = int(parts[9])
        aln_len = int(parts[10])
        if aln_len > 0:
            score = (matches / aln_len) * min(aln_len / query_len, 1.0)
            if contig_id not in scores or score > scores[contig_id]:
                scores[contig_id] = round(score, 4)

    return scores


# ---------------------------------------------------------------------------
# Ensemble
# ---------------------------------------------------------------------------

def _ensemble_score(
    kmer_scores: np.ndarray,
    aln_scores: dict[str, float],
    contig_ids: list[str],
    w_kmer: float = 0.4,
    w_aln: float = 0.6,
) -> np.ndarray:
    """
    Combina score k-mer y score de alineamiento con pesos.
    Si no hay hit de alineamiento, usa solo k-mer.
    """
    final = np.zeros(len(contig_ids), dtype=np.float32)
    for i, cid in enumerate(contig_ids):
        aln = aln_scores.get(cid, None)
        if aln is None:
            final[i] = kmer_scores[i]
        else:
            final[i] = w_kmer * kmer_scores[i] + w_aln * aln
    return final


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------

def _write_tsv(
    output_prefix: str,
    contig_ids: list[str],
    lengths: list[int],
    scores: np.ndarray,
    aln_scores: dict[str, float],
    db_name: str,
    threshold: float,
) -> Path:
    """Escribe el archivo TSV de predicciones."""
    out_path = Path(f"{output_prefix}.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "contig_id", "length_bp", "prediction",
            "confidence", "aln_score", "db"
        ])
        for cid, length, score in zip(contig_ids, lengths, scores):
            prediction = "plasmid" if score >= threshold else "chromosome"
            aln = aln_scores.get(cid, "NA")
            writer.writerow([cid, length, prediction, round(float(score), 4), aln, db_name])

    logger.info(f"TSV escrito: {out_path}")
    return out_path


def _write_json(
    output_prefix: str,
    contig_ids: list[str],
    lengths: list[int],
    scores: np.ndarray,
    aln_scores: dict[str, float],
    db_name: str,
    threshold: float,
) -> Path:
    """Escribe el archivo JSON de predicciones."""
    out_path = Path(f"{output_prefix}.json")
    results = []
    for cid, length, score in zip(contig_ids, lengths, scores):
        results.append({
            "contig_id": cid,
            "length_bp": length,
            "prediction": "plasmid" if score >= threshold else "chromosome",
            "confidence": round(float(score), 4),
            "aln_score": aln_scores.get(cid, None),
            "db": db_name,
        })

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"JSON escrito: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Punto de entrada principal
# ---------------------------------------------------------------------------

def run_prediction(
    fasta_path: Path,
    db_name: str,
    output_prefix: str,
    threads: int,
    min_length: int,
    build_network: bool,
) -> None:
    """Pipeline completo de predicción."""

    # 1. Verificar índice de DB
    logger.info(f"[1/4] Verificando índice de DB '{db_name}'...")
    index_path = get_index_path(db_name)

    # 2. Extraer features k-mer
    logger.info(f"[2/4] Extrayendo features k-mer (k={KMER_SIZE})...")
    contig_ids, lengths, features = extract_features(fasta_path)

    # 3. Clasificación ML
    logger.info(f"[3/4] Clasificando con modelo RF...")
    kmer_scores = _kmer_score(features)

    # 4. Alineamiento Minimap2
    logger.info(f"[4/4] Alineando contra '{db_name}' con Minimap2...")
    aln_scores = _run_minimap2(fasta_path, index_path, threads)

    # 5. Ensemble
    final_scores = _ensemble_score(kmer_scores, aln_scores, contig_ids)

    # Escribir salidas
    tsv_path = _write_tsv(
        output_prefix, contig_ids, lengths,
        final_scores, aln_scores, db_name, CONFIDENCE_THRESHOLD
    )
    _write_json(
        output_prefix, contig_ids, lengths,
        final_scores, aln_scores, db_name, CONFIDENCE_THRESHOLD
    )

    # Resumen
    n_plasmids = sum(1 for s in final_scores if s >= CONFIDENCE_THRESHOLD)
    n_chrom = len(final_scores) - n_plasmids

    print(f"\n  Contigs analizados : {len(contig_ids):,}")
    print(f"  Plasmidios         : {n_plasmids:,}")
    print(f"  Cromosoma          : {n_chrom:,}")
    print(f"  Salida TSV         : {tsv_path}")

    # Red de secuencias
    if build_network:
        from plasannotator.network.builder import build_network as _build_network
        from plasannotator.config import ANI_THRESHOLD
        _build_network(
            tsv_path=tsv_path,
            db_name=db_name,
            output_prefix=output_prefix.replace("predictions", "network"),
            ani_threshold=ANI_THRESHOLD,
            threads=threads,
        )