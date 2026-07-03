"""
Script de entrenamiento del modelo Random Forest para PlasAnnotatoR v2.0
Datos positivos: plásmidos PLSDB 2025 (todas las secuencias, sin filtro de tamaño)
Datos negativos: fragmentos cromosómicos de RefSeq (Chromosome level) + secuencias de fagos de RefSeq
Features: frecuencias de 5-mers canonicos
"""

import pickle
import random
import numpy as np
from pathlib import Path
from itertools import product
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, classification_report
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
log = logging.getLogger(__name__)

# Rutas
BASE = Path("/home/bionfo/Escritorio/Hayde/PlasAnnotatoR")
PLASMIDS_FASTA = BASE / "data/plsdb/sequences.fasta"
CHROMOSOMES_DIR = BASE / "data/chromosomes"
PHAGES_FASTA = BASE / "data/phages.fasta"
MODEL_OUTPUT = BASE / "data/models/rf_model.pkl"
TEST_SET_OUTPUT = BASE / "data/models/test_set.pkl"
TEST_FASTA_OUTPUT = BASE / "data/models/test_set.fasta"
REPORT_OUTPUT = BASE / "data/models/training_report.txt"

# Parámetros
K = 5
RANDOM_SEED = 42
N_JOBS = 8

# Rangos de fragmentos cromosómicos (pb)
FRAGMENT_SIZES = [
    (1000, 5000),
    (5000, 50000),
    (50000, 500000)
]


def get_all_kmers(k):
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=k)]


ALL_KMERS = get_all_kmers(K)
KMER_INDEX = {kmer: i for i, kmer in enumerate(ALL_KMERS)}


def sequence_to_vector(sequence):
    sequence = sequence.upper()
    vector = np.zeros(len(ALL_KMERS))
    total = 0
    for i in range(len(sequence) - K + 1):
        kmer = sequence[i:i+K]
        if 'N' not in kmer and kmer in KMER_INDEX:
            vector[KMER_INDEX[kmer]] += 1
            total += 1
    if total > 0:
        vector = vector / total
    return vector


def fragment_sequence(sequence, min_len, max_len):
    seq_len = len(sequence)
    if seq_len < min_len:
        return None
    frag_len = random.randint(min_len, min(max_len, seq_len))
    start = random.randint(0, seq_len - frag_len)
    return sequence[start:start + frag_len]


def load_plasmids(fasta_path):
    """Carga todos los plásmidos de PLSDB 2025 sin ningún filtro de tamaño."""
    log.info("Cargando plasmidos de {}".format(fasta_path))
    ids = []
    sequences = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        ids.append(rec.id)
        sequences.append(str(rec.seq))
    log.info("Plasmidos cargados: {}".format(len(sequences)))
    return ids, sequences


def load_chromosomes(chromosomes_dir, frags_per_file=150):
    log.info("Cargando cromosomas de {}".format(chromosomes_dir))
    ids = []
    fragments = []
    files = list(Path(chromosomes_dir).glob("*.fna.gz"))
    log.info("Archivos encontrados: {}".format(len(files)))

    frag_counter = 0
    for f in files:
        try:
            import gzip
            with gzip.open(f, 'rt') as handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    if "plasmid" in rec.description.lower():
                        continue
                    seq = str(rec.seq)
                    for min_len, max_len in FRAGMENT_SIZES:
                        for i in range(frags_per_file):
                            frag = fragment_sequence(seq, min_len, max_len)
                            if frag:
                                frag_id = "{}_frag_{}_{}".format(rec.id, min_len, frag_counter)
                                ids.append(frag_id)
                                fragments.append(frag)
                                frag_counter += 1
        except Exception as e:
            log.warning("Error procesando {}: {}".format(f, e))

    log.info("Fragmentos cromosomicos generados: {}".format(len(fragments)))
    return ids, fragments


def load_phages(fasta_path):
    """Carga secuencias de fagos de RefSeq."""
    log.info("Cargando fagos de {}".format(fasta_path))
    ids = []
    sequences = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        ids.append(rec.id)
        sequences.append(str(rec.seq))
    log.info("Fagos cargados: {}".format(len(sequences)))
    return ids, sequences


def compute_features(sequences, label, desc=""):
    log.info("Calculando features para {} {} secuencias...".format(len(sequences), desc))
    X = np.array([sequence_to_vector(seq) for seq in sequences])
    y = np.array([label] * len(sequences))
    return X, y


def train_model(X, y, all_ids, all_seqs):
    log.info("Dividiendo datos en train/test (80/20)...")
    X_train, X_test, y_train, y_test, ids_train, ids_test, seqs_train, seqs_test = train_test_split(
        X, y, all_ids, all_seqs, test_size=0.2, random_state=RANDOM_SEED, stratify=y
    )

    log.info("Entrenando Random Forest...")
    log.info("Train: {} | Test: {}".format(len(X_train), len(X_test)))

    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        min_samples_split=5,
        min_samples_leaf=2,
        n_jobs=N_JOBS,
        random_state=RANDOM_SEED,
        class_weight='balanced'
    )
    clf.fit(X_train, y_train)

    # Evaluación en test set completo
    y_prob_full = clf.predict_proba(X_test)[:, 1]
    auc_full = roc_auc_score(y_test, y_prob_full)
    log.info("AUC en test set completo: {:.4f}".format(auc_full))
    log.info("\n{}".format(classification_report(y_test, clf.predict(X_test),
             target_names=['chromosome', 'plasmid'])))

    log.info("Guardando test set...")
    test_set = {
        'X_test': X_test,
        'y_test': y_test,
        'y_prob_rf': y_prob_full,
        'ids_test': ids_test
    }
    with open(TEST_SET_OUTPUT, 'wb') as f:
        pickle.dump(test_set, f)
    log.info("Test set guardado en {}".format(TEST_SET_OUTPUT))

    log.info("Guardando test set FASTA...")
    with open(TEST_FASTA_OUTPUT, 'w') as f:
        for seq_id, seq in zip(ids_test, seqs_test):
            f.write(">{}\n{}\n".format(seq_id, seq))
    log.info("Test set FASTA guardado en {}".format(TEST_FASTA_OUTPUT))

    return clf, auc_full


def main():
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    MODEL_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    plasmid_ids, plasmids = load_plasmids(PLASMIDS_FASTA)
    chrom_ids, chromosomes = load_chromosomes(CHROMOSOMES_DIR, frags_per_file=150)
    phage_ids, phages = load_phages(PHAGES_FASTA)

    # Combinar cromosomas y fagos como negativos
    chrom_ids = chrom_ids + phage_ids
    chromosomes = chromosomes + phages
    log.info("Total negativos (cromosomas + fagos): {}".format(len(chromosomes)))

    n = min(len(plasmids), len(chromosomes))
    log.info("Balanceando: {} plasmidos vs {} negativos".format(n, n))

    plas_idx = random.sample(range(len(plasmids)), n)
    chrom_idx = random.sample(range(len(chromosomes)), n)

    plasmids = [plasmids[i] for i in plas_idx]
    plasmid_ids = [plasmid_ids[i] for i in plas_idx]
    chromosomes = [chromosomes[i] for i in chrom_idx]
    chrom_ids = [chrom_ids[i] for i in chrom_idx]

    X_plas, y_plas = compute_features(plasmids, label=1, desc="plasmidos")
    X_chrom, y_chrom = compute_features(chromosomes, label=0, desc="negativos")

    X = np.vstack([X_plas, X_chrom])
    y = np.concatenate([y_plas, y_chrom])
    all_ids = plasmid_ids + chrom_ids
    all_seqs = plasmids + chromosomes

    log.info("Dataset total: {} secuencias, {} features".format(X.shape[0], X.shape[1]))

    clf, auc = train_model(X, y, all_ids, all_seqs)

    with open(MODEL_OUTPUT, 'wb') as f:
        pickle.dump(clf, f)
    log.info("Modelo guardado en {}".format(MODEL_OUTPUT))

    with open(REPORT_OUTPUT, 'w') as f:
        f.write("AUC en test set completo: {:.4f}\n".format(auc))
    log.info("Reporte guardado en {}".format(REPORT_OUTPUT))


if __name__ == "__main__":
    main()