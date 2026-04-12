"""
Extracción de features k-mer a partir de secuencias FASTA.
"""

from itertools import product
from collections import Counter

import numpy as np
from Bio import SeqIO

from plasannotator.config import KMER_SIZE, MIN_CONTIG_LENGTH, logger


def get_kmer_vocab(k: int) -> list[str]:
    """Genera todos los k-mers posibles para un k dado."""
    return ["".join(p) for p in product("ACGT", repeat=k)]


def kmer_frequency(sequence: str, k: int, vocab: list[str]) -> np.ndarray:
    """
    Calcula la frecuencia relativa de cada k-mer en una secuencia.
    Retorna un vector numpy normalizado.
    """
    seq = sequence.upper().replace("N", "")
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
    total = sum(counts.values()) or 1
    return np.array([counts.get(kmer, 0) / total for kmer in vocab], dtype=np.float32)


def extract_features(fasta_path) -> tuple[list[str], list[int], np.ndarray]:
    """
    Lee un FASTA y extrae features k-mer para cada contig.

    Retorna:
        ids      — lista de IDs de contigs
        lengths  — longitud de cada contig en pb
        matrix   — matriz (n_contigs x n_features)
    """
    vocab = get_kmer_vocab(KMER_SIZE)
    ids, lengths, vectors = [], [], []

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq)
        if len(seq) < MIN_CONTIG_LENGTH:
            logger.debug(f"Contig omitido (muy corto): {record.id} ({len(seq)} pb)")
            continue

        ids.append(record.id)
        lengths.append(len(seq))
        vectors.append(kmer_frequency(seq, KMER_SIZE, vocab))

    if not ids:
        raise ValueError(
            f"No se encontraron contigs con longitud >= {MIN_CONTIG_LENGTH} pb."
        )

    logger.info(f"Features extraídas: {len(ids)} contigs · k={KMER_SIZE}")
    return ids, lengths, np.vstack(vectors)