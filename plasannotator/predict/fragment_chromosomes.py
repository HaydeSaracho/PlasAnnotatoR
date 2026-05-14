"""
Fragmentación de cromosomas bacterianos para dataset negativo de entrenamiento.
"""

import gzip
import random
from pathlib import Path
from Bio import SeqIO

INPUT_DIR = Path("/mnt/d/PlasAnnotatoR/data/chromosomes_raw")
OUTPUT = Path("/mnt/d/PlasAnnotatoR/data/chromosomes_fragments.fasta")

MIN_LEN = 1000
MAX_LEN = 500000
FRAGMENTS_PER_GENOME = 10
MAX_N_RATIO = 0.1
SEED = 42

random.seed(SEED)

total = 0
errors = 0

with open(OUTPUT, "w") as out:
    files = sorted(INPUT_DIR.glob("*.fna.gz"))
    print(f"Procesando {len(files)} genomas...")
    for gz_file in files:
        accession = gz_file.stem.replace("_genomic.fna", "")
        try:
            with gzip.open(gz_file, "rt") as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq = str(record.seq)
                    seq_len = len(seq)
                    if seq_len < MIN_LEN:
                        continue
                    for i in range(FRAGMENTS_PER_GENOME):
                        frag_len = random.randint(MIN_LEN, min(MAX_LEN, seq_len))
                        start = random.randint(0, seq_len - frag_len)
                        fragment = seq[start:start + frag_len]
                        if fragment.count("N") / frag_len > MAX_N_RATIO:
                            continue
                        frag_id = f"{accession}_frag{i+1}"
                        out.write(f">{frag_id}\n{fragment}\n")
                        total += 1
        except Exception as e:
            errors += 1

print(f"Fragmentos generados: {total:,}")
print(f"Errores: {errors}")
