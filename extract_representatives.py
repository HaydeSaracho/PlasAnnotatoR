"""
Extrae secuencias completas de los FASTAs concatenados
usando los IDs de los representantes nr.
"""
from Bio import SeqIO
from pathlib import Path

INDEX_DIR = Path("/mnt/d/PlasAnnotatoR/data/indexes")
DB_DIR = Path("/mnt/d/PlasAnnotatoR/data/databases")
OUTPUT = INDEX_DIR / "all_plasmids_nr_full.fasta"

# DBs con concat disponible en indexes/
CONCAT_DBS = {
    "plsdb":   INDEX_DIR / "plsdb_concat.fasta",
    "compass": INDEX_DIR / "compass_concat.fasta",
    "embl":    INDEX_DIR / "embl_concat.fasta",
    "refseq":  INDEX_DIR / "refseq_concat.fasta",
    "mmge":    INDEX_DIR / "mmge_concat.fasta",
}

# DBs pequeñas — leer directamente desde databases/
SMALL_DBS = {
    "ddbj":    (DB_DIR / "ddbj",    "*.fa"),
    "kraken2": (DB_DIR / "kraken2", "*.fasta"),
    "tpa":     (DB_DIR / "tpa",     "*.fa"),
}

def load_representatives(nr_fasta):
    ids = set()
    with open(nr_fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line.strip()[1:].split()[0])
    return ids

total = 0
with open(OUTPUT, "w") as out:
    for db, concat_path in CONCAT_DBS.items():
        nr_path = INDEX_DIR / f"{db}_nr.fasta"
        if not nr_path.exists() or not concat_path.exists():
            print(f"  Saltando {db} — archivos no encontrados")
            continue
        reps = load_representatives(nr_path)
        print(f"Extrayendo {db}: {len(reps):,} representantes...")
        count = 0
        for record in SeqIO.parse(str(concat_path), "fasta"):
            if record.id in reps and len(record.seq) >= 1000:
                out.write(f">{record.id}\n{str(record.seq)}\n")
                count += 1
        print(f"  {count:,} secuencias extraídas")
        total += count

    for db, (db_path, pattern) in SMALL_DBS.items():
        nr_path = INDEX_DIR / f"{db}_nr.fasta"
        if not nr_path.exists() or not db_path.exists():
            print(f"  Saltando {db} — archivos no encontrados")
            continue
        reps = load_representatives(nr_path)
        print(f"Extrayendo {db}: {len(reps):,} representantes...")
        count = 0
        for fasta in sorted(db_path.glob(pattern)):
            for record in SeqIO.parse(str(fasta), "fasta"):
                if record.id in reps and len(record.seq) >= 1000:
                    out.write(f">{record.id}\n{str(record.seq)}\n")
                    count += 1
        print(f"  {count:,} secuencias extraídas")
        total += count

print(f"\nTotal: {total:,} secuencias → {OUTPUT}")
