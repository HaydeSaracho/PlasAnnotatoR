"""
Configuración global: rutas, parámetros por defecto y logging.
"""

import logging
from pathlib import Path

# --- Rutas base ---
ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data"
DB_DIR = DATA_DIR / "databases"
INDEX_DIR = DATA_DIR / "indexes"
RESULTS_DIR = DATA_DIR / "results"

# Crear directorios si no existen
for _dir in [DATA_DIR, DB_DIR, INDEX_DIR, RESULTS_DIR]:
    _dir.mkdir(parents=True, exist_ok=True)

# --- Parámetros de predicción ---
KMER_SIZE = 5
MIN_CONTIG_LENGTH = 1000
CONFIDENCE_THRESHOLD = 0.7

# --- Parámetros de red ---
ANI_THRESHOLD = 0.95

# --- Parámetros de indexación ---
MINIMAP2_THREADS = 4
MINIMAP2_PRESET = "map-ont"  # map-ont | asm5 | asm10

# --- Bases de datos disponibles ---
AVAILABLE_DBS = {
    "plsdb":   {"dir": "plsdb",   "description": "PLSDB plasmid database",              "sequences": 50554},
    "imgpr":   {"dir": "imgpr",   "description": "IMG/PR integrated plasmid database",   "sequences": 699973},
    "compass": {"dir": "compass", "description": "COMPASS plasmid database",             "sequences": 12084},
    "genbank": {"dir": "genbank", "description": "NCBI GenBank plasmids",                "sequences": 108316},
    "refseq":  {"dir": "refseq",  "description": "NCBI RefSeq plasmids",                "sequences": 86009},
    "ddbj":    {"dir": "ddbj",    "description": "DDBJ plasmid sequences",               "sequences": 7794},
    "embl":    {"dir": "embl",    "description": "EMBL plasmid sequences",               "sequences": 22604},
    "tpa":     {"dir": "tpa",     "description": "Third Party Annotation",               "sequences": 8},
    "kraken2": {"dir": "kraken2", "description": "Kraken2 plasmid library",              "sequences": 905},
    "mmge":    {"dir": "mmge",    "description": "mMGE mobile genetic elements",         "sequences": 92492},
}

# --- Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger("plasannotator")