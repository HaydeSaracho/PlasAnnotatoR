"""
DB Manager: concatenación de FASTAs e indexación con Minimap2.
"""

import subprocess
from pathlib import Path

from plasannotator.config import (
    AVAILABLE_DBS,
    DB_DIR,
    INDEX_DIR,
    MINIMAP2_THREADS,
    logger,
)


def list_databases() -> None:
    """Muestra el estado de todas las bases de datos disponibles."""
    print(f"\n{'Nombre':<12} {'Secuencias':>12} {'Estado':<15} {'Descripción'}")
    print("-" * 72)

    for name, meta in AVAILABLE_DBS.items():
        db_path = DB_DIR / meta["dir"]
        idx_path = INDEX_DIR / f"{name}.mmi"

        if not db_path.exists():
            status = "no instalada"
        elif not idx_path.exists():
            status = "sin índice"
        else:
            status = "lista"

        print(
            f"{name:<12} {meta['sequences']:>12,} {status:<15} {meta['description']}"
        )
    print()


def concatenate_db(name: str) -> Path:
    """
    Concatena todos los FASTAs individuales de una DB en un solo archivo.
    Retorna la ruta del FASTA concatenado.
    """
    if name not in AVAILABLE_DBS:
        raise ValueError(f"Base de datos '{name}' no reconocida.")

    db_path = DB_DIR / AVAILABLE_DBS[name]["dir"]
    if not db_path.exists():
        raise FileNotFoundError(
            f"Carpeta de DB no encontrada: {db_path}\n"
            f"Asegúrate de que los FASTAs estén en data/databases/{name}/"
        )

    out_fasta = INDEX_DIR / f"{name}.fasta"

    if out_fasta.exists():
        logger.info(f"FASTA concatenado ya existe: {out_fasta}")
        return out_fasta

    logger.info(f"Concatenando FASTAs de '{name}'...")
    fastas = sorted(db_path.glob("*.fasta")) + sorted(db_path.glob("*.fa")) + sorted(db_path.glob("*.fna"))

    if not fastas:
        raise FileNotFoundError(f"No se encontraron archivos FASTA en {db_path}")

    with open(out_fasta, "wb") as out:
        for fasta in fastas:
            with open(fasta, "rb") as f:
                out.write(f.read())

    logger.info(f"Concatenados {len(fastas):,} archivos → {out_fasta}")
    return out_fasta


def build_index(name: str) -> Path:
    """
    Construye el índice Minimap2 (.mmi) para una DB.
    Retorna la ruta del índice.
    """
    idx_path = INDEX_DIR / f"{name}.mmi"

    if idx_path.exists():
        logger.info(f"Índice ya existe: {idx_path}")
        return idx_path

    fasta_path = concatenate_db(name)

    logger.info(f"Construyendo índice Minimap2 para '{name}'...")
    cmd = [
        "minimap2",
        "-d", str(idx_path),
        "-t", str(MINIMAP2_THREADS),
        str(fasta_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"Error al construir índice para '{name}':\n{result.stderr}"
        )

    logger.info(f"Índice construido: {idx_path}")
    return idx_path


def index_database(name: str) -> None:
    """
    Punto de entrada principal: concatena e indexa una DB.
    Equivale al comando: plasannotator db index <nombre>
    """
    try:
        build_index(name)
        print(f"\n[OK] '{name}' indexada y lista para usar.\n")
    except (ValueError, FileNotFoundError, RuntimeError) as e:
        print(f"\n[ERROR] {e}\n")


def get_index_path(name: str) -> Path:
    """
    Retorna la ruta del índice de una DB y verifica que exista.
    Usado por el predictor antes de correr alineamientos.
    """
    if name not in AVAILABLE_DBS:
        raise ValueError(f"Base de datos '{name}' no reconocida.")

    idx_path = INDEX_DIR / f"{name}.mmi"

    if not idx_path.exists():
        raise FileNotFoundError(
            f"Índice no encontrado para '{name}'.\n"
            f"Ejecuta primero: plasannotator db index {name}"
        )

    return idx_path