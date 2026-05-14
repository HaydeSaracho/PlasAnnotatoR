"""
Constructor de redes de secuencias plasmídicas usando ANI vía Minimap2.
"""

import csv
import subprocess
from pathlib import Path

import networkx as nx
from Bio import SeqIO

from plasannotator.config import (
    ANI_THRESHOLD,
    MINIMAP2_THREADS,
    logger,
)
from plasannotator.db.manager import get_index_path, INDEX_DIR
from plasannotator.network.export import export_graphml, export_json


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

METADATA_PATH = Path("/mnt/d/PlasAnnotatoR/ALL.plasmid_list.download.tsv")


def _load_metadata() -> dict[str, dict]:
    """
    Carga la tabla de metadata de plasmidios.
    Retorna dict {plasmid_id: {topology, completeness, size_bp, gc,
                               host, mob_type, mobility,
                               primary_cluster, secondary_cluster}}
    """
    metadata = {}
    if not METADATA_PATH.exists():
        logger.warning(f"Metadata no encontrada: {METADATA_PATH}")
        return metadata

    with open(METADATA_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row["Plasmid_ID"].strip()
            metadata[pid] = {
                "topology":          row.get("Topology", "-"),
                "completeness":      row.get("Completeness", "-"),
                "size_bp":           row.get("Size (bp)", "-"),
                "gc":                row.get("GC", "-"),
                "host":              row.get("Host", "-"),
                "mob_type":          row.get("MOB_type(s)", "-"),
                "mobility":          row.get("Predicted_Mobility", "-"),
                "primary_cluster":   row.get("Primary_Cluster_ID", "-"),
                "secondary_cluster": row.get("Secondary_Cluster_ID", "-"),
            }
    logger.info(f"Metadata cargada: {len(metadata):,} entradas")
    return metadata


# ---------------------------------------------------------------------------
# Carga de predicciones
# ---------------------------------------------------------------------------

def _load_plasmid_contigs(tsv_path: Path) -> list[dict]:
    """
    Lee el TSV de predicciones y retorna solo los contigs
    clasificados como plasmidio.
    """
    plasmids = []
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["prediction"] == "plasmid":
                plasmids.append({
                    "id": row["contig_id"],
                    "length_bp": int(row["length_bp"]),
                    "confidence": float(row["confidence"]),
                    "db": row["db"],
                })
    return plasmids


def _extract_plasmid_fasta(
    original_fasta: Path,
    plasmid_ids: set[str],
    out_path: Path,
) -> Path:
    """
    Extrae del FASTA original solo los contigs predichos como plasmidio.
    """
    records = [
        r for r in SeqIO.parse(str(original_fasta), "fasta")
        if r.id in plasmid_ids
    ]
    SeqIO.write(records, str(out_path), "fasta")
    logger.info(f"FASTA de plasmidios extraído: {len(records)} contigs → {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Alineamiento
# ---------------------------------------------------------------------------

def _run_all_vs_all_minimap2(
    query_fasta: Path,
    index_path: Path,
    threads: int,
) -> list[tuple[str, str, float]]:
    """
    Alinea los contigs plasmídicos contra la DB de referencia
    y también entre sí (all-vs-all) usando Minimap2.
    Retorna lista de (query_id, target_id, ani_score).
    """
    edges = []

    # Query vs referencia DB
    cmd = [
        "minimap2",
        "-c",
        "--secondary=no",
        "-t", str(threads),
        str(index_path),
        str(query_fasta),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        qid = parts[0]
        tid = parts[5]
        matches = int(parts[9])
        aln_len = int(parts[10])
        qlen = int(parts[1])
        if aln_len > 0:
            ani = (matches / aln_len) * min(aln_len / qlen, 1.0)
            edges.append((qid, tid, round(ani, 4)))

    # Query vs query (all-vs-all entre plasmidios del usuario)
    cmd_self = [
        "minimap2",
        "-c",
        "--secondary=no",
        "-t", str(threads),
        str(query_fasta),
        str(query_fasta),
    ]
    result_self = subprocess.run(cmd_self, capture_output=True, text=True)

    for line in result_self.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        qid = parts[0]
        tid = parts[5]
        if qid == tid:
            continue
        matches = int(parts[9])
        aln_len = int(parts[10])
        qlen = int(parts[1])
        if aln_len > 0:
            ani = (matches / aln_len) * min(aln_len / qlen, 1.0)
            edges.append((qid, tid, round(ani, 4)))

    return edges


# ---------------------------------------------------------------------------
# Construcción del grafo
# ---------------------------------------------------------------------------

def _build_graph(
    plasmids: list[dict],
    edges: list[tuple[str, str, float]],
    ani_threshold: float,
) -> nx.Graph:
    """
    Construye el grafo networkx con nodos (contigs + referencias)
    y aristas filtradas por umbral ANI.
    """
    G = nx.Graph()

    # Agregar nodos de los contigs del usuario
    for p in plasmids:
        G.add_node(
            p["id"],
            node_type="query",
            length_bp=p["length_bp"],
            confidence=p["confidence"],
            db=p["db"],
        )

    # Agregar aristas que superen el umbral ANI
    ref_nodes = set()
    for qid, tid, ani in edges:
        if ani >= ani_threshold:
            if tid not in G:
                G.add_node(tid, node_type="reference")
                ref_nodes.add(tid)
            G.add_edge(qid, tid, ani=ani, weight=ani)

    logger.info(
        f"Grafo construido: {G.number_of_nodes()} nodos · "
        f"{G.number_of_edges()} aristas · "
        f"{len(ref_nodes)} referencias"
    )
    return G


def _enrich_graph(G: nx.Graph, metadata: dict[str, dict]) -> nx.Graph:
    """
    Agrega atributos de metadata a los nodos de referencia del grafo.
    """
    enriched = 0
    for node in G.nodes():
        if G.nodes[node].get("node_type") == "reference":
            meta = metadata.get(node, {})
            if meta:
                G.nodes[node].update(meta)
                enriched += 1
    logger.info(f"Nodos enriquecidos con metadata: {enriched}")
    return G


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def build_network(
    tsv_path: Path,
    db_name: str,
    output_prefix: str,
    ani_threshold: float,
    threads: int,
    original_fasta: Path = None,
) -> None:
    """
    Pipeline completo de construcción de red plasmídica.

    Args:
        tsv_path        : TSV de predicciones de plasannotator predict
        db_name         : nombre de la DB de referencia
        output_prefix   : prefijo para los archivos de salida
        ani_threshold   : umbral ANI mínimo para conectar nodos
        threads         : hilos para Minimap2
        original_fasta  : FASTA original (necesario para extraer secuencias)
    """
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    # 1. Cargar contigs plasmídicos
    logger.info("[1/4] Cargando predicciones...")
    plasmids = _load_plasmid_contigs(tsv_path)

    if not plasmids:
        print("\n  No se encontraron contigs clasificados como plasmidio.\n")
        return

    logger.info(f"  {len(plasmids)} contigs plasmídicos encontrados.")

    # 2. Extraer FASTA de plasmidios
    logger.info("[2/4] Extrayendo secuencias plasmídicas...")
    plasmid_ids = {p["id"] for p in plasmids}

    if original_fasta is None:
        original_fasta = tsv_path.parent / (
            tsv_path.stem.replace("predictions", "input") + ".fasta"
        )

    plasmid_fasta = output_prefix.parent / "plasmids_query.fasta"

    if original_fasta.exists():
        _extract_plasmid_fasta(original_fasta, plasmid_ids, plasmid_fasta)
    else:
        print(
            f"\n  [AVISO] No se encontró el FASTA original en {original_fasta}.\n"
            f"  Pásalo con --fasta para construir la red correctamente.\n"
        )
        return

    # 3. Alineamiento
    logger.info(f"[3/4] Alineando contra DB '{db_name}' (ANI ≥ {ani_threshold})...")
    index_path = get_index_path(db_name)
    edges = _run_all_vs_all_minimap2(plasmid_fasta, index_path, threads)
    logger.info(f"  {len(edges)} pares de secuencias alineados.")

    # 3.5 Cargar metadata
    logger.info("[3.5/4] Cargando metadata de plasmidios...")
    metadata = _load_metadata()

    # 4. Construir y exportar grafo
    logger.info("[4/4] Construyendo grafo...")
    G = _build_graph(plasmids, edges, ani_threshold)
    G = _enrich_graph(G, metadata)

    graphml_path = export_graphml(G, output_prefix)
    json_path = export_json(G, output_prefix)

    # Resumen
    print(f"\n  Nodos totales    : {G.number_of_nodes():,}")
    print(f"  Nodos query      : {len(plasmids):,}")
    print(f"  Nodos referencia : {G.number_of_nodes() - len(plasmids):,}")
    print(f"  Aristas (ANI≥{ani_threshold}) : {G.number_of_edges():,}")
    print(f"  GraphML          : {graphml_path}")
    print(f"  JSON             : {json_path}\n")