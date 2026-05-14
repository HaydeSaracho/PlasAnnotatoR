"""
Exportación de grafos a GraphML y JSON.
"""

import json
from pathlib import Path

import networkx as nx

from plasannotator.config import logger


def export_graphml(G: nx.Graph, output_prefix: Path) -> Path:
    """
    Exporta el grafo a formato GraphML compatible con Gephi y Cytoscape.
    """
    out_path = Path(f"{output_prefix}.graphml")
    nx.write_graphml(G, str(out_path))
    logger.info(f"GraphML exportado: {out_path}")
    return out_path


def export_json(G: nx.Graph, output_prefix: Path) -> Path:
    """
    Exporta el grafo a JSON con formato nodos/aristas
    compatible con visualización web (D3.js, Cytoscape.js).
    """
    out_path = Path(f"{output_prefix}_network.json")

    data = {
        "nodes": [
            {"id": n, **G.nodes[n]}
            for n in G.nodes
        ],
        "edges": [
            {"source": u, "target": v, **G.edges[u, v]}
            for u, v in G.edges
        ],
        "metadata": {
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_query": sum(
                1 for n in G.nodes if G.nodes[n].get("node_type") == "query"
            ),
            "n_reference": sum(
                1 for n in G.nodes if G.nodes[n].get("node_type") == "reference"
            ),
        },
    }

    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"JSON exportado: {out_path}")
    return out_path