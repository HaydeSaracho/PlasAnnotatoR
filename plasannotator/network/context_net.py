"""
Layer 3: Taxonomic contextual network of plasmids
Shows taxonomic tree: Genus -> Plasmid
AMR, BGC, CAZy and unclassified plasmids are shown only in the report table
"""

import pandas as pd
from pathlib import Path
from pyvis.network import Network


def build_network(annotation_df, output_dir, config_path="config.yaml"):
    """
    Builds a taxonomic network of plasmids.
    Structure: Genus -> Plasmid (genus level only)
    Unclassified plasmids: shown only in the report table.
    Plasmid-genus distance based on PLSDB identity %.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    net = Network(
        height="750px",
        width="100%",
        bgcolor="#1a1a2e",
        font_color="white",
        notebook=False
    )
    net.barnes_hut(
        gravity=-8000,
        central_gravity=0.3,
        spring_length=120,
        spring_strength=0.08,
        damping=0.9
    )

    plasmids = annotation_df[annotation_df["ensemble_label"] == "plasmid"]

    if len(plasmids) == 0:
        print("[Network] No plasmids to visualize.")
        return None

    # Split with and without taxonomy
    if "TAXONOMY_genus" in plasmids.columns:
        with_tax = plasmids[plasmids["TAXONOMY_genus"].notna() &
                            (plasmids["TAXONOMY_genus"].astype(str).str.strip() != "") &
                            (plasmids["TAXONOMY_genus"].astype(str).str.strip() != "nan")]
    else:
        with_tax = pd.DataFrame()

    without_tax = plasmids[~plasmids.index.isin(with_tax.index)]

    print("[Network] Building network with {} plasmids ({} with taxonomy, {} without)...".format(
        len(plasmids), len(with_tax), len(without_tax)))

    if not without_tax.empty:
        print("[Network] {} plasmids without taxonomy (see report table)".format(len(without_tax)))

    if with_tax.empty:
        print("[Network] No plasmids with taxonomy to visualize.")
        return None

    added_nodes = set()

    def safe_add_node(node_id, **kwargs):
        if node_id not in added_nodes:
            net.add_node(node_id, **kwargs)
            added_nodes.add(node_id)

    # Fixed sizes
    SIZE_GENUS = 20
    SIZE_PLASMID = 10

    # Build network: Genus -> Plasmid
    for _, row in with_tax.iterrows():
        contig_id = row["contig_id"]
        score = round(row.get("ensemble_score", 0), 3)
        genus = str(row.get("TAXONOMY_genus", "")).strip()
        pident = row.get("plsdb_pident", None)

        # Plasmid-genus distance based on PLSDB identity %
        # Higher identity = shorter distance
        if pd.notna(pident):
            edge_length = max(30, int(200 - pident * 1.5))
        else:
            edge_length = 100

        # Tooltip
        amr = str(row.get("card_annotation", ""))[:40] if pd.notna(row.get("card_annotation")) else "-"
        bgc = str(row.get("mibig_annotation", ""))[:40] if pd.notna(row.get("mibig_annotation")) else "-"
        cazy = str(row.get("cazy_annotation", ""))[:40] if pd.notna(row.get("cazy_annotation")) else "-"
        pident_str = "{}%".format(round(pident, 1)) if pd.notna(pident) else "-"
        tooltip = "Contig: {}\nEnsemble score: {}\nPLSDB identity: {}\nGenus: {}\nAMR: {}\nBGC: {}\nCAZy: {}".format(
            contig_id, score, pident_str, genus, amr, bgc, cazy)

        # Plasmid node
        safe_add_node(
            contig_id,
            label="",
            title=tooltip,
            color="#e94560",
            size=SIZE_PLASMID,
            shape="dot"
        )

        # Genus node
        if genus and genus != "nan":
            genus_node = "G:{}".format(genus)
            safe_add_node(
                genus_node,
                label=genus,
                title="Genus: {}".format(genus),
                color="#ffb86c",
                size=SIZE_GENUS,
                shape="dot"
            )
            net.add_edge(genus_node, contig_id,
                        color="#ffb86c55", width=1,
                        length=edge_length)

    output_file = output_dir / "plasmid_network.html"
    net.save_graph(str(output_file))
    print("[Network] Network saved to {}".format(output_file))
    return output_file


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))

    test_df = pd.DataFrame({
        "contig_id": ["NODE_5_plasmid", "NODE_42_plasmid", "NODE_99_plasmid",
                      "NODE_7_plasmid", "NODE_8_plasmid"],
        "ensemble_score": [0.98, 0.91, 0.85, 0.77, 0.65],
        "ensemble_label": ["plasmid", "plasmid", "plasmid", "plasmid", "plasmid"],
        "card_annotation": ["tetA tetracycline efflux pump", None, "blaTEM-1", None, None],
        "mibig_annotation": [None, "erythromycin BGC", None, None, None],
        "cazy_annotation": ["GH163", None, "GT4", None, None],
        "TAXONOMY_order": ["Enterobacterales", "Enterobacterales", "Burkholderiales", None, None],
        "TAXONOMY_family": ["Enterobacteriaceae", "Enterobacteriaceae", "Burkholderiaceae", None, None],
        "TAXONOMY_genus": ["Escherichia", "Klebsiella", "Burkholderia", None, None],
        "plsdb_subject": ["NZ_CP012345.1", "NZ_CP098765.1", "NZ_CP111111.1", None, None],
        "plsdb_pident": [98.5, 91.2, 87.3, None, None]
    })

    output_file = build_network(
        annotation_df=test_df,
        output_dir="/tmp/network_test",
        config_path="config.yaml"
    )

    import subprocess
    subprocess.Popen(["xdg-open", str(output_file)])
    print("Network saved to: {}".format(output_file))