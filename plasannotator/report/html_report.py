"""
Layer 4: Integrated HTML report for PlasAnnotatoR
"""

import pandas as pd
import shutil
from pathlib import Path
from datetime import datetime


def generate_report(annotation_df, network_html_path, output_dir, input_fasta=""):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plasmids = annotation_df[annotation_df["ensemble_label"] == "plasmid"]
    n_total = len(annotation_df)
    n_plasmids = len(plasmids)
    n_chrom = len(annotation_df[annotation_df["ensemble_label"] == "chromosome"])

    def count(col): return int(plasmids[col].notna().sum()) if col in plasmids.columns else 0
    def pct(n): return round(100 * n / n_plasmids, 1) if n_plasmids > 0 else 0

    # Classification
    # Core modules
    n_replication = count("replication_gene")
    n_conjugation = count("conjugation_gene")
    n_mobility = count("dna_mobility_gene")
    n_transfer = int(plasmids[
        [c for c in ["conjugation_gene", "dna_mobility_gene"] if c in plasmids.columns]
    ].notna().any(axis=1).sum()) if n_plasmids > 0 else 0
    n_stability = count("stability_gene")

    # Accessory
    n_amr = count("card_subject")
    n_bgc = count("mibig_subject")
    n_cazy = count("cazy_subject")
    n_virulence = count("virulence_defense_gene")
    n_metal = count("metal_biocide_gene")
    n_stress = count("stress_response_gene")

    # Taxonomy
    n_tax = count("TAXONOMY_genus")

    # Network
    network_html = ""
    if network_html_path and Path(network_html_path).exists():
        network_filename = Path(network_html_path).name
        dest = output_dir / network_filename
        if Path(network_html_path).resolve() != dest.resolve():
            shutil.copy(str(network_html_path), str(dest))
        network_html = '<iframe src="{}" width="100%" height="580px" frameborder="0" style="border-radius:6px; border: 1px solid #dee2e6;"></iframe>'.format(network_filename)

    # Table rows
    table_rows = ""
    for i, (_, row) in enumerate(plasmids.iterrows()):
        score = round(row.get("ensemble_score", 0), 3)

        def badge(col_subject, col_annot, tag_class, tag_label):
            val = row.get(col_annot, "")
            val = str(val)[:50] if pd.notna(val) else "-"
            if pd.notna(row.get(col_subject)):
                return '<span class="tag tag-{}">{}</span> {}'.format(tag_class, tag_label, val)
            return val

        def plasann_badge(gene_col, product_col, tag_class, tag_label):
            gene = row.get(gene_col, "")
            product = row.get(product_col, "")
            if pd.notna(gene) and str(gene) != "nan":
                prod_str = str(product)[:40] if pd.notna(product) else ""
                return '<span class="tag tag-{}">{}</span> {}'.format(tag_class, str(gene), prod_str)
            return "-"

        amr = badge("card_subject", "card_annotation", "amr", "AMR")
        bgc = badge("mibig_subject", "mibig_annotation", "bgc", "BGC")
        cazy = badge("cazy_subject", "cazy_annotation", "cazy", "CAZy")
        replication = plasann_badge("replication_gene", "replication_product", "replication", "REP")
        transfer_conj = plasann_badge("conjugation_gene", "conjugation_product", "transfer", "CONJ")
        transfer_mob = plasann_badge("dna_mobility_gene", "dna_mobility_product", "transfer", "MOB")
        transfer = transfer_conj if transfer_conj != "-" else transfer_mob
        stability = plasann_badge("stability_gene", "stability_product", "stability", "TA")
        virulence = plasann_badge("virulence_defense_gene", "virulence_defense_product", "virulence", "VIR")
        metal = plasann_badge("metal_biocide_gene", "metal_biocide_product", "metal", "MET")
        stress = plasann_badge("stress_response_gene", "stress_response_product", "stress", "STR")
        genus = row.get("TAXONOMY_genus", "")
        genus = "<em>{}</em>".format(str(genus)) if pd.notna(genus) else "-"

        has_annotation = any([
            pd.notna(row.get("card_subject")),
            pd.notna(row.get("mibig_subject")),
            pd.notna(row.get("cazy_subject")),
            pd.notna(row.get("replication_gene")),
            pd.notna(row.get("conjugation_gene")),
            pd.notna(row.get("dna_mobility_gene")),
            pd.notna(row.get("stability_gene")),
        ])
        row_class = "annotated" if has_annotation else ""
        row_bg = 'style="background:#f8f9fa;"' if i % 2 == 0 else ""

        table_rows += """
        <tr class="{}" {}>
            <td style="font-family:monospace;font-size:0.82em;">{}</td>
            <td style="text-align:center;"><span class="score-badge">{}</span></td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
        </tr>""".format(row_class, row_bg, row["contig_id"], score,
                        replication, transfer, stability,
                        amr, bgc, cazy, virulence, metal, stress, genus)

    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>PlasAnnotatoR Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:'Segoe UI',Arial,sans-serif; background:#f4f6f9; color:#2c3e50; font-size:14px; }}
        header {{ background:#fff; padding:24px 48px; border-bottom:3px solid #c2185b;
                  display:flex; align-items:center; justify-content:space-between;
                  box-shadow:0 2px 8px rgba(0,0,0,0.08); }}
        .header-left h1 {{ color:#c2185b; font-size:1.8em; font-weight:700; }}
        .header-left p {{ color:#6c757d; margin-top:4px; font-size:0.88em; }}
        .header-right {{ text-align:right; color:#6c757d; font-size:0.82em; }}
        .header-right strong {{ color:#c2185b; }}
        .container {{ max-width:1600px; margin:0 auto; padding:32px 48px; }}

        /* PANEL SECTIONS */
        .panel-section {{ margin-bottom:12px; }}
        .panel-label {{ font-size:0.72em; font-weight:700; color:#6c757d;
                        text-transform:uppercase; letter-spacing:0.8px;
                        margin-bottom:8px; padding-left:4px;
                        border-left:3px solid #c2185b; padding-left:8px; }}
        .stats-grid {{ display:grid; gap:12px; margin-bottom:4px; }}
        .grid-3 {{ grid-template-columns: repeat(3,1fr); }}
        .grid-4 {{ grid-template-columns: repeat(4,1fr); }}
        .grid-6 {{ grid-template-columns: repeat(6,1fr); }}
        .grid-1 {{ grid-template-columns: repeat(1,1fr); width:33%; }}

        .stat-card {{ background:#fff; border-radius:8px; padding:16px 12px;
                      text-align:center; border:1px solid #e0e6ed;
                      box-shadow:0 1px 4px rgba(0,0,0,0.06); }}
        .stat-card .number {{ font-size:2em; font-weight:700; line-height:1; }}
        .stat-card .pct {{ font-size:0.75em; color:#6c757d; margin-top:4px; }}
        .stat-card .label {{ color:#6c757d; margin-top:6px; font-size:0.75em;
                             font-weight:500; text-transform:uppercase; letter-spacing:0.3px; }}

        /* Classification */
        .stat-card.total .number {{ color:#2c3e50; }}
        .stat-card.plasmid {{ border-top:3px solid #e63946; }}
        .stat-card.plasmid .number {{ color:#e63946; }}
        .stat-card.chrom {{ border-top:3px solid #6c757d; }}
        .stat-card.chrom .number {{ color:#6c757d; }}

        /* Core modules */
        .stat-card.replication {{ border-top:3px solid #1565c0; }}
        .stat-card.replication .number {{ color:#1565c0; }}
        .stat-card.transfer {{ border-top:3px solid #6a1b9a; }}
        .stat-card.transfer .number {{ color:#6a1b9a; }}
        .stat-card.stability {{ border-top:3px solid #00695c; }}
        .stat-card.stability .number {{ color:#00695c; }}

        /* Accessory */
        .stat-card.amr {{ border-top:3px solid #e76f51; }}
        .stat-card.amr .number {{ color:#e76f51; }}
        .stat-card.bgc {{ border-top:3px solid #2a9d8f; }}
        .stat-card.bgc .number {{ color:#2a9d8f; }}
        .stat-card.cazy {{ border-top:3px solid #457b9d; }}
        .stat-card.cazy .number {{ color:#457b9d; }}
        .stat-card.virulence {{ border-top:3px solid #b71c1c; }}
        .stat-card.virulence .number {{ color:#b71c1c; }}
        .stat-card.metal {{ border-top:3px solid #4527a0; }}
        .stat-card.metal .number {{ color:#4527a0; }}
        .stat-card.stress {{ border-top:3px solid #e65100; }}
        .stat-card.stress .number {{ color:#e65100; }}

        /* Taxonomy */
        .stat-card.tax {{ border-top:3px solid #9b72cf; }}
        .stat-card.tax .number {{ color:#9b72cf; }}

        /* SECTIONS */
        .section {{ background:#fff; border-radius:8px; padding:24px;
                    margin-bottom:24px; border:1px solid #e0e6ed;
                    box-shadow:0 1px 4px rgba(0,0,0,0.06); margin-top:20px; }}
        .section-header {{ display:flex; align-items:center; margin-bottom:16px;
                           padding-bottom:12px; border-bottom:1px solid #e0e6ed; }}
        .section-header h2 {{ color:#c2185b; font-size:1.05em; font-weight:600;
                              text-transform:uppercase; letter-spacing:0.5px; }}

        /* TABLE */
        table {{ width:100%; border-collapse:collapse; font-size:0.82em; }}
        th {{ background:#c2185b; color:#fff; padding:10px 12px; text-align:left;
              font-weight:600; font-size:0.80em; text-transform:uppercase; letter-spacing:0.4px; }}
        th.core {{ background:#37474f; }}
        th.accessory {{ background:#455a64; }}
        td {{ padding:8px 12px; border-bottom:1px solid #e0e6ed; vertical-align:middle; }}
        tr:hover td {{ background:#fce4ec !important; }}
        tr.annotated td:first-child {{ border-left:3px solid #e63946; }}

        /* BADGES */
        .score-badge {{ background:#fce4ec; color:#c2185b; padding:2px 8px;
                        border-radius:12px; font-size:0.85em; font-weight:600;
                        border:1px solid #f48fb1; }}
        .tag {{ display:inline-block; padding:1px 6px; border-radius:4px;
                font-size:0.72em; font-weight:700; margin-right:4px; vertical-align:middle; }}
        .tag-amr {{ background:#fde8e4; color:#e76f51; border:1px solid #f4b8a8; }}
        .tag-bgc {{ background:#d8f3f0; color:#2a9d8f; border:1px solid #a8ddd8; }}
        .tag-cazy {{ background:#dceaf5; color:#457b9d; border:1px solid #a8c8e8; }}
        .tag-replication {{ background:#e3f2fd; color:#1565c0; border:1px solid #bbdefb; }}
        .tag-transfer {{ background:#ede7f6; color:#6a1b9a; border:1px solid #d1c4e9; }}
        .tag-stability {{ background:#e0f2f1; color:#00695c; border:1px solid #b2dfdb; }}
        .tag-virulence {{ background:#ffebee; color:#b71c1c; border:1px solid #ffcdd2; }}
        .tag-metal {{ background:#ede7f6; color:#4527a0; border:1px solid #d1c4e9; }}
        .tag-stress {{ background:#fff3e0; color:#e65100; border:1px solid #ffe0b2; }}

        footer {{ text-align:center; padding:20px; color:#6c757d; font-size:0.80em;
                  border-top:1px solid #e0e6ed; margin-top:8px; background:#fff; }}
        footer strong {{ color:#c2185b; }}
    </style>
</head>
<body>
    <header>
        <div class="header-left">
            <h1>&#129516; PlasAnnotatoR</h1>
            <p>Plasmid classification, annotation and taxonomy report</p>
        </div>
        <div class="header-right">
            <div><strong>Input:</strong> {input}</div>
            <div><strong>Date:</strong> {date}</div>
        </div>
    </header>
    <div class="container">

        <!-- CLASSIFICATION -->
        <div class="panel-section">
            <div class="panel-label">Classification</div>
            <div class="stats-grid grid-3">
                <div class="stat-card total">
                    <div class="number">{n_total}</div>
                    <div class="label">Total contigs</div>
                </div>
                <div class="stat-card plasmid">
                    <div class="number">{n_plasmids}</div>
                    <div class="label">Plasmids detected</div>
                </div>
                <div class="stat-card chrom">
                    <div class="number">{n_chrom}</div>
                    <div class="label">Chromosomal</div>
                </div>
            </div>
        </div>

        <!-- CORE MODULES -->
        <div class="panel-section">
            <div class="panel-label">Plasmid Core Modules</div>
            <div class="stats-grid grid-3">
                <div class="stat-card replication">
                    <div class="number">{n_replication}</div>
                    <div class="pct">{pct_replication}% of plasmids</div>
                    <div class="label">Replication</div>
                </div>
                <div class="stat-card transfer">
                    <div class="number">{n_transfer}</div>
                    <div class="pct">{pct_transfer}% of plasmids</div>
                    <div class="label">Transfer</div>
                </div>
                <div class="stat-card stability">
                    <div class="number">{n_stability}</div>
                    <div class="pct">{pct_stability}% of plasmids</div>
                    <div class="label">Stability</div>
                </div>
            </div>
        </div>

        <!-- ACCESSORY GENES -->
        <div class="panel-section">
            <div class="panel-label">Accessory Genes</div>
            <div class="stats-grid grid-6">
                <div class="stat-card amr">
                    <div class="number">{n_amr}</div>
                    <div class="pct">{pct_amr}% of plasmids</div>
                    <div class="label">AMR genes</div>
                </div>
                <div class="stat-card bgc">
                    <div class="number">{n_bgc}</div>
                    <div class="pct">{pct_bgc}% of plasmids</div>
                    <div class="label">BGCs</div>
                </div>
                <div class="stat-card cazy">
                    <div class="number">{n_cazy}</div>
                    <div class="pct">{pct_cazy}% of plasmids</div>
                    <div class="label">CAZymes</div>
                </div>
                <div class="stat-card virulence">
                    <div class="number">{n_virulence}</div>
                    <div class="pct">{pct_virulence}% of plasmids</div>
                    <div class="label">Virulence</div>
                </div>
                <div class="stat-card metal">
                    <div class="number">{n_metal}</div>
                    <div class="pct">{pct_metal}% of plasmids</div>
                    <div class="label">Metal/Biocide</div>
                </div>
                <div class="stat-card stress">
                    <div class="number">{n_stress}</div>
                    <div class="pct">{pct_stress}% of plasmids</div>
                    <div class="label">Stress Response</div>
                </div>
            </div>
        </div>

        <!-- TAXONOMY -->
        <div class="panel-section">
            <div class="panel-label">Assigned Taxonomy</div>
            <div class="stats-grid grid-1">
                <div class="stat-card tax">
                    <div class="number">{n_tax}</div>
                    <div class="pct">{pct_tax}% of plasmids</div>
                    <div class="label">PLSDB taxonomy</div>
                </div>
            </div>
        </div>

        <!-- NETWORK -->
        <div class="section">
            <div class="section-header"><h2>Contextual Network</h2></div>
            {network_html}
        </div>

        <!-- TABLE -->
        <div class="section">
            <div class="section-header"><h2>Detected Plasmids</h2></div>
            <table>
                <thead>
                    <tr>
                        <th rowspan="2">Contig ID</th>
                        <th rowspan="2" style="text-align:center;">Score</th>
                        <th colspan="3" class="core" style="text-align:center;">Core Modules</th>
                        <th colspan="6" class="accessory" style="text-align:center;">Accessory Genes</th>
                        <th rowspan="2">Genus (PLSDB)</th>
                    </tr>
                    <tr>
                        <th class="core">Replication</th>
                        <th class="core">Transfer</th>
                        <th class="core">Stability</th>
                        <th class="accessory">AMR (CARD)</th>
                        <th class="accessory">BGC (MIBiG)</th>
                        <th class="accessory">CAZyme</th>
                        <th class="accessory">Virulence</th>
                        <th class="accessory">Metal/Biocide</th>
                        <th class="accessory">Stress</th>
                    </tr>
                </thead>
                <tbody>{table_rows}</tbody>
            </table>
        </div>
    </div>
    <footer>
        <strong>PlasAnnotatoR</strong> &nbsp;|&nbsp;
        Ensemble: RF(0.9899) + PLASMe(0.9519) + PlasmidHunter(0.9702) + PlasClass(0.9552)
    </footer>
</body>
</html>""".format(
        input=Path(input_fasta).name if input_fasta else "N/A",
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_total=n_total, n_plasmids=n_plasmids, n_chrom=n_chrom,
        n_replication=n_replication, pct_replication=pct(n_replication),
        n_transfer=n_transfer, pct_transfer=pct(n_transfer),
        n_stability=n_stability, pct_stability=pct(n_stability),
        n_amr=n_amr, pct_amr=pct(n_amr),
        n_bgc=n_bgc, pct_bgc=pct(n_bgc),
        n_cazy=n_cazy, pct_cazy=pct(n_cazy),
        n_virulence=n_virulence, pct_virulence=pct(n_virulence),
        n_metal=n_metal, pct_metal=pct(n_metal),
        n_stress=n_stress, pct_stress=pct(n_stress),
        n_tax=n_tax, pct_tax=pct(n_tax),
        network_html=network_html,
        table_rows=table_rows
    )

    output_file = output_dir / "plasannotator_report.html"
    with open(output_file, 'w') as f:
        f.write(html)
    print("[Report] Report saved to {}".format(output_file))
    return output_file
