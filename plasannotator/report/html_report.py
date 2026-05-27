"""
Layer 4: Integrated HTML report for PlasAnnotatoR
Generates a complete report with statistics, tables and interactive network
"""

import pandas as pd
import shutil
from pathlib import Path
from datetime import datetime


PLASANN_CATEGORIES = [
    ("conjugation", "Conjugation", "#6d4c41", "#efebe9", "#d7ccc8"),
    ("virulence_defense", "Virulence", "#b71c1c", "#ffebee", "#ffcdd2"),
    ("plasmid_backbone", "Backbone", "#1565c0", "#e3f2fd", "#bbdefb"),
    ("stress_response", "Stress", "#e65100", "#fff3e0", "#ffe0b2"),
    ("dna_mobility", "DNA Mobility", "#558b2f", "#f1f8e9", "#dcedc8"),
    ("metal_biocide", "Metal/Biocide", "#4527a0", "#ede7f6", "#d1c4e9"),
    ("toxin_antitoxin", "Toxin-AT", "#00695c", "#e0f2f1", "#b2dfdb"),
]


def generate_report(annotation_df, network_html_path, output_dir, input_fasta=""):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plasmids = annotation_df[annotation_df["ensemble_label"] == "plasmid"]
    chromosomes = annotation_df[annotation_df["ensemble_label"] == "chromosome"]

    n_total = len(annotation_df)
    n_plasmids = len(plasmids)
    n_chrom = len(chromosomes)
    n_amr = int(plasmids["card_subject"].notna().sum()) if "card_subject" in plasmids.columns else 0
    n_bgc = int(plasmids["mibig_subject"].notna().sum()) if "mibig_subject" in plasmids.columns else 0
    n_cazy = int(plasmids["cazy_subject"].notna().sum()) if "cazy_subject" in plasmids.columns else 0
    n_tax = int(plasmids["TAXONOMY_genus"].notna().sum()) if "TAXONOMY_genus" in plasmids.columns else 0

    pct_amr = round(100 * n_amr / n_plasmids, 1) if n_plasmids > 0 else 0
    pct_bgc = round(100 * n_bgc / n_plasmids, 1) if n_plasmids > 0 else 0
    pct_cazy = round(100 * n_cazy / n_plasmids, 1) if n_plasmids > 0 else 0
    pct_tax = round(100 * n_tax / n_plasmids, 1) if n_plasmids > 0 else 0

    # PlasAnn stats
    plasann_stats = {}
    for label, display, color, bg, border in PLASANN_CATEGORIES:
        col = "{}_gene".format(label)
        n = int(plasmids[col].notna().sum()) if col in plasmids.columns else 0
        pct = round(100 * n / n_plasmids, 1) if n_plasmids > 0 else 0
        plasann_stats[label] = (n, pct, display, color, bg, border)

    # Network
    network_html = ""
    if network_html_path and Path(network_html_path).exists():
        network_filename = Path(network_html_path).name
        dest = output_dir / network_filename
        if Path(network_html_path).resolve() != dest.resolve():
            shutil.copy(str(network_html_path), str(dest))
        network_html = '<iframe src="{}" width="100%" height="580px" frameborder="0" style="border-radius:6px; border: 1px solid #dee2e6;"></iframe>'.format(network_filename)

    # PlasAnn tag CSS
    plasann_tag_css = ""
    for label, display, color, bg, border in PLASANN_CATEGORIES:
        plasann_tag_css += ".tag-{} {{ background: {}; color: {}; border: 1px solid {}; }}\n".format(
            label, bg, color, border)

    # PlasAnn stat cards CSS
    plasann_card_css = ""
    for label, display, color, bg, border in PLASANN_CATEGORIES:
        plasann_card_css += ".stat-card.{} {{ border-top: 3px solid {}; }}\n".format(label, color)
        plasann_card_css += ".stat-card.{} .number {{ color: {}; }}\n".format(label, color)

    # PlasAnn stat cards HTML
    plasann_cards_html = ""
    for label, (n, pct, display, color, bg, border) in plasann_stats.items():
        plasann_cards_html += """
        <div class="stat-card {}">
            <div class="number">{}</div>
            <div class="pct">{}% of plasmids</div>
            <div class="label">{}</div>
        </div>""".format(label, n, pct, display)

    # Table header extra columns
    plasann_th = ""
    for label, display, color, bg, border in PLASANN_CATEGORIES:
        plasann_th += "<th>{}</th>".format(display)

    # Table rows
    table_rows = ""
    for i, (_, row) in enumerate(plasmids.iterrows()):
        amr = row.get("card_annotation", "")
        amr = str(amr)[:50] if pd.notna(amr) else "-"
        bgc = row.get("mibig_annotation", "")
        bgc = str(bgc)[:50] if pd.notna(bgc) else "-"
        cazy = row.get("cazy_annotation", "")
        cazy = str(cazy)[:50] if pd.notna(cazy) else "-"
        genus = row.get("TAXONOMY_genus", "")
        genus = "<em>{}</em>".format(str(genus)) if pd.notna(genus) else "-"
        score = round(row.get("ensemble_score", 0), 3)

        amr_badge = '<span class="tag tag-amr">AMR</span> {}'.format(amr) if pd.notna(row.get("card_subject")) else amr
        bgc_badge = '<span class="tag tag-bgc">BGC</span> {}'.format(bgc) if pd.notna(row.get("mibig_subject")) else bgc
        cazy_badge = '<span class="tag tag-cazy">CAZy</span> {}'.format(cazy) if pd.notna(row.get("cazy_subject")) else cazy

        # PlasAnn columns
        plasann_tds = ""
        for label, display, color, bg, border in PLASANN_CATEGORIES:
            gene_col = "{}_gene".format(label)
            product_col = "{}_product".format(label)
            gene = row.get(gene_col, "")
            product = row.get(product_col, "")
            if pd.notna(gene) and str(gene) != "nan":
                product_str = str(product)[:40] if pd.notna(product) else ""
                plasann_tds += '<td><span class="tag tag-{}">{}</span> {}</td>'.format(
                    label, str(gene), product_str)
            else:
                plasann_tds += "<td>-</td>"

        has_annotation = any([
            pd.notna(row.get("card_subject")),
            pd.notna(row.get("mibig_subject")),
            pd.notna(row.get("cazy_subject")),
            any(pd.notna(row.get("{}_gene".format(l))) for l, *_ in PLASANN_CATEGORIES)
        ])

        row_class = "annotated" if has_annotation else ""
        row_bg = 'style="background:#f8f9fa;"' if i % 2 == 0 else ""

        table_rows += """
        <tr class="{}" {}>
            <td style="font-family:monospace; font-size:0.82em;">{}</td>
            <td style="text-align:center;"><span class="score-badge">{}</span></td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            {}
            <td>{}</td>
        </tr>""".format(row_class, row_bg, row["contig_id"], score,
                        amr_badge, bgc_badge, cazy_badge, plasann_tds, genus)

    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PlasAnnotatoR - Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f4f6f9;
                color: #2c3e50; font-size: 14px; }}
        header {{ background: #ffffff; padding: 24px 48px;
                  border-bottom: 3px solid #c2185b;
                  display: flex; align-items: center; justify-content: space-between;
                  box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
        .header-left h1 {{ color: #c2185b; font-size: 1.8em; font-weight: 700; letter-spacing: -0.5px; }}
        .header-left p {{ color: #6c757d; margin-top: 4px; font-size: 0.88em; }}
        .header-right {{ text-align: right; color: #6c757d; font-size: 0.82em; }}
        .header-right strong {{ color: #c2185b; }}
        .container {{ max-width: 1600px; margin: 0 auto; padding: 32px 48px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(7, 1fr);
                       gap: 16px; margin-bottom: 20px; }}
        .stats-grid-plasann {{ display: grid; grid-template-columns: repeat(7, 1fr);
                               gap: 16px; margin-bottom: 32px; }}
        .stats-section-label {{ font-size: 0.75em; font-weight: 600; color: #6c757d;
                                text-transform: uppercase; letter-spacing: 0.5px;
                                margin-bottom: 8px; }}
        .stat-card {{ background: #ffffff; border-radius: 8px; padding: 18px 12px;
                      text-align: center; border: 1px solid #e0e6ed;
                      box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
        .stat-card .number {{ font-size: 2em; font-weight: 700; line-height: 1; }}
        .stat-card .pct {{ font-size: 0.78em; color: #6c757d; margin-top: 4px; }}
        .stat-card .label {{ color: #6c757d; margin-top: 6px; font-size: 0.78em;
                             font-weight: 500; text-transform: uppercase; letter-spacing: 0.3px; }}
        .stat-card.total .number {{ color: #2c3e50; }}
        .stat-card.plasmid {{ border-top: 3px solid #e63946; }}
        .stat-card.plasmid .number {{ color: #e63946; }}
        .stat-card.chrom {{ border-top: 3px solid #6c757d; }}
        .stat-card.chrom .number {{ color: #6c757d; }}
        .stat-card.amr {{ border-top: 3px solid #e76f51; }}
        .stat-card.amr .number {{ color: #e76f51; }}
        .stat-card.bgc {{ border-top: 3px solid #2a9d8f; }}
        .stat-card.bgc .number {{ color: #2a9d8f; }}
        .stat-card.cazy {{ border-top: 3px solid #457b9d; }}
        .stat-card.cazy .number {{ color: #457b9d; }}
        .stat-card.tax {{ border-top: 3px solid #9b72cf; }}
        .stat-card.tax .number {{ color: #9b72cf; }}
        {}
        .section {{ background: #ffffff; border-radius: 8px; padding: 24px;
                    margin-bottom: 24px; border: 1px solid #e0e6ed;
                    box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
        .section-header {{ display: flex; align-items: center; margin-bottom: 16px;
                           padding-bottom: 12px; border-bottom: 1px solid #e0e6ed; }}
        .section-header h2 {{ color: #c2185b; font-size: 1.05em; font-weight: 600;
                              text-transform: uppercase; letter-spacing: 0.5px; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.84em; }}
        th {{ background: #c2185b; color: #ffffff; padding: 10px 12px;
              text-align: left; font-weight: 600; font-size: 0.82em;
              text-transform: uppercase; letter-spacing: 0.4px; }}
        td {{ padding: 8px 12px; border-bottom: 1px solid #e0e6ed; vertical-align: middle; }}
        tr:hover td {{ background: #fce4ec !important; }}
        tr.annotated td:first-child {{ border-left: 3px solid #e63946; }}
        .score-badge {{ background: #fce4ec; color: #c2185b; padding: 2px 8px;
                        border-radius: 12px; font-size: 0.85em; font-weight: 600;
                        border: 1px solid #f48fb1; }}
        .tag {{ display: inline-block; padding: 1px 6px; border-radius: 4px;
                font-size: 0.72em; font-weight: 700; margin-right: 4px; vertical-align: middle; }}
        .tag-amr {{ background: #fde8e4; color: #e76f51; border: 1px solid #f4b8a8; }}
        .tag-bgc {{ background: #d8f3f0; color: #2a9d8f; border: 1px solid #a8ddd8; }}
        .tag-cazy {{ background: #dceaf5; color: #457b9d; border: 1px solid #a8c8e8; }}
        {}
        .network-container {{ width: 100%; height: 580px; border-radius: 6px; overflow: hidden; }}
        footer {{ text-align: center; padding: 20px; color: #6c757d;
                  font-size: 0.80em; border-top: 1px solid #e0e6ed;
                  margin-top: 8px; background: #ffffff; }}
        footer strong {{ color: #c2185b; }}
    </style>
</head>
<body>
    <header>
        <div class="header-left">
            <h1>&#129516; PlasAnnotatoR</h1>
            <p>Plasmid classification and annotation report</p>
        </div>
        <div class="header-right">
            <div><strong>Input:</strong> {}</div>
            <div><strong>Date:</strong> {}</div>
        </div>
    </header>
    <div class="container">
        <div class="stats-section-label">Classification &amp; Functional Annotation</div>
        <div class="stats-grid">
            <div class="stat-card total">
                <div class="number">{}</div>
                <div class="label">Total contigs</div>
            </div>
            <div class="stat-card plasmid">
                <div class="number">{}</div>
                <div class="label">Plasmids detected</div>
            </div>
            <div class="stat-card chrom">
                <div class="number">{}</div>
                <div class="label">Chromosomal</div>
            </div>
            <div class="stat-card amr">
                <div class="number">{}</div>
                <div class="pct">{}% of plasmids</div>
                <div class="label">AMR genes</div>
            </div>
            <div class="stat-card bgc">
                <div class="number">{}</div>
                <div class="pct">{}% of plasmids</div>
                <div class="label">BGCs</div>
            </div>
            <div class="stat-card cazy">
                <div class="number">{}</div>
                <div class="pct">{}% of plasmids</div>
                <div class="label">CAZymes</div>
            </div>
            <div class="stat-card tax">
                <div class="number">{}</div>
                <div class="pct">{}% of plasmids</div>
                <div class="label">PLSDB taxonomy</div>
            </div>
        </div>
        <div class="stats-section-label">PlasAnn Functional Categories</div>
        <div class="stats-grid-plasann">
            {}
        </div>
        <div class="section">
            <div class="section-header"><h2>Contextual Network</h2></div>
            <div class="network-container">{}</div>
        </div>
        <div class="section">
            <div class="section-header"><h2>Detected Plasmids</h2></div>
            <table>
                <thead>
                    <tr>
                        <th>Contig ID</th>
                        <th style="text-align:center;">Ensemble Score</th>
                        <th>Resistance (CARD)</th>
                        <th>BGC (MIBiG)</th>
                        <th>CAZyme</th>
                        {}
                        <th>Genus (PLSDB)</th>
                    </tr>
                </thead>
                <tbody>{}</tbody>
            </table>
        </div>
    </div>
    <footer>
        <strong>PlasAnnotatoR</strong> &nbsp;|&nbsp;
        Ensemble: RF(0.9872) + PLASMe(0.9748) + PlasmidHunter(0.9178) + PlasClass(0.9017)
    </footer>
</body>
</html>""".format(
        plasann_card_css,
        plasann_tag_css,
        Path(input_fasta).name if input_fasta else "N/A",
        datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_total, n_plasmids, n_chrom,
        n_amr, pct_amr,
        n_bgc, pct_bgc,
        n_cazy, pct_cazy,
        n_tax, pct_tax,
        plasann_cards_html,
        network_html,
        plasann_th,
        table_rows
    )

    output_file = output_dir / "plasannotator_report.html"
    with open(output_file, 'w') as f:
        f.write(html)
    print("[Report] Report saved to {}".format(output_file))
    return output_file
