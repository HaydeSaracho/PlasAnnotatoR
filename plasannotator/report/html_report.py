"""
Layer 4: Integrated HTML report for PlasAnnotatoR
Generates a complete report with statistics, tables and interactive network
"""

import pandas as pd
import shutil
from pathlib import Path
from datetime import datetime


def generate_report(annotation_df, network_html_path, output_dir, input_fasta=""):
    """
    Generates a complete HTML report.

    Args:
        annotation_df: DataFrame with ensemble + annotation results
        network_html_path: path to pyvis network HTML
        output_dir: output directory
        input_fasta: input file name

    Returns:
        path to HTML report
    """
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

    # Copy network to output directory and embed as iframe
    network_html = ""
    if network_html_path and Path(network_html_path).exists():
        network_filename = Path(network_html_path).name
        dest = output_dir / network_filename
        if Path(network_html_path).resolve() != dest.resolve():
            shutil.copy(str(network_html_path), str(dest))
        network_html = '<iframe src="{}" width="100%" height="580px" frameborder="0" style="border-radius:6px; border: 1px solid #dee2e6;"></iframe>'.format(
            network_filename
        )

    # Generate plasmid table
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

        has_annotation = any([
            pd.notna(row.get("card_subject")),
            pd.notna(row.get("mibig_subject")),
            pd.notna(row.get("cazy_subject"))
        ])

        amr_badge = '<span class="tag tag-amr">AMR</span> {}'.format(amr) if pd.notna(row.get("card_subject")) else amr
        bgc_badge = '<span class="tag tag-bgc">BGC</span> {}'.format(bgc) if pd.notna(row.get("mibig_subject")) else bgc
        cazy_badge = '<span class="tag tag-cazy">CAZy</span> {}'.format(cazy) if pd.notna(row.get("cazy_subject")) else cazy

        row_class = "annotated" if has_annotation else ""
        row_bg = 'style="background:#f8f9fa;"' if i % 2 == 0 else ""

        table_rows += """
        <tr class="{}" {}>
            <td style="font-family:monospace; font-size:0.82em;">{}</td>
            <td style="text-align:center;"><span class="score-badge">{}</span></td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
        </tr>""".format(row_class, row_bg, row["contig_id"], score,
                        amr_badge, bgc_badge, cazy_badge, genus)

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

        /* HEADER */
        header {{ background: #ffffff; padding: 24px 48px;
                  border-bottom: 3px solid #c2185b;
                  display: flex; align-items: center; justify-content: space-between;
                  box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
        .header-left h1 {{ color: #c2185b; font-size: 1.8em; font-weight: 700;
                           letter-spacing: -0.5px; }}
        .header-left p {{ color: #6c757d; margin-top: 4px; font-size: 0.88em; }}
        .header-right {{ text-align: right; color: #6c757d; font-size: 0.82em; }}
        .header-right strong {{ color: #c2185b; }}

        /* CONTAINER */
        .container {{ max-width: 1400px; margin: 0 auto; padding: 32px 48px; }}

        /* STAT CARDS */
        .stats-grid {{ display: grid;
                       grid-template-columns: repeat(7, 1fr);
                       gap: 16px; margin-bottom: 32px; }}
        .stat-card {{ background: #ffffff; border-radius: 8px; padding: 18px 12px;
                      text-align: center; border: 1px solid #e0e6ed;
                      box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
        .stat-card .number {{ font-size: 2em; font-weight: 700; line-height: 1; }}
        .stat-card .pct {{ font-size: 0.78em; color: #6c757d; margin-top: 4px; }}
        .stat-card .label {{ color: #6c757d; margin-top: 6px; font-size: 0.78em;
                             font-weight: 500; text-transform: uppercase;
                             letter-spacing: 0.3px; }}
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

        /* SECTIONS */
        .section {{ background: #ffffff; border-radius: 8px; padding: 24px;
                    margin-bottom: 24px; border: 1px solid #e0e6ed;
                    box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
        .section-header {{ display: flex; align-items: center; margin-bottom: 16px;
                           padding-bottom: 12px; border-bottom: 1px solid #e0e6ed; }}
        .section-header h2 {{ color: #c2185b; font-size: 1.05em; font-weight: 600;
                              text-transform: uppercase; letter-spacing: 0.5px; }}

        /* TABLE */
        table {{ width: 100%; border-collapse: collapse; font-size: 0.84em; }}
        th {{ background: #c2185b; color: #ffffff; padding: 10px 12px;
              text-align: left; font-weight: 600; font-size: 0.82em;
              text-transform: uppercase; letter-spacing: 0.4px; }}
        td {{ padding: 8px 12px; border-bottom: 1px solid #e0e6ed;
              vertical-align: middle; }}
        tr:hover td {{ background: #fce4ec !important; }}
        tr.annotated td:first-child {{ border-left: 3px solid #e63946; }}

        /* BADGES & TAGS */
        .score-badge {{ background: #fce4ec; color: #c2185b; padding: 2px 8px;
                        border-radius: 12px; font-size: 0.85em; font-weight: 600;
                        border: 1px solid #f48fb1; }}
        .tag {{ display: inline-block; padding: 1px 6px; border-radius: 4px;
                font-size: 0.72em; font-weight: 700; margin-right: 4px;
                vertical-align: middle; }}
        .tag-amr {{ background: #fde8e4; color: #e76f51; border: 1px solid #f4b8a8; }}
        .tag-bgc {{ background: #d8f3f0; color: #2a9d8f; border: 1px solid #a8ddd8; }}
        .tag-cazy {{ background: #dceaf5; color: #457b9d; border: 1px solid #a8c8e8; }}

        /* NETWORK */
        .network-container {{ width: 100%; height: 580px; border-radius: 6px;
                              overflow: hidden; }}

        /* FOOTER */
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

        <div class="section">
            <div class="section-header">
                <h2>Contextual Network</h2>
            </div>
            <div class="network-container">
                {}
            </div>
        </div>

        <div class="section">
            <div class="section-header">
                <h2>Detected Plasmids</h2>
            </div>
            <table>
                <thead>
                    <tr>
                        <th>Contig ID</th>
                        <th style="text-align:center;">Ensemble Score</th>
                        <th>Resistance (CARD)</th>
                        <th>BGC (MIBiG)</th>
                        <th>CAZyme</th>
                        <th>Genus (PLSDB)</th>
                    </tr>
                </thead>
                <tbody>
                    {}
                </tbody>
            </table>
        </div>

    </div>

    <footer>
        <strong>PlasAnnotatoR</strong> &nbsp;|&nbsp;
        Ensemble: RF(0.9872) + PLASMe(0.9748) + PlasmidHunter(0.9178) + PlasClass(0.9017)
    </footer>
</body>
</html>""".format(
        Path(input_fasta).name if input_fasta else "N/A",
        datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_total, n_plasmids, n_chrom,
        n_amr, pct_amr,
        n_bgc, pct_bgc,
        n_cazy, pct_cazy,
        n_tax, pct_tax,
        network_html,
        table_rows
    )

    output_file = output_dir / "plasannotator_report.html"
    with open(output_file, 'w') as f:
        f.write(html)
    print("[Report] Report saved to {}".format(output_file))
    return output_file


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from plasannotator.network.context_net import build_network

    test_df = pd.DataFrame({
        "contig_id": ["NODE_5_plasmid", "NODE_42_plasmid", "NODE_99_plasmid",
                      "NODE_1_chrom", "NODE_2_chrom"],
        "ensemble_score": [0.98, 0.91, 0.85, 0.12, 0.08],
        "ensemble_label": ["plasmid", "plasmid", "plasmid", "chromosome", "chromosome"],
        "card_subject": ["ARO:3000805", None, "ARO:3002645", None, None],
        "card_pident": [95.2, None, 88.4, None, None],
        "card_annotation": ["tetA tetracycline efflux pump", None, "blaTEM-1 beta-lactamase", None, None],
        "mibig_subject": [None, "BGC0000535", None, None, None],
        "mibig_pident": [None, 76.3, None, None, None],
        "mibig_annotation": [None, "erythromycin biosynthetic gene cluster", None, None, None],
        "cazy_subject": ["GH163|CAZy", None, "GT4|CAZy", None, None],
        "cazy_pident": [53.7, None, 61.2, None, None],
        "cazy_annotation": ["AEU37665.1|GH163", None, "GT4 glycosyltransferase", None, None],
        "TAXONOMY_genus": ["Escherichia", None, "Klebsiella", None, None],
        "plsdb_subject": ["NZ_CP012345.1", None, "NZ_CP098765.1", None, None],
        "plsdb_pident": [98.5, None, 91.2, None, None]
    })

    network_file = build_network(test_df, "/tmp/report_test")

    output_file = generate_report(
        annotation_df=test_df,
        network_html_path=network_file,
        output_dir="/tmp/report_test",
        input_fasta="test_metagenome.fasta"
    )

    import subprocess
    subprocess.Popen(["xdg-open", str(output_file)])
    print("Report: {}".format(output_file))