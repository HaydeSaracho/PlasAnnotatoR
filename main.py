"""
PlasAnnotatoR v1.0 - Main pipeline
Plasmid classification and annotation in metagenomes and bacterial genomes

Usage:
    python main.py -i input.fasta -o results/ -t 8
"""

import argparse
import sys
import yaml
import pandas as pd
from pathlib import Path
from datetime import datetime


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def print_banner():
    print("""
+-------------------------------------------------------+
|           PlasAnnotatoR v1.0                          |
|   Ensemble-based plasmid classification & annotation  |
|   RF(0.881) + PLASMe(0.890) + PlasmidHunter(0.790)    |
|   + PlasClass(0.842)                                  |
+-------------------------------------------------------+
    """)


def parse_args(default_threshold=0.6):
    parser = argparse.ArgumentParser(
        description="PlasAnnotatoR v1.0 - Plasmid classification and annotation"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file (contigs)")
    parser.add_argument("-o", "--output", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("-c", "--config", default="config.yaml",
                        help="Config file (default: config.yaml)")
    parser.add_argument("--threshold", type=float, default=default_threshold,
                        help="Ensemble score threshold for plasmid classification (default: {})".format(default_threshold))
    parser.add_argument("--skip-annotation", action="store_true",
                        help="Skip annotation layer (faster)")
    parser.add_argument("--skip-network", action="store_true",
                        help="Skip network visualization")
    return parser.parse_args()


def run_pipeline(input_fasta, output_dir, threads=8, config_path="config.yaml",
                 skip_annotation=False, skip_network=False, threshold=0.6):
    """
    Runs the complete PlasAnnotatoR pipeline.

    Args:
        input_fasta: path to input FASTA file
        output_dir: output directory
        threads: number of threads
        config_path: path to config.yaml
        skip_annotation: skip Layer 2 annotation
        skip_network: skip Layer 3 network
        threshold: ensemble score threshold (default: read from config.yaml)

    Returns:
        path to HTML report
    """
    start_time = datetime.now()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(Path(__file__).parent))

    print("\n[1/4] Layer 1: Ensemble classifier...")
    from plasannotator.ensemble.voter import run_ensemble
    ensemble_df = run_ensemble(
        input_fasta=input_fasta,
        output_dir=str(output_dir / "ensemble"),
        config_path=config_path,
        threads=threads,
        threshold=threshold
    )

    n_plasmids = (ensemble_df["ensemble_label"] == "plasmid").sum()
    n_total = len(ensemble_df)
    print("[1/4] Done: {}/{} contigs classified as plasmids".format(n_plasmids, n_total))

    ensemble_df.to_csv(output_dir / "ensemble_results.tsv", sep="\t", index=False)

    annotation_df = ensemble_df.copy()

    if not skip_annotation:
        print("\n[2/4] Layer 2: Annotation (CARD + MIBiG + CAZy + PLSDB)...")
        from plasannotator.annotation.annotator import run_annotation
        annotation_df = run_annotation(
            ensemble_df=ensemble_df,
            input_fasta=input_fasta,
            output_dir=str(output_dir / "annotation"),
            config_path=config_path,
            threads=threads
        )
        print("[2/4] Done.")
    else:
        print("\n[2/4] Annotation skipped (--skip-annotation).")

    network_file = None
    if not skip_network and n_plasmids > 0:
        print("\n[3/4] Layer 3: Contextual network...")
        from plasannotator.network.context_net import build_network
        network_file = build_network(
            annotation_df=annotation_df,
            output_dir=str(output_dir / "network"),
            config_path=config_path
        )
        print("[3/4] Done.")
    else:
        print("\n[3/4] Network skipped.")

    print("\n[4/4] Layer 4: Generating HTML report...")
    from plasannotator.report.html_report import generate_report
    report_file = generate_report(
        annotation_df=annotation_df,
        network_html_path=network_file,
        output_dir=str(output_dir),
        input_fasta=input_fasta
    )
    print("[4/4] Done.")

    elapsed = datetime.now() - start_time
    print("\n✅ Pipeline completed in {}".format(str(elapsed).split(".")[0]))
    print("📄 Report: {}".format(report_file))
    print("📁 Results: {}".format(output_dir))

    return report_file


if __name__ == "__main__":
    print_banner()

    # Load config first to get default threshold
    config_path = "config.yaml"
    for i, arg in enumerate(sys.argv):
        if arg in ("-c", "--config") and i + 1 < len(sys.argv):
            config_path = sys.argv[i + 1]

    if not Path(config_path).exists():
        print("ERROR: Config file not found: {}".format(config_path))
        sys.exit(1)

    config = load_config(config_path)
    default_threshold = config.get("ensemble", {}).get("final_threshold", 0.6)

    args = parse_args(default_threshold=default_threshold)

    if not Path(args.input).exists():
        print("ERROR: Input file not found: {}".format(args.input))
        sys.exit(1)

    run_pipeline(
        input_fasta=args.input,
        output_dir=args.output,
        threads=args.threads,
        config_path=args.config,
        skip_annotation=args.skip_annotation,
        skip_network=args.skip_network,
        threshold=args.threshold
    )
