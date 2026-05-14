"""
Layer 2: Annotation of plasmids detected by the ensemble
Searches against CARD (AMR), MIBiG (BGCs), CAZy and PLSDB (identity + taxonomy)
"""

import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import yaml


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def extract_plasmid_fasta(ensemble_df, input_fasta, output_fasta):
    """Extracts only contigs classified as plasmids."""
    plasmid_ids = set(ensemble_df[ensemble_df["ensemble_label"] == "plasmid"]["contig_id"])
    count = 0
    with open(output_fasta, 'w') as out:
        for rec in SeqIO.parse(str(input_fasta), "fasta"):
            if rec.id in plasmid_ids:
                out.write(">{}\n{}\n".format(rec.id, str(rec.seq)))
                count += 1
    print("[Annotator] {} plasmids extracted for annotation".format(count))
    return count


def run_blast(query_fasta, db_fasta, output_file, threads=8, task="blastn", evalue=1e-5):
    """Runs BLAST against a database."""
    db_path = Path(db_fasta)
    db_index = Path(str(db_fasta) + ".nin") if task == "blastn" else Path(str(db_fasta) + ".pin")

    if not db_index.exists():
        print("[Annotator] Building BLAST index for {}...".format(db_path.name))
        dbtype = "nucl" if task == "blastn" else "prot"
        subprocess.run([
            "makeblastdb", "-in", str(db_fasta),
            "-dbtype", dbtype, "-out", str(db_fasta)
        ], check=True, capture_output=True)

    cmd = [
        task,
        "-query", str(query_fasta),
        "-db", str(db_fasta),
        "-out", str(output_file),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-max_target_seqs", "5"
    ]
    print("[Annotator] Running {} against {}...".format(task, db_path.name))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[Annotator] Warning {}: {}".format(task, result.stderr[:200]))
    return output_file


def parse_blast_results(blast_file, label):
    """Parses BLAST output format 6."""
    cols = ["contig_id", "subject", "pident", "length", "evalue", "bitscore", "stitle"]
    try:
        df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)
        df = df.sort_values("bitscore", ascending=False).drop_duplicates("contig_id")
        df = df.rename(columns={
            "subject": "{}_subject".format(label),
            "pident": "{}_pident".format(label),
            "evalue": "{}_evalue".format(label),
            "stitle": "{}_annotation".format(label)
        })
        return df[["contig_id", "{}_subject".format(label),
                   "{}_pident".format(label), "{}_evalue".format(label),
                   "{}_annotation".format(label)]]
    except Exception:
        return pd.DataFrame(columns=["contig_id"])


def annotate_with_plsdb(plasmid_fasta, config, output_dir, threads=8):
    """
    BLASTs against PLSDB and enriches with taxonomy.
    Returns DataFrame with columns: contig_id, plsdb_subject, plsdb_pident,
    TAXONOMY_order, TAXONOMY_family, TAXONOMY_genus, TAXONOMY_species
    """
    try:
        plsdb_fasta = Path("data/plsdb/sequences.fasta")
        plsdb_out = output_dir / "plsdb_blast.tsv"

        run_blast(plasmid_fasta, plsdb_fasta, plsdb_out, threads=threads, task="blastn", evalue=1e-10)
        blast_df = parse_blast_results(plsdb_out, "plsdb")

        if blast_df.empty:
            print("[Annotator] PLSDB: no hits found")
            return pd.DataFrame(columns=["contig_id"])

        # Load taxonomy
        tax = pd.read_csv("data/plsdb/meta/taxonomy.csv")
        nuccore = pd.read_csv("data/plsdb/meta/nuccore.csv")

        # Map NUCCORE_ACC to taxonomy via TAXONOMY_UID
        if "NUCCORE_ACC" in nuccore.columns and "TAXONOMY_UID" in nuccore.columns:
            blast_df = blast_df.merge(
                nuccore[["NUCCORE_ACC", "TAXONOMY_UID"]].rename(
                    columns={"NUCCORE_ACC": "plsdb_subject"}),
                on="plsdb_subject", how="left"
            )
            blast_df = blast_df.merge(
                tax[["TAXONOMY_UID", "TAXONOMY_order", "TAXONOMY_family",
                     "TAXONOMY_genus", "TAXONOMY_species"]],
                on="TAXONOMY_UID", how="left"
            )

        print("[Annotator] PLSDB: {} hits with taxonomy".format(blast_df["plsdb_subject"].notna().sum()))
        return blast_df

    except Exception as e:
        print("[Annotator] Warning PLSDB: {}".format(e))
        return pd.DataFrame(columns=["contig_id"])


def predict_proteins(input_fasta, output_dir):
    """Predicts proteins with Prodigal."""
    proteins_fasta = output_dir / "plasmids_proteins.faa"
    cmd = [
        "prodigal",
        "-i", str(input_fasta),
        "-a", str(proteins_fasta),
        "-p", "meta",
        "-q"
    ]
    print("[Annotator] Predicting proteins with Prodigal...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[Annotator] Warning Prodigal: {}".format(result.stderr[:200]))
        return None
    count = sum(1 for line in open(proteins_fasta) if line.startswith(">"))
    print("[Annotator] {} proteins predicted".format(count))
    return proteins_fasta


def run_diamond(query_fasta, db_fasta, output_file, threads=8, evalue=1e-5):
    """Runs DIAMOND blastp against CAZy."""
    db_path = Path(str(db_fasta) + ".dmnd")

    if not db_path.exists():
        print("[Annotator] Building DIAMOND index for CAZy...")
        subprocess.run([
            "diamond", "makedb",
            "--in", str(db_fasta),
            "--db", str(db_fasta),
            "--threads", str(threads),
            "--quiet"
        ], check=True, capture_output=True)

    cmd = [
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(db_fasta),
        "--out", str(output_file),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--max-target-seqs", "5",
        "--quiet"
    ]
    print("[Annotator] Running DIAMOND blastp against CAZy...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[Annotator] Warning DIAMOND: {}".format(result.stderr[:200]))
    return output_file


def run_annotation(ensemble_df, input_fasta, output_dir, config_path="config.yaml", threads=8):
    """
    Runs the complete Layer 2 annotation.

    Args:
        ensemble_df: DataFrame with ensemble results
        input_fasta: original input FASTA
        output_dir: output directory
        config_path: path to config.yaml
        threads: number of threads

    Returns:
        DataFrame with complete annotations
    """
    config = load_config(config_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract plasmids
    plasmid_fasta = output_dir / "plasmids.fasta"
    n_plasmids = extract_plasmid_fasta(ensemble_df, input_fasta, plasmid_fasta)

    if n_plasmids == 0:
        print("[Annotator] No plasmids detected.")
        return ensemble_df

    base = Path(config_path).parent

    # BLAST against CARD
    card_db = base / config['databases']['card']
    card_out = output_dir / "card_blast.tsv"
    run_blast(plasmid_fasta, card_db, card_out, threads=threads, task="blastn")
    card_df = parse_blast_results(card_out, "card")

    # BLAST against MIBiG
    mibig_db = base / config['databases']['mibig']
    mibig_out = output_dir / "mibig_blast.tsv"
    run_blast(plasmid_fasta, mibig_db, mibig_out, threads=threads, task="blastn")
    mibig_df = parse_blast_results(mibig_out, "mibig")

    # DIAMOND blastp against CAZy
    cazy_db = base / config['databases']['cazy']
    proteins_fasta = predict_proteins(plasmid_fasta, output_dir)
    cazy_df = pd.DataFrame(columns=["contig_id"])
    if proteins_fasta and proteins_fasta.exists():
        cazy_out = output_dir / "cazy_diamond.tsv"
        run_diamond(proteins_fasta, cazy_db, cazy_out, threads=threads)
        cazy_df = parse_blast_results(cazy_out, "cazy")
        if not cazy_df.empty:
            cazy_df["contig_id"] = cazy_df["contig_id"].str.rsplit("_", n=1).str[0]
            cazy_df = cazy_df.drop_duplicates("contig_id")

    # BLAST against PLSDB + taxonomy
    plsdb_df = annotate_with_plsdb(plasmid_fasta, config, output_dir, threads=threads)

    # Merge results
    result = ensemble_df.copy()
    if not card_df.empty:
        result = result.merge(card_df, on="contig_id", how="left")
    if not mibig_df.empty:
        result = result.merge(mibig_df, on="contig_id", how="left")
    if not cazy_df.empty and "cazy_subject" in cazy_df.columns:
        result = result.merge(cazy_df, on="contig_id", how="left")
    if not plsdb_df.empty and "plsdb_subject" in plsdb_df.columns:
        result = result.merge(
            plsdb_df[["contig_id", "plsdb_subject", "plsdb_pident",
                      "TAXONOMY_order", "TAXONOMY_family",
                      "TAXONOMY_genus", "TAXONOMY_species"]],
            on="contig_id", how="left"
        )

    output_file = output_dir / "annotation_results.tsv"
    result.to_csv(output_file, sep="\t", index=False)
    print("[Annotator] Results saved to {}".format(output_file))

    return result


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))

    test_fasta = "/home/bionfo/PLASMe/test.fasta"
    real_ids = [rec.id for rec in SeqIO.parse(test_fasta, "fasta")][:2]
    print("Test IDs: {}".format(real_ids))

    test_df = pd.DataFrame({
        "contig_id": real_ids,
        "ensemble_score": [0.95, 0.87],
        "ensemble_label": ["plasmid", "plasmid"]
    })

    result = run_annotation(
        ensemble_df=test_df,
        input_fasta=test_fasta,
        output_dir="/tmp/annotation_test",
        config_path="config.yaml"
    )
    print(result.head())