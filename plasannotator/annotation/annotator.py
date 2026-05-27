"""
Layer 2: Annotation of plasmids detected by the ensemble
Searches against CARD (AMR), MIBiG (BGCs), CAZy, PLSDB (taxonomy)
and PlasAnn functional databases (conjugation, virulence, backbone,
stress response, DNA mobility, metal/biocide resistance, toxin-antitoxin)
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


def parse_blast_results(blast_file, label, min_pident=70.0, min_coverage=70.0):
    """Parses BLAST output format 6."""
    cols = ["contig_id", "subject", "pident", "length", "evalue", "bitscore", "stitle"]
    try:
        df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)
        df = df[df["pident"] >= min_pident]
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
    """BLASTs against PLSDB and enriches with taxonomy."""
    try:
        plsdb_fasta = Path("data/plsdb/sequences.fasta")
        plsdb_out = output_dir / "plsdb_blast.tsv"

        run_blast(plasmid_fasta, plsdb_fasta, plsdb_out, threads=threads, task="blastn", evalue=1e-10)
        blast_df = parse_blast_results(plsdb_out, "plsdb")

        if blast_df.empty:
            print("[Annotator] PLSDB: no hits found")
            return pd.DataFrame(columns=["contig_id"])

        tax = pd.read_csv("data/plsdb/meta/taxonomy.csv")
        nuccore = pd.read_csv("data/plsdb/meta/nuccore.csv")

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

        print("[Annotator] PLSDB: {} hits with taxonomy".format(
            blast_df["plsdb_subject"].notna().sum()))
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


def run_diamond(query_fasta, db_path, output_file, threads=8, evalue=1e-5):
    """Runs DIAMOND blastp."""
    dmnd_path = Path(str(db_path) + ".dmnd") if not str(db_path).endswith(".dmnd") else Path(db_path)

    if not dmnd_path.exists():
        print("[Annotator] Building DIAMOND index for {}...".format(dmnd_path.name))
        subprocess.run([
            "diamond", "makedb",
            "--in", str(db_path),
            "--db", str(db_path),
            "--quiet"
        ], check=True, capture_output=True)

    cmd = [
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(dmnd_path),
        "--out", str(output_file),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--max-target-seqs", "5",
        "--quiet"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[Annotator] Warning DIAMOND: {}".format(result.stderr[:200]))
    return output_file


def parse_diamond_plasann(blast_file, label, min_pident=70.0):
    """
    Parses DIAMOND output against PlasAnn databases.
    Subject ID format: ProteinID|GeneName|Product
    Returns best hit gene name per contig.
    """
    cols = ["protein_id", "subject", "pident", "length", "evalue", "bitscore", "stitle"]
    try:
        df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)
        df = df[df["pident"] >= min_pident]
        if df.empty:
            return pd.DataFrame(columns=["contig_id"])

        # Extract gene name from subject ID (format: ProteinID|GeneName|Product)
        df["gene_name"] = df["subject"].str.split("|").str[1]
        df["product"] = df["subject"].str.split("|").str[2].str.replace("_", " ")

        # Map protein_id back to contig_id (Prodigal adds _1, _2, etc.)
        df["contig_id"] = df["protein_id"].str.rsplit("_", n=1).str[0]

        df = df.sort_values("bitscore", ascending=False).drop_duplicates("contig_id")
        df = df.rename(columns={
            "gene_name": "{}_gene".format(label),
            "product": "{}_product".format(label),
            "pident": "{}_pident".format(label),
            "evalue": "{}_evalue".format(label)
        })
        return df[["contig_id", "{}_gene".format(label),
                   "{}_product".format(label), "{}_pident".format(label),
                   "{}_evalue".format(label)]]
    except Exception as e:
        print("[Annotator] Warning parsing {}: {}".format(label, e))
        return pd.DataFrame(columns=["contig_id"])


def run_plasann_annotation(proteins_fasta, output_dir, threads=8):
    """
    Runs DIAMOND blastp against all PlasAnn functional databases.
    Returns dict of DataFrames keyed by category label.
    """
    plasann_dir = Path("data/plasann")
    categories = {
        "conjugation": "conjugation",
        "virulence_defense": "virulence_defense",
        "plasmid_backbone": "plasmid_backbone",
        "stress_response": "stress_response",
        "dna_mobility": "dna_mobility",
        "metal_biocide": "metal_biocide",
        "toxin_antitoxin": "toxin_antitoxin"
    }

    results = {}
    for label, filename in categories.items():
        db_path = plasann_dir / "plasann_{}.dmnd".format(filename)
        if not db_path.exists():
            print("[Annotator] PlasAnn DB not found: {}".format(db_path))
            continue
        out_file = output_dir / "plasann_{}.tsv".format(label)
        print("[Annotator] Running DIAMOND against PlasAnn {}...".format(label))
        run_diamond(proteins_fasta, db_path, out_file, threads=threads)
        df = parse_diamond_plasann(out_file, label)
        if not df.empty:
            results[label] = df
            print("[Annotator] PlasAnn {}: {} hits".format(label, len(df)))

    return results


def run_annotation(ensemble_df, input_fasta, output_dir, config_path="config.yaml", threads=8):
    """Runs the complete Layer 2 annotation."""
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

    # Predict proteins for DIAMOND searches
    proteins_fasta = predict_proteins(plasmid_fasta, output_dir)

    # DIAMOND against CAZy
    cazy_df = pd.DataFrame(columns=["contig_id"])
    if proteins_fasta and proteins_fasta.exists():
        cazy_db = base / config['databases']['cazy']
        cazy_out = output_dir / "cazy_diamond.tsv"
        run_diamond(proteins_fasta, cazy_db, cazy_out, threads=threads)
        cazy_df = parse_blast_results(cazy_out, "cazy")
        if not cazy_df.empty:
            cazy_df["contig_id"] = cazy_df["contig_id"].str.rsplit("_", n=1).str[0]
            cazy_df = cazy_df.drop_duplicates("contig_id")

    # DIAMOND against PlasAnn functional databases
    plasann_results = {}
    if proteins_fasta and proteins_fasta.exists():
        plasann_results = run_plasann_annotation(proteins_fasta, output_dir, threads=threads)

    # BLAST against PLSDB + taxonomy
    plsdb_df = annotate_with_plsdb(plasmid_fasta, config, output_dir, threads=threads)

    # Merge all results
    result = ensemble_df.copy()
    if not card_df.empty:
        result = result.merge(card_df, on="contig_id", how="left")
    if not mibig_df.empty:
        result = result.merge(mibig_df, on="contig_id", how="left")
    if not cazy_df.empty and "cazy_subject" in cazy_df.columns:
        result = result.merge(cazy_df, on="contig_id", how="left")
    for label, df in plasann_results.items():
        if not df.empty and "contig_id" in df.columns:
            result = result.merge(df, on="contig_id", how="left")
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
