#!/bin/bash
set -e

echo "================================================"
echo "  PlasAnnotatoR - Database Download"
echo "================================================"

mkdir -p data/plsdb/meta data/indexes data/models

# RF model
echo ""
echo "[1/5] Downloading RF model..."
wget -O data/models/rf_model.zip "https://zenodo.org/records/20348780/files/rf_model.zip"
unzip data/models/rf_model.zip -d data/models/
rm data/models/rf_model.zip
echo "Done."

# PLSDB 2025 sequences
echo ""
echo "[2/5] Downloading PLSDB 2025 sequences (~1 GB)..."
wget -O data/plsdb/sequences.fasta "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta"
echo "Done."

# PLSDB 2025 metadata
echo ""
echo "[3/5] Downloading PLSDB 2025 metadata..."
wget -O data/plsdb/meta/meta.tar.gz "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz"
tar -xzf data/plsdb/meta/meta.tar.gz -C data/plsdb/meta/
rm data/plsdb/meta/meta.tar.gz
echo "Done."

# CARD
echo ""
echo "[4/5] Downloading CARD..."
wget -O data/indexes/card.tar.bz2 "https://card.mcmaster.ca/latest/data"
tar -xjf data/indexes/card.tar.bz2 -C data/indexes/
rm data/indexes/card.tar.bz2
echo "Done."

# MIBiG 4.0
echo ""
echo "[5/5] Downloading MIBiG 4.0..."
wget -O data/indexes/mibig_4.0.zip "https://zenodo.org/records/20350302/files/mibig_4.0.zip"
    unzip data/indexes/mibig_4.0.zip -d data/indexes/
    rm data/indexes/mibig_4.0.zip
echo "Done."

echo ""
echo "================================================"
echo "  NOTE: CAZy database must be downloaded manually."
echo "  Go to: https://pro.unl.edu/dbCAN2/browse_download.php"
echo "  Download: CAZyDB.07242025.fa"
echo "  Rename to: cazy.fasta"
echo "  Place in: data/indexes/"
echo "================================================"
echo ""
echo "  All databases downloaded successfully!"
echo "  You can now run the pipeline:"
echo "  micromamba activate plasannotator"
echo "  python main.py -i input.fasta -o results/"
echo "================================================"