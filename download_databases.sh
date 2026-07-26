#!/bin/bash
set -e
echo "================================================"
echo "  PlasAnnotatoR - Database Download"
echo "================================================"
mkdir -p data/plsdb/meta data/indexes data/models data/plasann data/plasme

# RF model
echo ""
echo "[1/7] Downloading RF model..."
wget -O data/models/rf_model.pkl "https://zenodo.org/records/21148201/files/rf_model.pkl?download=1"
echo "Done."

# PLSDB 2025 sequences
echo ""
echo "[2/7] Downloading PLSDB 2025 sequences (~7 GB)..."
wget -O data/plsdb/sequences.fasta "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta"
echo "Done."

# PLSDB 2025 metadata
echo ""
echo "[3/7] Downloading PLSDB 2025 metadata (~3.4 GB)..."
wget -O data/plsdb/meta/meta.tar.gz "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz"
tar -xzf data/plsdb/meta/meta.tar.gz -C data/plsdb/meta/
rm data/plsdb/meta/meta.tar.gz
echo "Done."

# CARD
echo ""
echo "[4/7] Downloading CARD..."
wget -O data/indexes/card.tar.bz2 "https://card.mcmaster.ca/latest/data"
tar -xjf data/indexes/card.tar.bz2 -C data/indexes/
rm data/indexes/card.tar.bz2
echo "Done."

# MIBiG 4.0
echo ""
echo "[5/7] Downloading MIBiG 4.0 nucleotide sequences..."
wget -O data/indexes/mibig_4.0.zip "https://zenodo.org/records/20350302/files/mibig_4.0.zip"
unzip data/indexes/mibig_4.0.zip -d data/indexes/
rm data/indexes/mibig_4.0.zip
echo "Done."

# CAZy
echo ""
echo "[6/7] Downloading CAZy database (CAZyDB.07242025, 1.2 GB)..."
echo "Note: CAZy is updated regularly. This version corresponds to July 2025."
echo "To use a more recent version, download manually from:"
echo "https://pro.unl.edu/dbCAN2/browse_download.php"
echo "and rename the file to cazy.fasta in data/indexes/"
wget -O data/indexes/cazy.zip "https://zenodo.org/records/20350742/files/cazy.zip"
unzip data/indexes/cazy.zip -d data/indexes/
rm data/indexes/cazy.zip
echo "Done."

# PlasAnn functional databases
echo ""
echo "[7/7] Downloading PlasAnn functional databases..."
wget -O data/plasann/plasann_databases.zip "https://zenodo.org/records/20501577/files/plasann_databases.zip"
unzip data/plasann/plasann_databases.zip -d data/plasann/
rm data/plasann/plasann_databases.zip
echo "Done."

# PLASMe database
echo ""
echo "Downloading PLASMe database (~12.4 GB)..."
wget -O data/plasme/DB.zip "https://zenodo.org/records/8046934/files/DB.zip"
unzip ~/PLASMe/DB.zip -d ~/PLASMe/
rm ~/PLASMe/DB.zip
echo "Done."

echo ""
echo "================================================"
echo "  All databases downloaded successfully!"
echo "  You can now run the pipeline:"
echo "  micromamba activate plasannotator"
echo "  python main.py -i input.fasta -o results/"
echo "================================================"
