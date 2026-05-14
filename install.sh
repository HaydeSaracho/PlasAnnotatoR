#!/bin/bash
set -e

echo "================================================"
echo "  PlasAnnotatoR - Installation"
echo "================================================"

# Check micromamba
if ! command -v micromamba &> /dev/null; then
    echo "ERROR: micromamba not found. Please install micromamba first."
    echo "https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    exit 1
fi

# Check conda/anaconda for PLASMe
if ! command -v conda &> /dev/null; then
    echo "WARNING: conda not found. PLASMe environment will be skipped."
    echo "Install Anaconda/Miniconda to enable PLASMe support."
    SKIP_PLASME=true
fi

echo ""
echo "[1/4] Installing PlasClass environment..."
micromamba env create -f envs/plasclass_env.yml --yes
echo "Done."

echo ""
echo "[2/4] Installing PLASMe environment..."
if [ "$SKIP_PLASME" = true ]; then
    echo "Skipped (conda not found)."
else
    conda env create -f envs/plasme_env.yml
    echo "Installing PLASMe..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate plasme
    pip install git+https://github.com/HubertTang/PLASMe.git
    conda deactivate
    echo "Done."
fi

echo ""
echo "[3/4] Installing PlasmidHunter environment..."
micromamba env create -f envs/plasmidhunter_env.yml --yes
echo "Done."

echo ""
echo "[4/4] Installing PlasAnnotatoR environment..."
micromamba env create -f envs/plasannotator_env.yml --yes
echo "Done."

echo ""
echo "================================================"
echo "  Installation complete!"
echo ""
echo "  Next steps:"
echo "  1. Edit config.yaml with your environment paths"
echo "  2. Download databases (see README.md)"
echo "  3. Run: micromamba activate plasannotator"
echo "  4. Run: python main.py -i input.fasta -o results/"
echo "================================================"
