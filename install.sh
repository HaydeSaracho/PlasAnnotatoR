#!/bin/bash
set -e

echo "================================================"
echo "  PlasAnnotatoR - Installation"
echo "================================================"

# Check micromamba
if ! command -v micromamba &> /dev/null; then
    echo "ERROR: micromamba not found. Please install micromamba first."
    echo "See README.md for installation instructions."
    exit 1
fi

# Check conda for PLASMe
if ! command -v conda &> /dev/null; then
    echo "WARNING: conda not found. PLASMe environment will be skipped."
    echo "Install Miniconda to enable PLASMe support."
    SKIP_PLASME=true
fi

echo ""
echo "[1/4] Installing PlasClass environment..."
micromamba env create -f envs/plasclass_env.yml --yes
echo "Installing PlasClass..."
~/micromamba/envs/plasclass/bin/pip install --no-deps git+https://github.com/Shamir-Lab/PlasClass.git
echo "Done."

echo ""
echo "[2/4] Installing PLASMe environment..."
if [ "$SKIP_PLASME" = true ]; then
    echo "Skipped (conda not found)."
else
    conda env create -f envs/plasme_env.yml
    echo "Cloning PLASMe..."
    git clone https://github.com/HubertTang/PLASMe.git ~/PLASMe
    echo "Downloading PLASMe database..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate plasme
    cd ~/PLASMe
    python PLASMe.py download
    cd -
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
echo "  1. cp config.yaml.example config.yaml"
echo "  2. Edit config.yaml with your environment paths"
echo "  3. Download databases (see README.md)"
echo "  4. micromamba activate plasannotator"
echo "  5. python main.py -i input.fasta -o results/"
echo "================================================"
