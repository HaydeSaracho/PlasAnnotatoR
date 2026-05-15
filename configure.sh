#!/bin/bash

echo "================================================"
echo "  PlasAnnotatoR - Configuration"
echo "================================================"

# Detect micromamba envs directory
MAMBA_ENVS=$(micromamba info | grep "envs directories" | awk '{print $NF}')
if [ -z "$MAMBA_ENVS" ]; then
    echo "ERROR: micromamba not found. Please run install.sh first."
    exit 1
fi

# Detect conda base
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "WARNING: conda not found. PLASMe paths will not be configured."
fi

# Detect PLASMe location
PLASME_DIR=$(find ~ -name "PLASMe.py" 2>/dev/null | head -1 | xargs dirname 2>/dev/null)
if [ -z "$PLASME_DIR" ]; then
    echo "WARNING: PLASMe not found. PLASMe paths will not be configured."
fi

# Generate config.yaml from example
cp config.yaml.example config.yaml

# Replace plasclass paths
sed -i "s|/path/to/envs/plasclass|${MAMBA_ENVS}/plasclass|g" config.yaml

# Replace plasmidhunter paths
sed -i "s|/path/to/envs/plasmidhunter|${MAMBA_ENVS}/plasmidhunter|g" config.yaml

# Replace plasme paths
if [ -n "$CONDA_BASE" ]; then
    sed -i "s|/path/to/envs/plasme|${CONDA_BASE}/envs/plasme|g" config.yaml
fi

# Replace PLASMe script and database paths
if [ -n "$PLASME_DIR" ]; then
    sed -i "s|/path/to/PLASMe/PLASMe.py|${PLASME_DIR}/PLASMe.py|g" config.yaml
    sed -i "s|/path/to/PLASMe/DB|${PLASME_DIR}/DB|g" config.yaml
fi

echo ""
echo "config.yaml generated with the following paths:"
echo ""
cat config.yaml
echo ""
echo "Please verify the paths above before running the pipeline."
echo "If any path is incorrect, edit config.yaml manually."
echo "================================================"
