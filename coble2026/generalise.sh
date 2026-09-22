#!/bin/bash
# Run this from the folder coble2026
# >$ bash d_generalise.sh singularity

ENV_MODE=$1
source ~/.bashrc
cd Generalise

echo "### 1. Downloading data from 10xgenomics and saving to data/ directory"
mkdir -p data
if [[ ! -f "data/pbmc68k_data.rds" ]]; then
    wget -O data/pbmc68k_data.rds https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/pbmc68k_data.rds
fi
if [[ ! -f "data/all_pure_select_11types.rds" ]]; then
    wget -O data/all_pure_select_11types.rds https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/all_pure_select_11types.rds
fi

sbatch code/extend.sh $ENV_MODE

# Completed!
echo "### 3. The batch has been submitted, the figures will be saved in the folder Fig2_PanelD/results/ ###"