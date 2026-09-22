#!/bin/bash
# Run this from the folder coble2026
# >$ bash b_reproduce.sh singularity

ENV_NAME="scherer-rep"
SING_NAME="../cbl-scherer-rep.sif"
ENV_MODE=$1
source ~/.bashrc
cd Reproduce

# Download data needed for the scripts
echo "### 1a. Downloading data from figshare and saving to data/ directory"
mkdir -p data
if [[ ! -f "data/larry_seurat.rds" ]]; then
    curl -L -o data/larry_seurat.rds 'https://api.figshare.com/v2/file/download/42479346'
else
    echo "data/larry_seurat.rds already exists, skipping"
fi
echo "### 1b. Downloading data from github and saving to data/ directory"
if [[ ! -f "data/panel_info_dropout_pwm.tsv" ]]; then
    curl -L -o data/panel_info_dropout_pwm.tsv 'https://raw.githubusercontent.com/veltenlab/EPI-clone/refs/heads/main/infos/panel_info_dropout_pwm.tsv'
else
    echo "data/panel_info_dropout_pwm.tsv already exists, skipping"
fi

# Running in either singularity or conda mode
if [ "$ENV_MODE" = "singularity" ] || [ "$ENV_MODE" = "apptainer" ]; then
    echo "### 2. Running in $ENV_MODE env ###"
    $ENV_MODE exec $SING_NAME Rscript code/fig1.R
else
    echo "### 2. Running in conda env ###"
    conda activate $ENV_NAME
    Rscript code/fig1.R
fi

# Completed!
echo "### 3. The figures should now be saved in the folder Fig2_PanelB/results/ ###"