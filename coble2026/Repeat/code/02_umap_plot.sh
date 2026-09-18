#!/bin/bash
#SBATCH --job-name=02_umap_plot
#SBATCH --output=_log02_umap_plot.out
#SBATCH --error=_log02_umap_plot.err
#SBATCH --ntasks=1
#SBATCH --time=2:30:00
#SBATCH --mem-per-cpu=16084
#SBATCH --cpus-per-task=2
#SBATCH --partition=smp

source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate tnbc-new

mkdir -p results

echo "## Run umap anmd plot ################"
Rscript code/02_umap_plot.R --tag icr --outpath ./results

