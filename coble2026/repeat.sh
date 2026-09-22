#!/bin/bash
###########################################
# RAlcraft 21/05/2026 - coble repeatability
# RAlcraft 11/06/2026 - update for public repo
###########################################
# This sets up the dependency between the
# jobs as they need different resources
###########################################
cd Repeat

echo "Submitting job 1 (01_sct_pca.R)..."
# Submit job 1, capture its job ID
JOB1=$(sbatch --parsable code/01_sct_pca.sh)

echo "Submitting job 2 (02_umap_plot.R) with dependency on job1=$JOB1..."
# Submit job 2, only runs if job 1 succeeded
sbatch --dependency=afterok:$JOB1 code/02_umap_plot.sh

echo "All jobs submitted to slurm, check with 'squeue --me'"
squeue --me