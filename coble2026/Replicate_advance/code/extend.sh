#!/bin/bash
#SBATCH --job-name=zheng-gen
#SBATCH --output=_log-slurm.out
#SBATCH --error=_log-slurm.err
#SBATCH --ntasks=1
#SBATCH --time=2:30:00
#SBATCH --mem-per-cpu=8042
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute

ENV_NAME="zheng-gen"
SING_NAME="../cbl-zheng-gen.sif"
ENV_MODE=$1

# Running in either singularity or conda mode
if [ "$ENV_MODE" = "singularity" ] || [ "$ENV_MODE" = "apptainer" ]; then
    echo "### 2. Running in $ENV_MODE env ###"
    $ENV_MODE exec $SING_NAME Rscript code/extend.R --celltype=all
    $ENV_MODE exec $SING_NAME Rscript code/extend.R --celltype=tcells
else
    echo "### 2. Running in conda env ###"
    source ~/.bashrc
    conda activate $ENV_NAME
    Rscript code/extend.R --celltype=all
    Rscript code/extend.R --celltype=tcells
fi