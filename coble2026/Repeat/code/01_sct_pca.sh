#!/bin/bash
#SBATCH --job-name=01_sct_pca
#SBATCH --output=_log01_sct_pca.out
#SBATCH --error=_log01_sct_pca.err
#SBATCH --ntasks=1
#SBATCH --time=2:30:00
#SBATCH --mem-per-cpu=16084
#SBATCH --cpus-per-task=16
#SBATCH --partition=smp

source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate tnbc-new

in_path="/data/scratch/DBC/UBCN/BCRBIOIN/SHARED/analysis/RachaelNatrajan/TNBC_PDX_residual/snrseq/data/merge_G1_allgene_sct2_scthvg3sample_ward_spearman/seu_CIDER.RDS"

echo "Current working directory: $(pwd)"
echo "Input path: $in_path"

mkdir -p results

# only run if the object does not already exist
if [[ ! -f "results/seu_CIDER_icr.RDS" ]]; then
    echo "## Run sct and pca ################"
    Rscript code/01_sct_pca.R --tag icr --outpath ./results --seuratin $in_path
else
    echo "## Object already exists, skipping ################"
fi

