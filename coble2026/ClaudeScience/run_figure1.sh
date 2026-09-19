#!/usr/bin/env bash
# =============================================================================
# Reproduce Figure 1 of Bach et al. 2021 (Nat Commun 12:1502) in the COBLE
# environment built from cbl/.coble.cbl.
#
#   bash run_figure1.sh              # conda (default)
#   bash run_figure1.sh singularity  # cbl-bach-rep.sif
#   bash run_figure1.sh docker       # cbl-bach-rep
#
# Build the environment first:
#   coble build --recipe cbl/.coble.cbl --env bach-rep --rebuild --validate cbl/validate.sh
# =============================================================================
# NB: no `set -u`. Conda's own activation hooks (activate-gfortran_*.sh) read
# unset variables, so `conda activate` aborts under nounset.
set -eo pipefail

ENV_NAME="bach-rep"
SIF_NAME="cbl-${ENV_NAME}.sif"
IMG_NAME="cbl-${ENV_NAME}"
ENV_MODE="${1:-conda}"

cd "$(dirname "$0")"
export SCE_PATH="data/BRCA1_SCE.rds"
export OUT_DIR="results/Figure1"
export CODE_DIR="code"
mkdir -p data results

# --- 1. input data -----------------------------------------------------------
# The authors' processed SingleCellExperiment, as distributed by the Marioni
# lab. Every script in src/DownstreamAnalysis reads this object.
if [[ ! -f "$SCE_PATH" ]]; then
    echo "### 1. Downloading BRCA1_SCE.rds (2.0 GB) ###"
    curl -L --fail -o "$SCE_PATH" \
      'https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/BRCA1_SCE.rds'
else
    echo "### 1. $SCE_PATH present, skipping download ###"
fi

# --- 2. run ------------------------------------------------------------------
echo "### 2. Running Figure1.R in $ENV_MODE ###"
case "$ENV_MODE" in
  singularity|apptainer)
    "$ENV_MODE" exec \
      --env SCE_PATH="$SCE_PATH",OUT_DIR="$OUT_DIR",CODE_DIR="$CODE_DIR" \
      "$SIF_NAME" Rscript code/Figure1.R
    ;;
  docker)
    docker run --rm -v "$PWD":/workspace -w /workspace \
      -e SCE_PATH -e OUT_DIR -e CODE_DIR \
      "$IMG_NAME" Rscript code/Figure1.R
    ;;
  conda)
    # shellcheck disable=SC1091
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    coble validate || echo "(coble validate reported differences - see above)"
    Rscript code/Figure1.R
    ;;
  *)
    echo "Unknown mode '$ENV_MODE' (expected conda, docker or singularity)" >&2
    exit 1
    ;;
esac

echo "### 3. Figures written to $OUT_DIR/ ###"
ls -la "$OUT_DIR"
