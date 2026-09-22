#!/usr/bin/env bash

#####################################################
# COBLE:recipe, (c) ICR 2026
# Capture date: 2026-09-21
# Capture time: 16:31:42 BST
# Captured by: rachel.alcraft
# Captured on: 35605BRKLT
#####################################################
# source bashrc for conda
if [ -f ~/.bash_profile ]; then source ~/.bash_profile; elif [ -f ~/.bashrc ]; then source ~/.bashrc; elif command -v conda > /dev/null 2>&1; then eval "$(conda shell.bash hook)"; fi
# Using conda executable conda: /Users/rachel.alcraft/miniforge3/bin/conda
# Using conda alias conda: /Users/rachel.alcraft/miniforge3/bin/conda
#####################################################

conda env remove --name scherer-rep -y 2>/dev/null || true
conda create --no-default-packages --name scherer-rep -y
export PYTHONNOUSERSITE=1
unset PYTHONPATH
# clean up conda cache first
conda  clean --all -y --force-pkgs-dirs
# deactivate environment
conda deactivate | true
conda deactivate | true
# activate environment
conda activate scherer-rep

export PYTHONNOUSERSITE=1
export | grep PYTHONNOUSERSITE
# Channels section
conda config --env --remove-key channels 2>/dev/null || true
conda config --env --set channel_priority strict
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

# INSTALL SECTION FOR CONDA
#^######################################
#^ coble build --recipe Figure02/cbl/recipe.cbl --env scherer-rep --rebuild
#^ coble build --recipe Figure02/cbl/recipe.cbl --env scherer-rep --containers docker,singularity --validate validate.sh
#^######################################
# coble:
# coble: environment: scherer-rep
# comments:


# flags:
# Flag: Directive: cran-repo, Value: https://packagemanager.posit.co/cran/2025-10-10
Rscript -e 'options(repos=c(CRAN="https://packagemanager.posit.co/cran/2025-10-10"))'

# conda:
conda install -y --solver=libmamba --no-update-deps \
r-base=4.2.2 
# r-conda:
conda install -y --solver=libmamba --no-update-deps \
'r-ggplot2=3.4.1' \
'r-Seurat=4.3.0' \
'r-ROCR=1.0' \
'r-fossil=0.4.0' \
'r-reshape2=1.4.4' \
'r-RCurl' \
'r-xml2' \
'r-XML' \
'r-locfit' \
'r-corrplot' \
'r-harmony' \
'r-infotheo' \
'r-caret' \
'r-randomForest' \
'r-uwot' 
# bioc-conda:
conda install -y --solver=libmamba --no-update-deps \
'bioconductor-biomaRt' \
'bioconductor-GenomicFeatures' \
'bioconductor-rtracklayer' \
'bioconductor-HDF5Array' \
'bioconductor-Rhdf5lib' \
'bioconductor-rhdf5filters' \
'bioconductor-rhdf5' \
'bioconductor-ComplexHeatmap=2.14.0' \
'bioconductor-GenomeInfoDb' \
'bioconductor-XVector' 
# r-package:
Rscript -e 'install.packages("BiocManager", repos="https://packagemanager.posit.co/cran/2025-10-10", dependencies=NA, Ncpus=1)'
# bioc-package:
Rscript -e 'BiocManager::install("methylumi", dependencies=NA, Ncpus=1)'
Rscript -e 'BiocManager::install("RnBeads", dependencies=NA, Ncpus=1)'
Rscript -e 'BiocManager::install("GenomicRanges", dependencies=NA, Ncpus=1)'
# r-conda:
conda install -y --solver=libmamba --no-update-deps \
'r-pheatmap=1.0.12' \
'r-viridis=0.6.2' 



# End of recipe
# Validation script setup
echo "#!/usr/bin/env bash" > ${CONDA_PREFIX}/bin/validate.sh
echo 'echo "COBLE validation: No script has been specified for scherer-rep environment."' >> ${CONDA_PREFIX}/bin/validate.sh
chmod +x ${CONDA_PREFIX}/bin/validate.sh
chmod +x ${CONDA_PREFIX}/bin/validate.sh
mkdir -p ${CONDA_PREFIX}/coble-recipe
cp Reproduce/cbl/recipe.cbl ${CONDA_PREFIX}/coble-recipe
cp /Users/rachel.alcraft/miniforge3/bin/coble ${CONDA_PREFIX}/bin/
cp /Users/rachel.alcraft/miniforge3/bin/coble-* ${CONDA_PREFIX}/bin/

