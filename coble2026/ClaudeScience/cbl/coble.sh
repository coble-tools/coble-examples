#!/usr/bin/env bash

#####################################################
# COBLE:recipe, (c) ICR 2026
# Capture date: 2026-09-19
# Capture time: 20:34:27 BST
# Captured by: rachel.alcraft
#####################################################
# source bashrc for conda
if [ -f ~/.bash_profile ]; then source ~/.bash_profile; elif [ -f ~/.bashrc ]; then source ~/.bashrc; elif command -v conda > /dev/null 2>&1; then eval "$(conda shell.bash hook)"; fi
# Using conda executable conda: /Users/rachel.alcraft/miniforge3/bin/conda
# Using conda alias conda: /Users/rachel.alcraft/miniforge3/bin/conda
#####################################################

conda env remove --name bach-rep -y 2>/dev/null || true
conda create --no-default-packages --name bach-rep -y
export PYTHONNOUSERSITE=1
unset PYTHONPATH
# clean up conda cache first
conda  clean --all -y --force-pkgs-dirs
# deactivate environment
conda deactivate | true
conda deactivate | true
# activate environment
conda activate bach-rep

export PYTHONNOUSERSITE=1
export | grep PYTHONNOUSERSITE
# Channels section
conda config --env --show channels | grep -q 'channels:' && conda config --env --remove-key channels || true
conda config --env --set channel_priority strict
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

# INSTALL SECTION FOR CONDA
#######################################
# COBLE recipe -- Bach et al. 2021, Nat Commun 12:1502
# "Time-resolved single-cell analysis of Brca1 associated mammary tumourigenesis
#  reveals aberrant differentiation of luminal progenitors"
# https://doi.org/10.1038/s41467-021-21783-3
#
# coble build --recipe cbl/.coble.cbl --env bach-rep --rebuild --validate cbl/validate.sh
# coble build --recipe cbl/.coble.cbl --env bach-rep --containers docker,singularity --validate cbl/validate.sh
#######################################
# comments:


# flags:
# Flag: Directive: dependencies, Value: na
conda config --env --set channel_priority strict

# languages:
CONDA_BASE=$(conda info --base)
ARCH=$(uname -m)

conda install -y --solver=libmamba  -c conda-forge 'r-base=3.6.3'
conda install -y --solver=libmamba  r-remotes r-biocmanager r-renv

# compilers:
# Flag: Directive: cran-repo, Value: strict
# The bioconda r36 binaries are linked against libopenblasp-r0.3.7.dylib by
# soname. Left free, the solver takes the current libopenblas and edgeR.dylib
# (hence scran) fails to dyn.load. Pin it before anything Bioconductor lands.

# conda:
conda install -y --solver=libmamba --no-update-deps \
libopenblas=0.3.7 
# bioconda's bioconductor-genomeinfodbdata post-link script fetches its data
# tarball with `curl` and no -L. bioconductor.org has since moved to 302
# redirects for package downloads, so that fetch quietly stores the redirect
# page, the md5 check fails and the whole conda transaction aborts. Give curl
# a config that follows redirects; scoped to this environment and this build.

# bash:
mkdir -p $CONDA_PREFIX/etc/coble
echo location > $CONDA_PREFIX/etc/coble/.curlrc
export CURL_HOME=$CONDA_PREFIX/etc/coble
# Bioconductor 3.10 -- the release that pairs with R 3.6.x, as in the
# authors' Container/ListOfPackages.txt

# bioc-conda:
conda install -y --solver=libmamba --no-update-deps \
'bioconductor-SingleCellExperiment=1.8.0' \
'bioconductor-scater=1.14.0' \
'bioconductor-scran=1.14.1' \
'bioconductor-edgeR=3.28.0' 
# CRAN, at the versions current when the paper was published

# r-conda:
conda install -y --solver=libmamba --no-update-deps \
'r-ggplot2=3.3.3' \
'r-cowplot=1.1.1' \
'r-dplyr=1.0.6' \
'r-plyr=1.8.6' \
'r-ggrepel=0.9.1' \
'r-viridis=0.6.1' \
'r-RColorBrewer=1.1_2' \
'r-ggthemes=4.2.4' \
'r-Matrix=1.3_3' \
'r-svglite=2.0.0' 
# ggrastr's dependencies, taken from conda so that nothing is compiled

# r-conda:
conda install -y --solver=libmamba --no-update-deps \
'r-Cairo=1.5_12.2' \
'r-ggbeeswarm=0.6.0' \
'r-png=0.1_7' \
'r-ragg=0.4.0' 
# ggrastr itself has no conda build for r36; taken from the CRAN archive

# r-package:
Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ggrastr/ggrastr_0.2.3.tar.gz", repos=NULL, type="source")'

# End of recipe
# Validation script setup

# Validate script available in environment at CONDA PREFIX: validate.sh
cp cbl/validate.sh ${CONDA_PREFIX}/bin/validate.sh
chmod +x ${CONDA_PREFIX}/bin/validate.sh
mkdir -p ${CONDA_PREFIX}/coble-recipe
cp cbl/.coble.cbl ${CONDA_PREFIX}/coble-recipe
cp /Users/rachel.alcraft/miniforge3/bin/coble ${CONDA_PREFIX}/bin/
cp /Users/rachel.alcraft/miniforge3/bin/coble-* ${CONDA_PREFIX}/bin/

