#######################################
# COBLE recipe -- Bach et al. 2021, Nat Commun 12:1502
# "Time-resolved single-cell analysis of Brca1 associated mammary tumourigenesis
#  reveals aberrant differentiation of luminal progenitors"
# https://doi.org/10.1038/s41467-021-21783-3
#
# coble build --recipe cbl/.coble.cbl --env bach-rep --rebuild --validate cbl/validate.sh
# coble build --recipe cbl/.coble.cbl --env bach-rep --containers docker,singularity --validate cbl/validate.sh
#######################################
coble:
  - environment: bach-rep

comments:
  Environment to reproduce Figure 1 of Bach et al. 2021 from the authors' own
  code (MarioniLab/Tumorigenesis2018, src/DownstreamAnalysis/Figure1.Rmd) and
  their distributed SingleCellExperiment (BRCA1_SCE.rds).
  .
  INTENT. The authors' analysis container (Container/Tumorigenesis.def) builds
  R 3.6.3 from source on Ubuntu and then takes whatever BiocManager::install()
  offered at the time, which for R 3.6.3 is Bioconductor 3.10. This recipe
  reconstructs that release pair explicitly rather than by date of build:
  R 3.6.3 + Bioconductor 3.10 (scran 1.14.1, scater 1.14.0,
  SingleCellExperiment 1.8.0, edgeR 3.28.0). The CRAN packages are pinned to
  the versions current when the paper was published (March 2021), which is the
  "local version of R and packages" the authors state was used for the
  downstream plotting scripts.
  .
  ggrastr is the one package with no conda build for R 3.6, so it is taken from
  the CRAN archive at 0.2.3, the release current at publication. It has no
  compiled code; its Cairo and ggbeeswarm dependencies come from conda above,
  so nothing is compiled during the build.
  .
  Not included, deliberately: gganimate/gifski. Figure1.Rmd ends with a block
  that renders the decorative repository animation fancyUMAP.gif. It is not a
  panel of Figure 1 and is omitted from the replication script.


channels:
  - bioconda
  - conda-forge

flags:
  - dependencies: NA
  - priority: strict
  # Requires Rosetta installed once beforehand: softwareupdate --install-rosetta
  - export<os=darwin,arch=arm64>: CONDA_SUBDIR=osx-64


languages:
  - r-base=3.6.3@conda-forge

flags:
  - cran-repo: https://packagemanager.posit.co/cran/2021-03-25

# The bioconda r36 binaries are linked against libopenblasp-r0.3.7.dylib by
# soname. Left free, the solver takes the current libopenblas and edgeR.dylib
# (hence scran) fails to dyn.load. Pin it before anything Bioconductor lands.
conda:
  - libopenblas=0.3.7

# bioconda's bioconductor-genomeinfodbdata post-link script fetches its data
# tarball with `curl` and no -L. bioconductor.org has since moved to 302
# redirects for package downloads, so that fetch quietly stores the redirect
# page, the md5 check fails and the whole conda transaction aborts. Give curl
# a config that follows redirects; scoped to this environment and this build.
bash:
  mkdir -p $CONDA_PREFIX/etc/coble
  echo location > $CONDA_PREFIX/etc/coble/.curlrc
  export CURL_HOME=$CONDA_PREFIX/etc/coble

# Bioconductor 3.10 -- the release that pairs with R 3.6.x, as in the
# authors' Container/ListOfPackages.txt
bioc-conda:
  - SingleCellExperiment=1.8.0
  - scater=1.14.0
  - scran=1.14.1
  - edgeR=3.28.0

# CRAN, at the versions current when the paper was published
r-conda:
  - ggplot2=3.3.3
  - cowplot=1.1.1
  - dplyr=1.0.6
  - plyr=1.8.6
  - ggrepel=0.9.1
  - viridis=0.6.1
  - RColorBrewer=1.1_2
  - ggthemes=4.2.4
  - Matrix=1.3_3
  - svglite=2.0.0

# ggrastr's dependencies, taken from conda so that nothing is compiled
r-conda:
  - Cairo=1.5_12.2
  - ggbeeswarm=0.6.0
  - png=0.1_7
  - ragg=0.4.0

# ggrastr itself has no conda build for r36; taken from the CRAN archive
r-package:
  - ggrastr=0.2.3
