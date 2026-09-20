#######################################
# coble build --recipe Figure04/cbl/recipe.cbl --env zheng-gen --rebuild
# coble build --recipe Figure04/cbl/recipe.cbl --env zheng-en --containers docker,singularity --validate validate.sh
#######################################
coble:
  - environment: zheng02
comments:
Zheng2017 with modern toolchain for extended analysis
channels:
  - bioconda
  - conda-forge
compilers:
  - cran-repo: https://packagemanager.posit.co/cran/2025-10-10
flags:
  - ncpus: 1
conda:
  - r-base=4.5.2
r-package:
  - lattice
  - chron
  - DBI
  - BH
  - assertthat
  - MASS
  - digest
  - RColorBrewer
  - dichromat
  - munsell
  - labeling
  - stringi
  - stringr
  - R6
  - lazyeval
  - codetools
r-conda:
  - Rcpp
  - plyr
  - data.table
  - svd
  - SCpubr
  - Seurat
r-package:
  - Matrix
  - gtable
  - scales
  - reshape2
r-package:
  - ggplot2
  - dplyr
  - Rtsne
r-package:
  - pheatmap
conda:
  - cairo
  - pango
  - fonts-anaconda
  - fontconfig
bioc-conda:
  - SingleR
  - celldex
r-conda:
  - remotes
  - pbkrtest
  - lme4
  - car
  - systemfonts
  - textshaping
  - svglite
  - viridis
bioc-conda:
  - ComplexHeatmap
  - glmGamPoi
r-url:
  - https://github.com/jinworks/CellChat/archive/refs/heads/main.tar.gz
r-conda:
  - data.table
r-github:
  - immunogenomics/presto

