#######################################
# coble build --recipe Figure03/cbl/recipe.cbl --env zheng-rep --rebuild
# coble build --recipe Figure03/cbl/recipe.cbl --env zheng-rep --containers docker,singularity --validate validate.sh
#######################################
coble:
  - environment: zheng2017
comments:
Zheng2017 with modern toolchain

channels:
  - bioconda
  - conda-forge

flags:
  - cran-repo: https://packagemanager.posit.co/cran/2025-10-10
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