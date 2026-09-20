# COBLE:export, (c) ICR 2026
# Capture date: 2026-09-19
# Capture time: 20:39:47 BST
# Captured by: rachel.alcraft

coble:

  - environment: bach-rep

channels:
  - bioconda
  - conda-forge

languages:
  - r-base=3.6.3@conda-forge
flags:
  - compile-tools: true
  - dependencies: false
  - priority: flexible

conda:
  - clangxx=14.0.6@conda-forge
  - icu=67.1@conda-forge
  - libcblas=3.8.0@conda-forge
  - libclang-cpp14=14.0.6@conda-forge
  - libcurl=7.76.1@conda-forge
  - libcxx=23.1.1@conda-forge
  - libgcc=16.2.0@conda-forge
  - libzlib=1.2.13@conda-forge
  - zlib=1.2.13@conda-forge
  - bwidget=1.10.1@conda-forge
  - bzip2=1.0.8@conda-forge
  - c-ares=1.34.8@conda-forge
  - ca-certificates=2026.7.22@conda-forge
  - cairo=1.16.0@conda-forge
  - clang-14=14.0.6@conda-forge
  - clang=14.0.6@conda-forge
  - compiler-rt=14.0.6@conda-forge
  - curl=7.76.1@conda-forge
  - font-ttf-dejavu-sans-mono=2.37@conda-forge
  - font-ttf-inconsolata=3.000@conda-forge
  - font-ttf-source-code-pro=2.038@conda-forge
  - font-ttf-ubuntu=0.83@conda-forge
  - fontconfig=2.13.1@conda-forge
  - fonts-conda-forge=1@conda-forge
  - freetype=2.12.1@conda-forge
  - fribidi=1.0.16@conda-forge
  - gmp=6.3.0@conda-forge
  - graphite2=1.3.15@conda-forge
  - gsl=2.6@conda-forge
  - harfbuzz=2.7.2@conda-forge
  - isl=0.22.1@conda-forge
  - jpeg=9e@conda-forge
  - krb5=1.17.2@conda-forge
  - lerc=4.2.0@conda-forge
  - libblas=3.8.0@conda-forge
  - libdeflate=1.25@conda-forge
  - libedit=3.1.20250104@conda-forge
  - libev=4.33@conda-forge
  - libffi=3.7.0@conda-forge
  - libgfortran=4.0.0@conda-forge
  - libgfortran4=7.5.0@conda-forge
  - libglib=2.80.2@conda-forge
  - libiconv=1.18@conda-forge
  - libintl=0.25.1@conda-forge
  - liblapack=3.8.0@conda-forge
  - libllvm14=14.0.6@conda-forge
  - liblzma-devel=5.8.3@conda-forge
  - liblzma=5.8.3@conda-forge
  - libnghttp2=1.51.0@conda-forge
  - libopenblas=0.3.7@conda-forge
  - libpng=1.6.43@conda-forge
  - libssh2=1.10.0@conda-forge
  - libtiff=4.4.0@conda-forge
  - libxcb=1.17.0@conda-forge
  - libxml2=2.9.10@conda-forge
  - llvm-openmp=23.1.1@conda-forge
  - llvm-tools=14.0.6@conda-forge
  - make=4.4.1@conda-forge
  - mpc=1.4.0@conda-forge
  - mpfr=4.2.2@conda-forge
  - ncurses=6.6@conda-forge
  - openssl=1.1.1w@conda-forge
  - pango=1.42.4@conda-forge
  - pcre=8.45@conda-forge
  - pcre2=10.43@conda-forge
  - pixman=0.46.4@conda-forge
  - pthread-stubs=0.4@conda-forge
  - readline=8.3@conda-forge
  - sigtool=0.1.3@conda-forge
  - tapi=1100.0.11@conda-forge
  - tk=8.6.13@conda-forge
  - tktable=2.10@conda-forge
  - xorg-libice=1.1.2@conda-forge
  - xorg-libsm=1.2.6@conda-forge
  - xorg-libx11=1.8.13@conda-forge
  - xorg-libxau=1.0.12@conda-forge
  - xorg-libxdmcp=1.1.5@conda-forge
  - xorg-libxt=1.3.1@conda-forge
  - xz-gpl-tools=5.8.3@conda-forge
  - xz-tools=5.8.3@conda-forge
  - xz=5.8.3@conda-forge
  - zstd=1.5.6@conda-forge

r-conda:
  - assertthat=0.2.1@conda-forge
  - backports=1.2.1@conda-forge
  - beeswarm=0.3.1@conda-forge
  - bh=1.75.0_0@conda-forge
  - biocmanager=1.30.15@conda-forge
  - bitops=1.0_7@conda-forge
  - brio=1.1.2@conda-forge
  - callr=3.7.0@conda-forge
  - cli=2.5.0@conda-forge
  - colorspace=2.0_1@conda-forge
  - cowplot=1.1.1@conda-forge
  - cpp11=0.2.7@conda-forge
  - crayon=1.4.1@conda-forge
  - desc=1.3.0@conda-forge
  - diffobj=0.3.4@conda-forge
  - digest=0.6.27@conda-forge
  - dplyr=1.0.6@conda-forge
  - dqrng=0.3.0@conda-forge
  - ellipsis=0.3.2@conda-forge
  - evaluate=0.14@conda-forge
  - fansi=0.4.2@conda-forge
  - farver=2.1.0@conda-forge
  - formatr=1.9@conda-forge
  - futile.logger=1.4.3@conda-forge
  - futile.options=1.0.1@conda-forge
  - gdtools=0.2.2@conda-forge
  - generics=0.1.0@conda-forge
  - ggbeeswarm=0.6.0@conda-forge
  - ggplot2=3.3.3@conda-forge
  - ggrepel=0.9.1@conda-forge
  - ggthemes=4.2.4@conda-forge
  - glue=1.4.2@conda-forge
  - gridextra=2.3@conda-forge
  - gtable=0.3.0@conda-forge
  - igraph=1.2.6@conda-forge
  - irlba=2.3.3@conda-forge
  - isoband=0.2.4@conda-forge
  - jsonlite=1.7.2@conda-forge
  - labeling=0.4.2@conda-forge
  - lambda.r=1.2.4@conda-forge
  - lattice=0.20_44@conda-forge
  - lifecycle=1.0.0@conda-forge
  - locfit=1.5_9.4@conda-forge
  - magrittr=2.0.1@conda-forge
  - mass=7.3_54@conda-forge
  - matrix=1.3_3@conda-forge
  - matrixstats=0.58.0@conda-forge
  - mgcv=1.8_35@conda-forge
  - munsell=0.5.0@conda-forge
  - nlme=3.1_150@conda-forge
  - pillar=1.6.1@conda-forge
  - pkgconfig=2.0.3@conda-forge
  - pkgload=1.2.1@conda-forge
  - plyr=1.8.6@conda-forge
  - png=0.1_7@conda-forge
  - praise=1.0.0@conda-forge
  - processx=3.5.2@conda-forge
  - ps=1.6.0@conda-forge
  - purrr=0.3.4@conda-forge
  - r6=2.5.0@conda-forge
  - ragg=0.4.0@conda-forge
  - rcolorbrewer=1.1_2@conda-forge
  - rcpp=1.0.6@conda-forge
  - rcppannoy=0.0.18@conda-forge
  - rcpphnsw=0.3.0@conda-forge
  - rcurl=1.98_1.3@conda-forge
  - rematch2=2.1.2@conda-forge
  - remotes=2.3.0@conda-forge
  - renv=0.13.2@conda-forge
  - rlang=0.4.11@conda-forge
  - rprojroot=2.0.2@conda-forge
  - rstudioapi=0.13@conda-forge
  - rsvd=1.0.3@conda-forge
  - scales=1.1.1@conda-forge
  - sitmo=2.0.1@conda-forge
  - snow=0.4_3@conda-forge
  - statmod=1.4.35@conda-forge
  - stringi=1.5.3@conda-forge
  - stringr=1.4.0@conda-forge
  - svglite=2.0.0@conda-forge
  - testthat=3.0.2@conda-forge
  - textshaping=0.1.2@conda-forge
  - tibble=3.1.2@conda-forge
  - tidyselect=1.1.1@conda-forge
  - utf8=1.2.1@conda-forge
  - vctrs=0.3.8@conda-forge
  - vipor=0.4.5@conda-forge
  - viridis=0.6.1@conda-forge
  - viridislite=0.4.0@conda-forge
  - waldo=0.2.5@conda-forge
  - withr=2.4.2@conda-forge

bioc-conda:
  - beachmat=2.2.0@bioconda
  - biocgenerics=0.32.0@bioconda
  - biocneighbors=1.4.0@bioconda
  - biocparallel=1.20.0@bioconda
  - biocsingular=1.2.0@bioconda
  - delayedarray=0.12.0@bioconda
  - delayedmatrixstats=1.8.0@bioconda
  - edger=3.28.0@bioconda
  - genomeinfodb=1.22.0@bioconda
  - genomeinfodbdata=1.2.2@bioconda
  - genomicranges=1.38.0@bioconda
  - hdf5array=1.14.0@bioconda
  - iranges=2.20.0@bioconda
  - limma=3.42.0@bioconda
  - rhdf5=2.30.0@bioconda
  - rhdf5lib=1.8.0@bioconda
  - s4vectors=0.24.0@bioconda
  - scater=1.14.0@bioconda
  - scran=1.14.1@bioconda
  - singlecellexperiment=1.8.0@bioconda
  - summarizedexperiment=1.16.0@bioconda
  - xvector=0.26.0@bioconda
  - zlibbioc=1.32.0@bioconda

r-package:

r-package:
  - ggrastr=0.2.3

pip:
  - -e git+https://oauth2:glpat-vTcELeeNMhOKFBdsB4BV0G86MQp1OmFoCA.01.0y1bnsxzf@git.icr.ac.uk/rae-private/bcds/ai-widget.git@678eb9cf34838c59025550986a4cc378026f50b4#egg==ai_streamlit
  - altair==6.2.2
  - anthropic==1.3.0
  - attrs==26.1.0
  - blinker==1.9.0
  - click==8.4.2
  - contourpy==1.3.3
  - cycler==0.12.1
  - docstring_parser==0.18.0
  - fonttools==4.63.0

pip:
  - git+https://github.com/ICR-RSE-Group/gitalma.git@71ee59c0a415837ff306514dd08e1b3bb4548ee2@git+https://github.com/ICR-RSE-Group/gitalma.git@71ee59c0a415837ff306514dd08e1b3bb4548ee2

pip:
  - gitdb==4.0.12
  - GitPython==3.1.56
  - h5py==3.16.0
  - httpcore2==2.12.0
  - httptools==0.8.0
  - httpx-oauth==0.15.1
  - httpx2==2.12.0
  - idna==3.19
  - itsdangerous==2.2.0
  - jiter==0.16.0
  - jsonschema-specifications==2025.9.1
  - jsonschema==4.26.0
  - kiwisolver==1.5.0
  - matplotlib==3.11.1
  - narwhals==2.24.0
  - numpy==2.5.1
  - openslide-bin==4.0.1.2
  - openslide-python==1.4.6
  - pandas==3.0.5
  - pillow==12.3.0
  - protobuf==7.35.1
  - py-spy==0.4.2
  - pyarrow==24.0.0
  - pydeck==0.9.3
  - PyJWT==2.13.0
  - pyparsing==3.3.2
  - python-dotenv==1.0.1
  - python-multipart==0.0.32
  - PyYAML==6.0.3
  - referencing==0.37.0
  - rpds-py==2026.6.3
  - setuptools==82.0.1
  - smmap==5.0.3
  - starlette==1.3.1
  - streamlit-oauth==0.1.14
  - streamlit==1.60.0
  - tenacity==9.1.4
  - toml==0.10.2
  - uvicorn==0.51.0
  - watchdog==6.0.0
  - websockets==16.1.1
  - zstandard==0.25.0

# r-package(unknown source):
#  - compiler=3.6.3
#  - datasets=3.6.3
#  - graphics=3.6.3
#  - grDevices=3.6.3
#  - grid=3.6.3
#  - methods=3.6.3
#  - parallel=3.6.3
#  - splines=3.6.3
#  - stats=3.6.3
#  - stats4=3.6.3
#  - tcltk=3.6.3
#  - tools=3.6.3
#  - utils=3.6.3
