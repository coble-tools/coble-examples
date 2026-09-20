#####################################################
# COBLE:Reproducible environment: carbine, (c) ICR 2026
#####################################################
coble:
  - environment: r-443-conda
channels:
  - defaults  
  - bioconda
  - conda-forge
flags:  
  - compile: tools=true  
languages:  
  - r-base=4.4.3
  - python=3.12
flags:    
  - compile: system=True
  - export: QT_QPA_PLATFORM=offscreen  
  - export: OTEL_SDK_DISABLED=true
  - export: R_OTEL_DISABLED=true
conda:  
  - arviz