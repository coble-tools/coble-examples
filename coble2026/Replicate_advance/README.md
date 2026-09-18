# Generalisability of the results
`extend to new settings`  
We should be able to build on results with new analyses.

*Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017). https://doi.org/10.1038/ncomms14049*  

Publication link: https://www.nature.com/articles/ncomms14049




To run with singularity:
```bash
srun --pty --mem=8042 --cpus-per-task=8 --time=08:00:00 bash
singularity shell cbl/cbl-zheng02.sif
bash extend.sh
```

sbatch << 'EOF'
#!/bin/bash
#SBATCH --job-name=zheng
#SBATCH --output=_logzheng.out
#SBATCH --error=_logzheng.err
#SBATCH --ntasks=1
#SBATCH --time=2:30:00
#SBATCH --mem-per-cpu=16084
#SBATCH --cpus-per-task=16
#SBATCH --partition=smp
singularity exec cbl/cbl-zheng02.sif Rscript code/extend.R
EOF