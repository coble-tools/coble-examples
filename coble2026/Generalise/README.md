# Generalisability of the results
`extend to new settings`  
Take existing environment and analysis and see that the results generalise by extending it to a new dataset, and gain any additional insight.

Original dataset:  
*Bach, K., Pensa, S., Zarocsinceva, M. et al. Time-resolved single-cell analysis of Brca1 associated mammary tumourigenesis reveals aberrant differentiation of luminal progenitors. Nat Commun 12, 1502 (2021). https://doi.org/10.1038/s41467-021-21783-3*  
[Publication link: https://www.nature.com/articles/ncomms14049](https://www.nature.com/articles/s41467-021-21783-3)

Generalise-to dataset:  
*Pal B, Chen Y, Milevskiy MJG, Vaillant F, Prokopuk L, Dawson CA, Capaldo BD, Song X, Jackling F, Timpson P, Lindeman GJ, Smyth GK, Visvader JE. Single cell transcriptome atlas of mouse mammary epithelial cells across development. Breast Cancer Res. 2021 Jun 29;23(1):69. doi: 10.1186/s13058-021-01445-4. PMID: 34187545; PMCID: PMC8243869.*
[Publication link: ncbi-PMC8243869](https://pmc.ncbi.nlm.nih.gov/articles/PMC8243869/)

---  

## Environment
INSTALL container and data
```bash
wget https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/Tumorigenesis.sif
wget -P data https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/BRCA1_SCE.rds
# wget -r -np -nH --cut-dirs=6 ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10043/ # Not needed if not rerunning from the beginning this is fastq files
```

## Generalise
To open the container do: `singularity shell Tumorigenesis.sif` or alternatively run the analysies directly with the container.
Follow the steps here in bash. If memory is a problem there is a commented out slurm option to use more memory for the relevant stages:
```
# Produce the exact replication of bach's own figure 3i
sbatch -p smp --mem=64G --cpus-per-task=4 --wrap="singularity exec Tumorigenesis.sif Rscript code/FigSupp_strain.R"

# Download the Pal data
bash code/download_GSE164017.sh

# Prepare the Pal data into a processed rds object
singularity exec Tumorigenesis.sif Rscript code/prep_GSE164017.R
### sbatch -p smp --mem=64G --cpus-per-task=4 --wrap="singularity exec Tumorigenesis.sif Rscript code/prep_GSE164017.R"

# Replicate the coble publication figure 2i
singularity exec Tumorigenesis.sif Rscript code/Fig3i.R
### sbatch -p smp --mem=64G --cpus-per-task=4 --wrap="singularity exec Tumorigenesis.sif Rscript code/Fig3i.R"

# Show the gene expression in mouse strains for supplementary data
singularity exec Tumorigenesis.sif Rscript code/FigSupp_bach.R
### sbatch -p smp --mem=64G --cpus-per-task=4 --wrap="singularity exec Tumorigenesis.sif Rscript code/FigSupp_bach.R"
```

