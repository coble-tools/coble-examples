# Repeatability of the results
`same inputs, same team`  
We should be able to exactly repeat our own results in our lab on our hardware with our internal code even where not submitted to public realm.

*Zhang, Y., Moughari, F.A., Mavrommati, I. et al. Co-occurrence of transcriptionally distinct persister cell states underpins neoadjuvant therapy resistance in triple‑negative breast cancer. Genome Med 18, 65 (2026). https://doi.org/10.1186/s13073-026-01643-9*  

Publication link: https://link.springer.com/article/10.1186/s13073-026-01643-9

Target Figure: Figure 1f  
*f. UMAP for G1 epithelial cells coloured by MCs.*  
![alt text](figures/1f.png)

---  
This analysis takes substantial memory for the sctransform stage and needs to be run on a high-memory machine. The sct phase uses 23,789,461,504

The analysis scripts are split up so that the sctransform stage can be run separately on a high-memory machine and the rest of the analysis can be run on a regular machine. The results of sctransform are saved as an RDS file which can be loaded in the second stage of the analysis if preferred, avoiding the need to run sctransform again. 

S:\ralcraft\Syed\SingleCell\tumori\08_tnbc_atlas\tnbc_pdx_residual\code\snrseq\plots
Original source of results

` Rscript 01_sct.R`  
[ 2026-06-10 22:40:19 ] reading...data/seu_CIDER_rna.RDS 
[ 2026-06-10 22:41:05 ] run sctransform... 
[ 2026-06-10 22:52:50 ] saving interim rds... 
[ 2026-06-10 23:12:04 ] ...saved to data/seu_CIDER_sct.RDS 

` Rscript 02_pca.R`  
[ 2026-06-10 23:20:28 ] reading...data/seu_CIDER_sct.rds 
[ 2026-06-10 23:22:07 ] run pca... 
[ 2026-06-10 23:22:22 ] saving interim rds... 
[ 2026-06-10 23:41:32 ] ...saved to data/seu_CIDER_pca.rds 

`Rscript 03_umap.R`  
