#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(plyr)
library(dplyr)
library(data.table)
library(pheatmap)

options(bitmapType='cairo')

# -------------------------------------
# 1. normalize by RNA content (umi counts) and select the top 1000 most variable genes
# 2. use top 1000 variable genes for PCA
# 3. generate 2-D tSNE embedding
# 4. use k-means clustering to specify populations
# 5. identify sub-clusters in cluster #9
# -------------------------------------
# -------------------------------------
# specify paths and load functions
# -------------------------------------
cat("[", format(Sys.time()), "] ### specify paths to data, program and results directories in the code and source the utility functions\n")
DATA_DIR <- "data"        # SPECIFY HERE
PROG_DIR <- "code"     # SPECIFY HERE
RES_DIR  <- "results"      # SPECIFY HERE
source(file.path(PROG_DIR,'util.R'))
# -------------------------------------
# Load data to take into account diffferent runs and checkpointing
# -------------------------------------
# 11 purified PBMC data
cat("[", format(Sys.time()), "] ### load 11 purified PBMC the local run for purified PBMC data\n")
pure_11 <- readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))

# ------------------------------------------------------------
# load 68k PBMC data, and meta-data
# ------------------------------------------------------------
cat("[", format(Sys.time()), "] ### load 68k PBMC data, 11 purified PBMC data and meta-data\n")
pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
all_data <- pbmc_68k$all_data
purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$ens_genes)
# --------------------------------------------------------------------------------------
# normalize by RNA content (umi counts) and select the top 1000 most variable genes
# --------------------------------------------------------------------------------------
cat("[", format(Sys.time()), "] ### normalize by RNA content (umi counts) and select the top 1000 most variable genes\n")
m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)
m_n<-l$m
df<-.get_variable_gene(m_n)
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
# --------------------------------------------------
# plot dispersion vs. mean for the genes
# this produces Supp. Fig. 5c in the manuscript
# --------------------------------------------------
cat("[", format(Sys.time()), "] ### plot dispersion vs. mean for the genes\n")
p <- ggplot(df,aes(mean,dispersion,col=used))+geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
  scale_color_manual(values=c("grey","black"))+theme_classic()
ggsave(
  file.path(RES_DIR, "Fig5c.pdf"), p,
  width=7,
  height=7,
  dpi=600
)
# --------------------------------------------
# use top 1000 variable genes for PCA
# --------------------------------------------
cat("[", format(Sys.time()), "] ### use top 1000 variable genes for PCA\n")
set.seed(0)
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)
# --------------------------------------------
# generate 2-D tSNE embedding
# this step may take a long time
# --------------------------------------------
cat("[", format(Sys.time()), "] ### generate 2-D tSNE embedding\n")
tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
tdf_n_1000<-data.frame(tsne_n_1000$Y)
# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations
# this produces Fig. 3j in the manuscript
# ---------------------------------------------------------------------------------------------------------------------------
cat("[", format(Sys.time()), "] ### assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations\n")
m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11)
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])
tdf_n_1000$cls_id<-cls_id
# adjust ordering of cells for plotting aesthetics
tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
p <- ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()
ggsave(
  file.path(RES_DIR, "Fig3j.pdf"), p,
  width=9,
  height=7,
  dpi=600
)
# ----------------------------------
cat("[", format(Sys.time()), "] ### done\n")
si <- sessionInfo()
cat("[", format(Sys.time()), "] Session info:\n")
print(si)
# ----------------------------------



