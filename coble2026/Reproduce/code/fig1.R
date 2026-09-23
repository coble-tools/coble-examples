# from folder 06_velton
# https://github.com/veltenlab/EPI-clone/blob/main/figures/Figure1/Figure1.html

# Rscript Fig2_PanelB/code/fig1.R
#
# change history
# 2026-09-21	claude	added source data csv saves for figs 1c, 1g, 1i, 1j

source("code/source.R")

prefix <- paste0(getRversion(),"_")
suffix <- paste0("")
plot_folder <- paste0("results/")

# -- Write amplicon csv ---------------------------------------------
panel <- read.table('data/panel_info_dropout_pwm.tsv')
to_write <- panel[panel$DropoutUncutLKs>0.9, 'Type']
to_write <- to_write[to_write%in%c('HSC_high', 'MPP_high', 'MPPI_high', 'MPPII_high', 'WSH', 'IMR')]
to_write[to_write%in%c('HSC_high', 'MPP_high', 'MPPI_high', 'MPPII_high')] <- 'DMC'

# Differentiation uMAP (Figure 1c) ----------------------------------------------
# We simply load the data and visualize the integrated uMAP that shows three differentiation trajectories.
seurat_obj <- readRDS('data/larry_seurat.rds') #readRDS(url('https://figshare.com/ndownloader/files/42479346'))

###### Save the Fig 1c source data ###########
out_path_c <- "sourcedata/coble_fig1c_sourcedata.csv"
umap_coords_c <- Embeddings(seurat_obj, "umap")
out_df_c <- data.frame(umap_coords_c, CellType = seurat_obj$CellType)
colnames(out_df_c) <- c("umap_1", "umap_2", "CellType")
write.csv(out_df_c, out_path_c, row.names = TRUE)
cat("saved to", out_path_c, "\n")
###############################################

p <- DimPlot(seurat_obj,
  group.by = "CellType",reduction="umap", pt.size=1.2) +
  ggtitle("") +
  NoAxes() +
  NoLegend() +
  scale_color_manual(values = color_map)
ggsave(
  paste0(plot_folder, "Fig.c", suffix, ".pdf"), plot = p,
  width = 7,
  height = 7,
  dpi = 600,
)

# Surface protein expression per cell type (Figure 1f)---------------------------------
seurat_ab <- subset(seurat_obj,  ProcessingBatch%in%c('LK_LSK_stained', 'LARRY_mouse3', 'LARRY_mouse4'))

# Unsupervised uMAP (Figure 1g,h) ----------------------------
require(Seurat)
require(ggplot2)
require(ROCR)
require(fossil)
require(reshape2)
require(pheatmap)
require(GenomicRanges)
full_seurat <- readRDS('data/larry_seurat.rds') #(url('https://figshare.com/ndownloader/files/42479346'))
larry <- subset(full_seurat, Experiment == "LARRY main experiment")
usecpg <- rownames(larry)
larry <- ScaleData(larry, assay = "DNAm", features = usecpg, verbose = F)
larry <- RunPCA(larry, assay = "DNAm", features = usecpg, reduction.name = "pca", reduction.key = "PC_", npcs = 100, verbose = F)
larry <- RunUMAP(larry, reduction = "pca", dims = 1:50, verbose = F, seed.use = 1953, n.threads = 1)

###### Save the Fig 1g source data ###########
out_path_g <- "csv/coble_fig1g_sourcedata.csv"
umap_coords_g <- Embeddings(larry, "umap")
out_df_g <- data.frame(umap_coords_g, CellType = larry$CellType)
colnames(out_df_g) <- c("umap_1", "umap_2", "CellType")
write.csv(out_df_g, out_path_g, row.names = TRUE)
cat("saved to", out_path_g, "\n")
###############################################

celltypeColors <- c("HSC/MPP1"="maroon4",
                    "MPP2" ="darkgrey",
                    "MPP3" = "darkblue",
                    "MPP4" = "darkgreen",
                    "MEP"="#D5392C",
                    "EryP" = "#ED7950",
                    "MkP 1" = "#BDA8CB",
                    "MkP 2" = "#D36494",
                    "GMP" = "#C3C380",
                    "pre/pro-B" ="#94B1F9")
p <- DimPlot(larry, group.by = "CellType",reduction="umap", pt.size = 1.2) +
  ggtitle("") + NoAxes() + NoLegend() + scale_color_manual(values = celltypeColors)
ggsave(
  paste0(plot_folder, "Fig.g", suffix, ".pdf"), plot = p,
  width = 7,
  height = 7,
  dpi = 600
)

# Identification of static CpGs (Figure 1i) -------------------------------
#compute minimum pvalue for asociation with (any) protein
ks_test_safe <- function(x, y, ...) {
  result <- ks.test(x, y, ...)

  if (result$p.value == 0 && isFALSE(result$exact)) {
    # Recompute p-value in log space to recover the underflowed value
    n <- length(x) * length(y) / (length(x) + length(y))
    D <- as.numeric(result$statistic)
    arg <- sqrt(n) * D

    # KS asymptotic series: log(1 - CDF) directly
    # First term of the alternating series dominates: -2 * arg^2
    log_p <- log(2) + (-2 * arg^2)   # log(2 * exp(-2*arg^2))

    result$p.value <- exp(log_p)
    attr(result$p.value, "note") <- "recomputed via log-space approximation"
  }
  result
}
suppressWarnings({

pvals <- apply(larry@assays$DNAm@data,1, function(met) {
apply(larry@assays$AB@data, 1, function(prot) {
    use <- !is.na(prot)
    a <- prot[use][met[use]==1]
    b <- prot[use][met[use]==0]
    if (length(a) < 3 | length(b) < 3) return(1) else return(ks_test_safe(a,b)$p.value)
})
})
min_pval <- apply(pvals,2,min)

#establish bonferroni criterion
thr.protein.ass <- 1/(nrow(larry@assays$DNAm@data) * nrow(larry@assays$AB@data))

#determine average overall methylation level
avg_meth_rate <- apply(larry@assays$DNAm@data, 1, mean)
})

trueClone <- "LARRY"
ncells.bigClone <- 30
upper.thr.methrate <- 0.9
lower.thr.methrate <- 0.25
true_clone <- larry@meta.data[,trueClone]
for_prediction <- larry

for_prediction$use <- !is.na(true_clone)
for_prediction <- subset(for_prediction, use)
a <-table(for_prediction@meta.data[,trueClone])
use_for_prediction <- names(a)[a > ncells.bigClone]
for_prediction$use <- for_prediction@meta.data[,trueClone] %in% use_for_prediction
for_prediction <- subset(for_prediction,use)

cloneid <- factor(for_prediction@meta.data[,trueClone], levels = unique(for_prediction@meta.data[,trueClone]))
suppressWarnings({
  pvals_cloneass <- p.adjust(apply(for_prediction@assays[["DNAm"]]@data,1,function(met) {
  chisq.test(table(met,cloneid))$p.value
}),method = "bonferroni")
})

CpGSelection <- data.frame(CpG = names(pvals_cloneass), avg_meth_rate, min_pval, pvals_cloneass)
selected_not_protein <- names(min_pval)[min_pval > thr.protein.ass & avg_meth_rate < upper.thr.methrate & avg_meth_rate > lower.thr.methrate]
CpGSelection$Type <- ifelse(row.names(CpGSelection)%in%selected_not_protein, 'Static', 'Dynamic')
write.csv(CpGSelection, 'cpg_selection.csv')
panel <- read.table('data/panel_info_dropout_pwm.tsv',
                    sep='\t',
                    header=TRUE)

to_plot <- data.frame(PVal=min_pval,
                      AvgMeth=avg_meth_rate,
                      PValClone=pvals_cloneass)

###### Save the Fig 1i source data ###########
out_path_i <- "csv/coble_fig1i_sourcedata.csv"
write.csv(to_plot, out_path_i, row.names = TRUE)
cat("saved to", out_path_i, "\n")
###############################################

p <- ggplot(to_plot, aes(x = AvgMeth, y = log10(ifelse(PVal<1e-21, 1e-21, PVal)), color = -log10(PValClone+1e-50)))+
  geom_point(size=2.0, stroke=.5)+
  geom_hline(yintercept = log10(thr.protein.ass)) + geom_vline(xintercept = c(lower.thr.methrate,upper.thr.methrate)) +
  plot_theme + xlab("Average methylation") + ylab(ifelse(is.null(thr.protein.ass), "p value cell state association", "Association with surface\nprotein [log10]")) +
  scale_color_gradientn(colours = c("black","blue","red"), name = "-log10 p-val\nClone association")+
  scale_y_continuous(breaks=c(0, -7.5, -15), limits = c(-22,5))
ggsave(
  paste0(plot_folder, "Fig.i", suffix, ".pdf"),  plot = p,
  width = 7,
  height = 7,
  dpi = 600
)

# Check the localization of static vs. dynamic CpGs (Figure 1j) -------------------------------------------
cols_chrom <- c('active/weak promoter'='#ff678c',
          'Enhancer'='#ffdc64',
          'transcription'='#008c64',
          'Heterochromatin'='#6e1e8c',
          'weak enhancer'='#ffff00',
          'H3K9me3-repressed'='#787878',
          'Other'='#aaaaaa',
          'H3K9me3-repressed'='#f0f0f0')
panel_info <- read.table('data/panel_info_dropout_pwm.tsv',
                    sep='\t',
                    header=TRUE)
plot_dat <- panel_info[row.names(CpGSelection), 'ChromState', drop=FALSE]
plot_dat$ChromState <- c("active/weak promoter"='Other',
                        "heterochromatin"="Heterochromatin",
                        "poised promoter"='Other',
                        "strong enhancer"='Enhancer',
                        "transcription"='Other',
                        "weak enhancer"='Enhancer')[plot_dat$ChromState]
plot_dat$Type <- ifelse(row.names(plot_dat)%in%selected_not_protein, 'static', 'dynamic')
to_plot_dynamic <- plyr::count(plot_dat[plot_dat$Type=='dynamic', ])
to_plot_static <- plyr::count(plot_dat[plot_dat$Type=='static', ])
to_plot_dynamic$freq <- to_plot_dynamic$freq/sum(to_plot_dynamic$freq)
to_plot_static$freq <- to_plot_static$freq/sum(to_plot_static$freq)
to_plot <- rbind(to_plot_dynamic, to_plot_static)

###### Save the Fig 1j source data ###########
out_path_j <- "csv/coble_fig1j_sourcedata.csv"
write.csv(to_plot, out_path_j, row.names = FALSE)
cat("saved to", out_path_j, "\n")
###############################################

fisher.p <- fisher.test(table(plot_dat[plot_dat$ChromState!='Other', c('Type', 'ChromState')]))
p <- ggplot(to_plot, aes(x=Type, y=freq*100, fill=ChromState))+geom_bar(stat = 'identity')+
    plot_theme_legend+scale_fill_manual(values=cols_chrom)+
    xlab('CpG Class')+ylab('%of CpGs in chromatin state')
ggsave(
  paste0(plot_folder, "Fig.j", suffix, ".pdf"), plot = p,
  width = 7,
  height = 7,
  dpi = 600
)