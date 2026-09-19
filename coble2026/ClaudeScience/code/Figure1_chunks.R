# =============================================================================
# Bach et al. 2021, Nat Commun 12:1502   --  FIGURE 1
# https://doi.org/10.1038/s41467-021-21783-3
#
# VERBATIM R chunks from MarioniLab/Tumorigenesis2018
#   src/DownstreamAnalysis/Figure1.Rmd
# concatenated in document order. This file is sourced by Figure1.R with
# print.eval=TRUE so that top-level expressions are auto-printed exactly as
# the R notebook did (ggsave() calls with no plot= argument depend on it).
#
# Deviations from the .Rmd, and nothing else:
#   1. the two hard-coded paths in the setup chunk read SCE_PATH / OUT_DIR,
#      set by the driver, instead of the author's working-tree paths;
#   2. source("../functions.R") -> source(FUNCTIONS_R); the same file is
#      vendored unchanged alongside this one;
#   3. one added ggsave for p0 (the faceted UMAP), which the notebook only
#      displayed inline and never wrote to disk;
#   4. the gganimate block building the repo animation fancyUMAP.gif is
#      omitted: it is not a panel of Figure 1 and needs a gifski/av backend
#      that is not part of the published analysis.
# =============================================================================

# ----- chunk 1/5  setup and data load ------------------------------------
library(scran)
library(scater)
library(ggplot2)
library(cowplot)
library(dplyr)
library(edgeR)
library(Matrix)
library(ggrepel)
library(ggrastr)
library(viridis)
theme_set(theme_cowplot())

sce.full <- readRDS(SCE_PATH)

sce <- sce.full[,sce.full$Experiment=="Tumorigenesis" & sce.full$Condition!="WTOld"]
pD <- data.frame(colData(sce))

outdr <- OUT_DIR


# ----- chunk 2/5  Fig 1 UMAP faceted by tumour-time bin ------------------
library(viridis)
p0 <- ggplot(pD, aes(x=UMAP1, y=UMAP2)) +
    geom_point(data=pD[,c("UMAP1","UMAP2")], color="grey80",size=.1,alpha=.2) +
    geom_point(size=.1, aes(color=Groups)) +
    facet_wrap(~ptimeBin,nrow=2) +
    scale_color_manual(values=levels(pD$GroupColors)) +
    guides(colour=guide_legend(title.position="top", override.aes = list(size=3))) +
    theme_void() +
    theme(legend.position="bottom",
	  legend.direction="horizontal",
	  strip.text=element_text(size=14))
p0

# [added] the notebook only displayed p0; write it out as well
ggsave(paste0(outdr,"UMAP_ptimeBin_Groups.pdf"), p0, width=9, height=6)


# ----- chunk 3/5  compartment composition barplot ------------------------
sce.sub <- sce.full[,sce.full$Condition %in% c("1","2","3","4","5","WTYoung","WTOld")]
pD <- data.frame(colData(sce.sub))
pD.add <- data.frame(colData(sce.sub))
pD.add$Condition <- "All"
pD <- droplevels(rbind(pD,pD.add))

sumry <- group_by(pD, Condition, Groups) %>%
    summarize(n=n()) %>%
    mutate(frac=n/sum(n)) %>%
    ungroup()

comb <- table(sumry$Groups,sumry$Condition) 
missing.comb <- which(comb==0, arr.ind=TRUE)
if (nrow(missing.comb)>0) {
    add.missing <- data.frame("Groups"=rownames(missing.comb),
		      "Condition"=colnames(comb)[missing.comb[,"col"]],
		      "frac"=0,
		      "n"=0)
    sumry <- rbind(sumry,add.missing)
}

sumry$Groups <- factor(sumry$Groups, levels=c("TumorEpithel","Basal","Luminal",
					      "Myeloid","Lymphoid","Fibroblast",
					      "Stroma"))

sumry$Colors <- plyr::mapvalues(as.character(sumry$Groups),
				as.character(unique(pD$Groups)),
				as.character(unique(pD$GroupColors)))

sumry <- arrange(sumry,Groups)

barplt <- ggplot(sumry, aes(x=Condition, y=frac, fill=Groups)) +
    geom_bar(stat="identity",color="black") +
    scale_fill_manual(values=unique(sumry$Colors)) +
    xlab("Condition") +
    ylab("Fraction of Cells") 
ggsave(paste0(outdr,"Barplot_Compartments.svg"))


# ----- chunk 4/5  cell-type UMAPs, rasterised, with labels ---------------
dpi <- 800

ggplot(pD, aes(x=UMAP1, y=UMAP2, fill=CellTypesFinal)) +
    geom_point_rast(size=1.8,pch=21,stroke=.2,color="black") +
    scale_fill_manual(values=levels(pD$Color)) +
    theme_void() +
    theme(legend.position="none")

ggsave(paste0(outdr,"CellTypesUMAP.png"),dpi=dpi)

pD.full <- data.frame(colData(sce.full))
ggplot(pD.full, aes(x=UMAP1, y=UMAP2, fill=CellTypesFinal)) +
    geom_point_rast(size=1.8,pch=21,stroke=.2,color="black") +
    scale_fill_manual(values=levels(pD$Color)) +
    theme_void() +
    theme(legend.position="none")
ggsave(paste0(outdr,"CellTypesUMAP_FullData.png"),dpi=dpi)

# ggplot(pD, aes(x=UMAP1, y=UMAP2, fill=Batch)) +
#     geom_point_rast(size=1.8,pch=21,stroke=.2,color="black") +
#     theme_void() +
#     theme(legend.position="none")


# With Annotation
X <- model.matrix(~0+pD$CellTypesFinal)
colnames(X) <- levels(pD$CellTypesFinal)
X <- t(t(X)/colSums(X))
celltyp_ctrs <- data.frame(crossprod(X, as.matrix(pD[,c("UMAP1","UMAP2")])))
celltyp_ctrs$CellType <- rownames(celltyp_ctrs)

source(FUNCTIONS_R)
celltyp_ctrs$CellType <- renameForPlot(celltyp_ctrs$CellType)

ggplot(pD, aes(x=UMAP1, y=UMAP2, color=CellTypesFinal)) +
    geom_point_rast(size=.3, dpi=dpi) +
    scale_color_manual(values=levels(pD$Color)) +
    geom_text_repel(data=celltyp_ctrs,aes(label=CellType,color=NULL)) +
    theme_void() +
    theme(legend.position="none")
ggsave(paste0(outdr,"CellTypesUMAP_WithLabels.pdf"),dpi=dpi)


# ----- chunk 5/5  per-condition split UMAPs ------------------------------
pltUMAP <- function(pD.sub,pD) {
    p <- ggplot(pD, aes(x=UMAP1, y=UMAP2)) +
	geom_point_rast(color="grey80",size=1.8,dpi=dpi) +
	geom_point_rast(data=pD.sub, aes(fill=Groups), size=1.8,dpi=900,pch=21,stroke=.2) +
	scale_fill_manual(values=levels(pD.sub$GroupColors)) +
	theme_void() +
	theme(legend.position="none")
    return(p)
}


pDWt <- data.frame(colData(sce.full[,sce.full$Condition %in% c("WTYoung","WTOld")]))

set.seed(42)
# Select each pD individually and downsample
pD1 <- pD[pD$Condition==1,]
pD1 <- pD1[sample(nrow(pD1),12000),]

pD2 <- pD[pD$Condition==2,]
pD2 <- pD2[sample(nrow(pD2),12000),]

pD3 <- pD[pD$Condition==3,]
pD3 <- pD3[sample(nrow(pD3),12000),]

pD4 <- pD[pD$Condition==4,]
pD4 <- pD4[sample(nrow(pD4),12000),]

pD5 <- pD[pD$Condition==5,]

pDWtO <- pDWt[pDWt$Condition == "WTOld",]
pDWtO <- pDWtO[sample(nrow(pDWtO),12000),]
pDWtY <- pDWt[pDWt$Condition == "WTYoung",]
pDWtY <- pDWtY[sample(nrow(pDWtY),12000),]


ggsave(paste0(outdr,"UMAP_AllTumorigenesis.png"),pltUMAP(pD,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_1.png"),pltUMAP(pD1,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_2.png"),pltUMAP(pD2,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_3.png"),pltUMAP(pD3,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_4.png"),pltUMAP(pD4,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_5.png"),pltUMAP(pD5,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_WT_Young.png"),pltUMAP(pDWtY,pD),dpi=dpi)
ggsave(paste0(outdr,"UMAP_Bin_WT_Old.png"),pltUMAP(pDWtO,pD),dpi=dpi)

# Save one for later figure w all cells
pDall <- data.frame(colData(sce.full))

p <- ggplot(pDall, aes(x=UMAP1, y=UMAP2)) +
    geom_point_rast(aes(color=Groups), size=.3,dpi=900) +
    scale_color_manual(values=levels(pDall$GroupColors)) +
    theme_void() +
    theme(legend.position="none")

ggsave(paste0(outdr,"UMAP_TumorAndPregnancyCellGroups.png"), p, dpi=dpi)

library(RColorBrewer)
col.t <- brewer.pal(n=6,"PuRd")[5]
col.p <- brewer.pal(n=5,"GnBu")[5]

pDall <- group_by(pDall, Condition) %>%
	sample_n(3000) %>%
	ungroup()
set.seed(123)
pDall <- pDall[sample(1:nrow(pDall),nrow(pDall)),]

p <- ggplot(pDall, aes(x=UMAP1, y=UMAP2)) +
    geom_point_rast(aes(fill=Experiment), size=1.8,dpi=900,pch=21,stroke=.2) +
    scale_fill_manual(values=c(col.p,col.t)) +
    theme_void() +
    theme(legend.position="none")

ggsave(paste0(outdr,"UMAP_TumorAndPregnancyExperiment.png"), p, dpi=dpi)

