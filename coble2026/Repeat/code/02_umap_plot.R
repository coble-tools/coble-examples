### AUTHORS ###################################################################
# Name		    Date
# RAlcraft		2026-06-13
#
### DESCRIPTION ###############################################################
# Replication of ICR code for publication demonstration of `recreate``
### LIBRARIES #################################################################
library("BCN.general.utilities");
library("Seurat");
library("FC14.plotting.lib");
library("FC14.genesets.DB");
library("ggplot2");
library("BCN.genesets.utilities");
library("BCN.singlecell.utilities");
### FUNCTIONS #################################################################
rae_print_msg <- function(msg){
	cat("\n[", format(Sys.time()), "]", msg, "\n");
}
get_named_arg <- function(args, name){
    key <- paste0("--", name);
    key_eq <- paste0(key, "=");

    # Supports both --name value and --name=value syntaxes.
    idx <- match(key, args);
    if(!is.na(idx) && idx < length(args)){
        return(args[idx + 1]);
    }

    eq_idx <- startsWith(args, key_eq);
    if(any(eq_idx)){
        return(sub(key_eq, "", args[which(eq_idx)[1]]));
    }

    return(NULL);
}
### MAIN ######################################################################
################# input arguments ######################
args <- commandArgs(trailingOnly = TRUE);
append_tag <- get_named_arg(args, "tag");
input_tag <- get_named_arg(args, "tag");
outpath <- get_named_arg(args, "outpath");
interim_rds_path <- file.path(outpath, paste0("seu_CIDER_", input_tag, ".RDS"))
rae_print_msg(paste0("input seurat: ", interim_rds_path));
rae_print_msg(paste0("out path: ", outpath));
rae_print_msg(paste0("append tag: ", append_tag));
################# file path (updated paths) ######################
############# set parameters ##################
#set the seed for random events
seed <- set.seed(3235235);
regress.vars <- c(
	"percent.mt",
	"percent.rpl",
	"nFeature_RNA",
	"nCount_RNA"
);
umapdim <- 30;
################### load completed rds ################
rae_print_msg(paste0("reading...", interim_rds_path));
seu <- readRDS(interim_rds_path);
rae_print_msg("...loaded RDS");
################### run umap ################
rae_print_msg("run umap...");
seu <- RunUMAP(seu, dims = 1:umapdim, reduction = "pca", n.neighbors = 30, min.dist = 0.3);
rae_print_msg("...done");
################### make plot ################
rae_print_msg("make.reduction.plot...");
make.reduction.plot(
    seu,
    projection = "umap",
    mdata = "seurat_clusters_chr",
    xpos = "bottom",
    legend.cex = 1.2,
    point.size = 0.8,
    height = 7.7,
    width = 7.7,
    resolution = 600,
    font.face = 1,
    file.name = file.path(outpath, paste0("umap_CIDER_cluster_", append_tag, ".pdf"))
);
rae_print_msg(paste("...plotted to", file.path(outpath, paste0("umap_CIDER_cluster_", append_tag, ".pdf"))));
#### SAVE SESSION ##############################################################
BCN.general.utilities::save.session.info(
	file.name = paste0("_02session", append_tag, ".txt"),
	directory = outpath
);
################################################################

