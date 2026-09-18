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
library("BCN.singlecell.utilities"); #### newly added line
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
in_path <- get_named_arg(args, "seuratin");
append_tag <- get_named_arg(args, "tag");
outpath <- get_named_arg(args, "outpath");
interim_rds_path <- file.path(outpath, paste0("seu_CIDER_", append_tag, ".RDS"))
rae_print_msg(paste0("input seurat: ", in_path));
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
rae_print_msg(paste0("reading...", in_path));
seu <- readRDS(in_path);
rae_print_msg("...loaded RDS");
################### cut calced results to re-run ################
# Remove SCTransform results
seu <- DietSeurat(seu, assays = "RNA")  # keeps only raw RNA, drops SCT assay entirely
################### run sctransform ################
rae_print_msg("run sctransform...");
seu <- SCTransform(seu, vst.flavor = "v2", vars.to.regress = regress.vars, verbose = FALSE);
rae_print_msg("...done");
################### run pca ################
rae_print_msg("run pca...");
seu <- RunPCA(seu, verbose = FALSE);
rae_print_msg("...done");
################### cut calced results to save ################
seu_slim <- DietSeurat(seu, assays = c("SCT"), dimreducs = "pca", graphs    = c("SCT_nn", "SCT_snn"))
################### save interim rds ################
rae_print_msg("saving interim rds...");
saveRDS(seu, interim_rds_path)
rae_print_msg(paste("...saved to", interim_rds_path));
#### SAVE SESSION ##############################################################
BCN.general.utilities::save.session.info(
	file.name = paste0("_01session", append_tag, ".txt"),
	directory = outpath
);
################################################################

