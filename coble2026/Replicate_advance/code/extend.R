# =============================================================================
# extend.R
# Main orchestrator script for modernised PBMC scRNA-seq pipeline
# 
# Usage:
#   Rscript code/extend.R --celltype=tcells
#   Rscript code/extend.R --celltype=all
#
# Options: all, tcells, bcells, myeloid, other
# =============================================================================

# -----------------------------------------------------------------------------
# Source function files
# -----------------------------------------------------------------------------
PROG_DIR <- "code/functions"
source(file.path(PROG_DIR, "setup.R"))
source(file.path(PROG_DIR, "preprocess.R"))
source(file.path(PROG_DIR, "annotate.R"))
source(file.path(PROG_DIR, "subset.R"))
source(file.path(PROG_DIR, "visualise.R"))
source(file.path(PROG_DIR, "cellchat.R"))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
DATA_DIR <- "data"
RES_DIR  <- "results"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Parse command line argument
# -----------------------------------------------------------------------------
CELL_SUBSET <- parse_args()

# Results go in subset-specific directory
SUBSET_RES_DIR <- file.path(RES_DIR, CELL_SUBSET)
dir.create(SUBSET_RES_DIR, showWarnings = FALSE, recursive = TRUE)

# Cell group definitions
cell_groups <- define_cell_groups()

# =============================================================================
# STAGE 1: Full dataset pipeline
# Check for checkpoint - skip if already done
# =============================================================================
full_object_path <- file.path(RES_DIR, "pbmc_annotated.rds")

if (!file.exists(full_object_path)) {
  cat("[", format(Sys.time()), "] ### STAGE 1: Full dataset pipeline ###\n")

  # Load data
  pbmc <- load_data(DATA_DIR)

  # Preprocess - SCTransform, PCA, UMAP, clustering on ALL cells
  pbmc <- run_preprocess(pbmc, n_pcs = 30, resolution = 0.5, label = "full")

  # Elbow plot to guide PC selection
  plot_elbow(pbmc, RES_DIR, label = "full")

  # Load reference and annotate
  ref              <- load_reference(DATA_DIR)
  annotation       <- run_annotation(pbmc, ref)
  pbmc             <- annotation$seu
  singler_results  <- annotation$results

  # Diagnostic plots
  plot_singler_diagnostics(singler_results, RES_DIR)

  # Write cell type summary
  write_cell_types(pbmc, RES_DIR)

  # Summarise groups
  summarise_groups(pbmc, cell_groups)

  # Plot full dataset UMAP
  plot_umap_all(pbmc, RES_DIR)

  # Find markers on full dataset
  Idents(pbmc) <- "seurat_clusters"
  markers <- FindAllMarkers(
    pbmc,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.5
  )
  write.csv(markers, file.path(RES_DIR, "cluster_markers_full.csv"), row.names = FALSE)
  plot_markers(pbmc, markers, RES_DIR, label = "full")

  # Checkpoint save
  cat("[", format(Sys.time()), "] Saving full annotated object...\n")
  saveRDS(pbmc, full_object_path)
  cat("[", format(Sys.time()), "] Checkpoint saved:", full_object_path, "\n")

} else {
  cat("[", format(Sys.time()), "] ### STAGE 1: Loading from checkpoint ###\n")
  pbmc <- readRDS(full_object_path)
  cat("[", format(Sys.time()), "] Loaded:", ncol(pbmc), "cells\n")
}

# =============================================================================
# STAGE 2: Subset pipeline
# =============================================================================
cat("[", format(Sys.time()), "] ### STAGE 2: Subset pipeline [", CELL_SUBSET, "] ###\n")

# Subset cells based on command line argument
pbmc_sub <- subset_cells(pbmc, CELL_SUBSET, cell_groups)

# Rerun SCTransform, PCA, UMAP, clustering on subset only
# Fewer PCs for subset - less cell type diversity
n_pcs_sub <- ifelse(CELL_SUBSET == "all", 30, 20)
pbmc_sub  <- run_preprocess(
  pbmc_sub,
  n_pcs      = n_pcs_sub,
  resolution = 0.3,     # lower resolution for subset - fewer expected clusters
  label      = CELL_SUBSET
)

# Elbow plot for subset
plot_elbow(pbmc_sub, SUBSET_RES_DIR, label = CELL_SUBSET)

# Plot subset UMAP - this will look much better than the full dataset UMAP
plot_umap_subset(pbmc_sub, CELL_SUBSET, SUBSET_RES_DIR)

# Marker genes for subset
Idents(pbmc_sub) <- "seurat_clusters"
markers_sub <- FindAllMarkers(
  pbmc_sub,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5
)
write.csv(
  markers_sub,
  file.path(SUBSET_RES_DIR, paste0("cluster_markers_", CELL_SUBSET, ".csv")),
  row.names = FALSE
)
plot_markers(pbmc_sub, markers_sub, SUBSET_RES_DIR, label = CELL_SUBSET)

# Save subset object
saveRDS(pbmc_sub, file.path(SUBSET_RES_DIR, paste0("pbmc_", CELL_SUBSET, ".rds")))

# =============================================================================
# STAGE 3: CellChat on stratified resample of subset
# =============================================================================
cat("[", format(Sys.time()), "] ### STAGE 3: CellChat [", CELL_SUBSET, "] ###\n")

pbmc_cc  <- resample_cells(pbmc_sub, prop = 0.15, min_cells = 10, max_cells = 500)
cellchat <- run_cellchat(pbmc_cc, SUBSET_RES_DIR)
plot_cellchat(cellchat, SUBSET_RES_DIR)

# =============================================================================
# Done
# =============================================================================
cat("[", format(Sys.time()), "] ### Pipeline complete for subset:", CELL_SUBSET, "###\n")
cat("[", format(Sys.time()), "] Results saved to:", SUBSET_RES_DIR, "\n")
si <- sessionInfo()
print(si)
