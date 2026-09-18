# =============================================================================
# preprocess.R
# SCTransform, PCA, UMAP, and clustering functions
# Designed to be called on both full dataset and subsets
# =============================================================================

# -----------------------------------------------------------------------------
# Load raw data from original 10x RDS format
# -----------------------------------------------------------------------------
load_data <- function(data_dir) {
  cat("[", format(Sys.time()), "] Loading data...\n")

  pbmc_68k    <- readRDS(file.path(data_dir, "pbmc68k_data.rds"))
  m           <- pbmc_68k$all_data[[1]]$hg19$mat
  gene_symbols <- pbmc_68k$all_data[[1]]$hg19$gene_symbols

  # Seurat requires genes x cells - transpose if needed
  if (nrow(m) > ncol(m)) {
    cat("[", format(Sys.time()), "] Transposing matrix to genes x cells...\n")
    m <- t(m)
  }

  rownames(m) <- make.unique(gene_symbols)

  pbmc <- CreateSeuratObject(
    counts       = m,
    min.cells    = 3,
    min.features = 200
  )

  cat("[", format(Sys.time()), "] Cells loaded:", ncol(pbmc), "\n")
  pbmc
}

# -----------------------------------------------------------------------------
# Run SCTransform, PCA, UMAP, and clustering
# Called on both full dataset and subsets
# n_pcs: number of PCs to use for UMAP/clustering
# resolution: Leiden clustering resolution
# -----------------------------------------------------------------------------
run_preprocess <- function(seu, n_pcs = 30, resolution = 0.5, label = "full") {
  cat("[", format(Sys.time()), "] Running SCTransform [", label, "]...\n")
  seu <- SCTransform(seu, variable.features.n = 3000, verbose = FALSE)

  cat("[", format(Sys.time()), "] Running PCA [", label, "]...\n")
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

  cat("[", format(Sys.time()), "] Running UMAP [", label, "]...\n")
  seu <- RunUMAP(seu, dims = 1:n_pcs, seed.use = 42)

  cat("[", format(Sys.time()), "] Clustering [", label, "]...\n")
  seu <- FindNeighbors(seu, dims = 1:n_pcs)
  seu <- FindClusters(seu, resolution = resolution, algorithm = 4, random.seed = 1)

  cat("[", format(Sys.time()), "] Clusters found:", nlevels(seu$seurat_clusters), "\n")
  seu
}
