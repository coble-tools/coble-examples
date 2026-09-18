# =============================================================================
# annotate.R
# SingleR cell type annotation functions
# =============================================================================

# -----------------------------------------------------------------------------
# Load or download Monaco immune reference
# Caches locally to avoid re-downloading on HPC compute nodes
# -----------------------------------------------------------------------------
load_reference <- function(data_dir) {
  ref_path <- file.path(data_dir, "monaco_immune_ref.rds")

  if (!file.exists(ref_path)) {
    cat("[", format(Sys.time()), "] Downloading Monaco reference...\n")
    ref <- celldex::MonacoImmuneData()
    saveRDS(ref, ref_path)
    cat("[", format(Sys.time()), "] Reference cached to:", ref_path, "\n")
  } else {
    cat("[", format(Sys.time()), "] Loading cached Monaco reference...\n")
    ref <- readRDS(ref_path)
  }

  ref
}

# -----------------------------------------------------------------------------
# Run SingleR annotation and attach labels to Seurat object
# -----------------------------------------------------------------------------
run_annotation <- function(seu, ref) {
  cat("[", format(Sys.time()), "] Running SingleR annotation...\n")

  sce_counts <- GetAssayData(seu, assay = "SCT", layer = "data")

  singler_results <- SingleR(
    test   = sce_counts,
    ref    = ref,
    labels = ref$label.fine
  )

  seu$singler_label  <- singler_results$labels
  seu$singler_pruned <- singler_results$pruned.labels  # NA = low confidence

  # Verify annotation worked
  if (!"singler_pruned" %in% colnames(seu@meta.data)) {
    stop("SingleR annotation failed - singler_pruned not found in metadata")
  }

  cat("[", format(Sys.time()), "] Annotation complete. Distribution:\n")
  print(table(seu$singler_pruned, useNA = "always"))

  list(seu = seu, results = singler_results)
}

# -----------------------------------------------------------------------------
# Write cell type counts to file and console
# -----------------------------------------------------------------------------
write_cell_types <- function(seu, res_dir) {
  cell_type_table <- as.data.frame(
    table(seu$singler_pruned, useNA = "always")
  )
  colnames(cell_type_table) <- c("cell_type", "count")
  cell_type_table <- cell_type_table[order(-cell_type_table$count), ]

  write.csv(
    cell_type_table,
    file.path(res_dir, "cell_types.csv"),
    row.names = FALSE
  )

  cat("\nCell type counts:\n")
  print(cell_type_table)
}
