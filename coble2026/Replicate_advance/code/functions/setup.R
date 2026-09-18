# =============================================================================
# setup.R
# Libraries, global options, and command line argument parsing
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(Seurat)
library(SingleR)
library(celldex)
library(CellChat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(ggsci)
library(future)

# -----------------------------------------------------------------------------
# Global options
# -----------------------------------------------------------------------------
options(bitmapType = "cairo")
options(future.globals.maxSize = Inf)
future::plan("sequential")

# -----------------------------------------------------------------------------
# Cell type group definitions
# These should be updated if SingleR labels change
# -----------------------------------------------------------------------------
define_cell_groups <- function() {
  t_cell_types <- c(
    "Naive CD4 T cells", "Naive CD8 T cells", "Th1 cells",
    "Terminal effector CD8 T cells", "Follicular helper T cells",
    "T regulatory cells", "Th2 cells", "Central memory CD8 T cells",
    "MAIT cells", "Th1/Th17 cells", "Th17 cells",
    "Effector memory CD8 T cells", "Terminal effector CD4 T cells",
    "Vd2 gd T cells", "Non-Vd2 gd T cells"
  )

  b_cell_types <- c(
    "Naive B cells", "Non-switched memory B cells",
    "Switched memory B cells", "Plasmablasts", "Exhausted B cells"
  )

  myeloid_types <- c(
    "Intermediate monocytes", "Non classical monocytes",
    "Classical monocytes", "Myeloid dendritic cells",
    "Plasmacytoid dendritic cells", "Low-density basophils"
  )

  list(
    tcells  = t_cell_types,
    bcells  = b_cell_types,
    myeloid = myeloid_types
  )
}

# -----------------------------------------------------------------------------
# Argument parsing
# Usage: Rscript code/extend.R --celltype=tcells
# -----------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  valid_subsets <- c("all", "tcells", "bcells", "myeloid", "other")

  if (length(args) == 0) {
    stop(
      "Usage: Rscript code/extend.R --celltype=<subset>\n",
      "Options: ", paste(valid_subsets, collapse = ", ")
    )
  }

  # Parse --key=value style arguments
  parsed <- list()
  for (arg in args) {
    parts <- strsplit(arg, "=")[[1]]
    if (length(parts) == 2) {
      key <- sub("^--", "", parts[1])
      parsed[[key]] <- parts[2]
    }
  }

  cell_subset <- parsed[["celltype"]]

  if (is.null(cell_subset)) {
    stop("Missing --celltype argument. Options: ", paste(valid_subsets, collapse = ", "))
  }

  if (!cell_subset %in% valid_subsets) {
    stop(
      "Invalid subset '", cell_subset, "'. Must be one of: ",
      paste(valid_subsets, collapse = ", ")
    )
  }

  cat("[", format(Sys.time()), "] Cell subset selected:", cell_subset, "\n")
  cell_subset
}
