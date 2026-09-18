# =============================================================================
# visualise.R
# All plotting functions
# =============================================================================

# -----------------------------------------------------------------------------
# Curated colour palette for T cells
# -----------------------------------------------------------------------------
t_cell_colours <- c(
  # CD4 helper -
  "Naive CD4 T cells"             = "#90d871",
  "Th1 cells"                     = "#FF7F00",
  "Th2 cells"                     = "#FFFF33",
  "Th17 cells"                    = "#F781BF",
  "Th1/Th17 cells"                = "#A65628",
  "Follicular helper T cells"     = "#984EA3",
  "T regulatory cells"            = "#418e3e",
  "Terminal effector CD4 T cells" = "#377EB8",
  # CD8 -
  "Naive CD8 T cells"             = "#00CED1",
  "Central memory CD8 T cells"    = "#1E90FF",
  "Effector memory CD8 T cells"   = "#000080",
  "Terminal effector CD8 T cells" = "#00FFFF",
  # Innate-like -
  "MAIT cells"                    = "#ed331e",
  "Vd2 gd T cells"                = "#98DF8A",
  "Non-Vd2 gd T cells"            = "#17BECF"
)

# -----------------------------------------------------------------------------
# Generate a palette for any set of cell types
# Uses Polychrome for maximum discrimination
# -----------------------------------------------------------------------------
make_palette <- function(cell_types) {
  n <- length(cell_types)
  cols <- ggsci::pal_d3("category20")(20)
  if (n > 20) cols <- rep(cols, ceiling(n / 20))  # recycle if somehow >20 types
  setNames(cols[1:n], cell_types)
}

# -----------------------------------------------------------------------------
# UMAP coloured by SingleR cell type - full dataset overview
# -----------------------------------------------------------------------------
plot_umap_all <- function(seu, res_dir) {
  cat("[", format(Sys.time()), "] Plotting full dataset UMAP...\n")

  cell_types <- unique(seu$singler_pruned[!is.na(seu$singler_pruned)])
  palette    <- make_palette(cell_types)

  p <- DimPlot(
    seu,
    reduction = "umap",
    group.by  = "singler_pruned",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = 0.3
  ) +
    scale_colour_manual(values = palette, na.value = "grey80") +
    ggtitle("All cells - SingleR annotation") +
    theme_classic() +
    theme(legend.text = element_text(size = 8))

  ggsave(
    file.path(res_dir, "umap_all_celltypes.pdf"), p,
    width = 12,
    height = 8,
    dpi = 600
  )
}

# -----------------------------------------------------------------------------
# UMAP coloured by SingleR cell type - subset
# Uses curated T cell palette where available, falls back to Polychrome
# -----------------------------------------------------------------------------
plot_umap_subset <- function(seu, cell_subset, res_dir) {
  cat("[", format(Sys.time()), "] Plotting subset UMAP [", cell_subset, "]...\n")

  cell_types <- unique(seu$singler_pruned[!is.na(seu$singler_pruned)])

  # Use curated T cell palette if available, otherwise generate one
  if (cell_subset == "tcells" &&
      all(cell_types %in% names(t_cell_colours))) {
    palette <- t_cell_colours[cell_types]
  } else {
    palette <- make_palette(cell_types)
  }

  p <- DimPlot(
    seu,
    reduction = "umap",
    group.by  = "singler_pruned",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = 0.5
  ) +
    scale_colour_manual(values = palette, na.value = "grey80") +
    ggtitle(paste0(cell_subset, " - recomputed UMAP")) +
    theme_classic() +
    theme(legend.text = element_text(size = 8))

  ggsave(
    file.path(res_dir, paste0("umap_", cell_subset, ".pdf")),
    p,
    width = 12,
    height = 8,
    dpi = 600
  )
}

# -----------------------------------------------------------------------------
# # Elbow plot for PC selection
# # -----------------------------------------------------------------------------
# plot_elbow <- function(seu, res_dir, label = "") {
#   p <- ElbowPlot(seu, ndims = 50) +
#     ggtitle(paste0("Elbow plot", ifelse(label != "", paste0(" - ", label), ""))) +
#     theme_classic()
#   ggsave(
#     file.path(res_dir, paste0("elbow_", label, ".pdf")), p,
#     width = 7,
#     height = 7,
#     dpi = 600
#   )
# }

# -----------------------------------------------------------------------------
# SingleR diagnostic plots
# -----------------------------------------------------------------------------
plot_singler_diagnostics <- function(singler_results, res_dir) {
  cat("[", format(Sys.time()), "] Plotting SingleR diagnostics...\n")

  png(file.path(res_dir, "singler_score_heatmap.png"), width = 1400, height = 800)
  plotScoreHeatmap(singler_results)
  dev.off()

  png(file.path(res_dir, "singler_delta_distribution.png"), width = 1200, height = 600)
  plotDeltaDistribution(singler_results, ncol = 4)
  dev.off()
}

# # -----------------------------------------------------------------------------
# # Marker gene dot plot
# # -----------------------------------------------------------------------------
# plot_markers <- function(seu, markers, res_dir, label = "") {
#   cat("[", format(Sys.time()), "] Plotting marker genes...\n")

#   top5 <- markers %>%
#     group_by(cluster) %>%
#     slice_max(order_by = avg_log2FC, n = 5)

#   p <- DotPlot(seu, features = unique(top5$gene)) +
#     RotatedAxis() +
#     ggtitle(paste0("Top 5 markers per cluster", ifelse(label != "", paste0(" - ", label), ""))) +
#     theme(axis.text.x = element_text(size = 7))

#   ggsave(
#     file.path(res_dir, paste0("dotplot_markers_", label, ".pdf")),
#     p,
#     width = 18,
#     height = 6,
#     dpi = 600
#   )
}
