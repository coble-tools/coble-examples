# =============================================================================
# cellchat.R
# CellChat cell-cell communication inference and visualisation
# =============================================================================

# -----------------------------------------------------------------------------
# Run full CellChat inference pipeline
# seu: subsetted and resampled Seurat object
# -----------------------------------------------------------------------------
run_cellchat <- function(seu, res_dir) {
  cat("[", format(Sys.time()), "] Setting up CellChat...\n")

  # Remove NA labels
  seu <- subset(
    seu,
    cells = colnames(seu)[!is.na(seu$singler_pruned)]
  )

  data_input <- GetAssayData(seu, assay = "SCT", layer = "data")

  cellchat <- createCellChat(
    object   = data_input,
    meta     = seu@meta.data,
    group.by = "singler_pruned"
  )

  cellchat@DB <- CellChatDB.human

  cat("[", format(Sys.time()), "] Running CellChat inference...\n")

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 10)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)

  saveRDS(cellchat, file.path(res_dir, "cellchat_object.rds"))
  cat("[", format(Sys.time()), "] CellChat object saved.\n")

  cellchat
}

# -----------------------------------------------------------------------------
# CellChat visualisations
# -----------------------------------------------------------------------------
plot_cellchat <- function(cellchat, res_dir) {
  cat("[", format(Sys.time()), "] Plotting CellChat results...\n")

  # 1. Interaction count circle plot
  pdf(file.path(res_dir, "cellchat_circle_count.pdf"),
    width = 12,
    height = 12)
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = as.numeric(table(cellchat@idents)),
    weight.scale  = TRUE,
    label.edge    = FALSE,
    title.name    = "Number of interactions"
  )
  dev.off()

  # 2. Interaction strength circle plot
  pdf(file.path(res_dir, "cellchat_circle_strength.pdf"),
    width = 12,
    height = 12)
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = as.numeric(table(cellchat@idents)),
    weight.scale  = TRUE,
    label.edge    = FALSE,
    title.name    = "Interaction strength"
  )
  dev.off()

  # # 3. Outgoing signalling heatmap
  # pdf(file.path(res_dir, "cellchat_pathway_heatmap.pdf"),
  #   width = 12,
  #   height = 12)
  # netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  # dev.off()

  # # 4. Bubble plot
  # p_bubble <- netVisual_bubble(
  #   cellchat,
  #   remove.isolate = FALSE,
  #   angle.x        = 45
  # )
  # ggsave(
  #   file.path(res_dir, "cellchat_bubble.pdf"),p_bubble,
  #   width = 14,
  #   height = 10,
  #   dpi = 600
  # )

  # 5. Top pathways - plot the most significant ones automatically
  pathways <- cellchat@netP$pathways
  cat("[", format(Sys.time()), "] Detected pathways:", paste(pathways, collapse = ", "), "\n")

  # Plot top 3 pathways
  for (pathway in head(pathways, 3)) {
    tryCatch({
      pdf(
        file.path(res_dir, paste0("cellchat_", pathway, ".pdf")),
        width = 12,
        height = 12
      )
      netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
      dev.off()
      cat("[", format(Sys.time()), "] Plotted pathway:", pathway, "\n")
    }, error = function(e) {
      cat("[", format(Sys.time()), "] Could not plot pathway:", pathway, "-", e$message, "\n")
    })
  }
}
