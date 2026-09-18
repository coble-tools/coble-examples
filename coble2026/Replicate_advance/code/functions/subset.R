# =============================================================================
# subset.R
# Cell subsetting and stratified resampling for CellChat
# =============================================================================

# -----------------------------------------------------------------------------
# Subset Seurat object by cell type group
# Returns subsetted Seurat object
# -----------------------------------------------------------------------------
subset_cells <- function(seu, cell_subset, cell_groups) {

  # Build the full group map including "other" and "all"
  all_assigned <- unlist(cell_groups)

  other_types <- unique(seu$singler_pruned[
    !is.na(seu$singler_pruned) & !seu$singler_pruned %in% all_assigned
  ])

  subset_map <- c(
    cell_groups,
    list(
      other = other_types,
      all   = unique(seu$singler_pruned[!is.na(seu$singler_pruned)])
    )
  )

  selected_types <- subset_map[[cell_subset]]

  if (is.null(selected_types) || length(selected_types) == 0) {
    stop("No cell types found for subset: ", cell_subset)
  }

  selected_cells <- colnames(seu)[
    !is.na(seu$singler_pruned) & seu$singler_pruned %in% selected_types
  ]

  if (length(selected_cells) == 0) {
    stop("No cells found after subsetting to: ", cell_subset)
  }

  cat("[", format(Sys.time()), "] Subsetting to", cell_subset, ":\n")
  print(table(seu$singler_pruned[colnames(seu) %in% selected_cells]))
  cat("Total cells:", length(selected_cells), "\n")

  subset(seu, cells = selected_cells)
}

# -----------------------------------------------------------------------------
# Stratified resampling for CellChat
# Takes proportional sample from each cell type with min/max bounds
# min_cells: minimum cells per type (protects rare populations)
# max_cells: maximum cells per type (prevents dominant types swamping)
# prop: proportion of each type to sample
# -----------------------------------------------------------------------------
resample_cells <- function(seu, prop = 0.15, min_cells = 10, max_cells = 500) {
  cat("[", format(Sys.time()), "] Stratified resampling for CellChat...\n")

  set.seed(42)

  cells_keep <- seu@meta.data %>%
    filter(!is.na(singler_pruned)) %>%
    mutate(cell_id = rownames(.)) %>%
    group_by(singler_pruned) %>%
    mutate(
      n_sample = pmax(
        as.integer(min_cells),
        pmin(as.integer(max_cells), ceiling(n() * prop))
      )
    ) %>%
    group_modify(~ slice_sample(.x, n = .x$n_sample[1])) %>%
    pull(cell_id)

  cat("Resampled cells per type:\n")
  print(table(seu$singler_pruned[colnames(seu) %in% cells_keep]))
  cat("Total resampled cells:", length(cells_keep), "\n")

  subset(seu, cells = cells_keep)
}

# -----------------------------------------------------------------------------
# Summarise cell counts per group - useful for QC before subsetting
# -----------------------------------------------------------------------------
summarise_groups <- function(seu, cell_groups) {
  all_assigned <- unlist(cell_groups)

  other_types <- unique(seu$singler_pruned[
    !is.na(seu$singler_pruned) & !seu$singler_pruned %in% all_assigned
  ])

  groups_with_other <- c(cell_groups, list(other = other_types))

  for (group_name in names(groups_with_other)) {
    types  <- groups_with_other[[group_name]]
    cells  <- colnames(seu)[
      !is.na(seu$singler_pruned) & seu$singler_pruned %in% types
    ]
    counts <- sort(
      table(seu$singler_pruned[colnames(seu) %in% cells]),
      decreasing = TRUE
    )

    cat("\n", strrep("=", 40), "\n", sep = "")
    cat(toupper(group_name), "- total cells:", length(cells), "\n")
    cat(strrep("-", 40), "\n", sep = "")
    for (ct in names(counts)) {
      cat(sprintf("  %-45s %d\n", ct, counts[ct]))
    }
  }

  # Sanity check
  total_valid    <- sum(!is.na(seu$singler_pruned))
  total_assigned <- sum(
    !is.na(seu$singler_pruned) &
    seu$singler_pruned %in% c(all_assigned, other_types)
  )

  cat("\n", strrep("=", 40), "\n", sep = "")
  cat(sprintf("Total valid (non-NA) : %d\n", total_valid))
  cat(sprintf("Total in groups      : %d\n", total_assigned))
  cat(sprintf("NAs                  : %d\n", sum(is.na(seu$singler_pruned))))

  if (total_assigned == total_valid) {
    cat("✓ All non-NA cells accounted for\n")
  } else {
    cat("✗ WARNING: some cells missing from groups!\n")
  }
}
