# compare_to_sourcedata.R
# compares reproduced umap embedding against fig1g published source data
# https://arxiv.org/pdf/1802.03426
#
# change history
# 2026-09-21	claude	write results to a text report instead of console

library(Seurat)
library(vegan)

sourcedata_path <- "sourcedata/origin_fig1g_sourcedata.csv"
cobledata_path <- "sourcedata/coble_fig1g_sourcedata.csv"
out_dir <- "sourcedata"
report_path <- file.path(out_dir, "comparison_report.txt")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- load published source data ----

source_data <- read.csv(sourcedata_path, row.names = 1)
colnames(source_data) <- c("umap_1", "umap_2", "CellType")
source_umap <- as.matrix(source_data[, c("umap_1", "umap_2")])

# ---- load reproduced object ----

repro_data <- read.csv(cobledata_path, row.names = 1)
colnames(repro_data) <- c("umap_1", "umap_2", "CellType")
repro_umap <- as.matrix(repro_data[, c("umap_1", "umap_2")])

# ---- match cells by name ----

common_cells <- intersect(rownames(source_umap), rownames(repro_umap))

source_umap <- source_umap[common_cells, ]
repro_umap <- repro_umap[common_cells, ]

# ---- procrustes comparison ----

proc_result <- procrustes(source_umap, repro_umap)
proc_test <- protest(source_umap, repro_umap, permutations = 999)

rmse <- sqrt(proc_result$ss / length(common_cells))
coord_range <- diff(range(source_umap))
rmse_pct <- 100 * rmse / coord_range

# ---- plots (unchanged) ----

png(file.path(out_dir, "procrustes_sourcedata.png"), width = 800, height = 800)
plot(proc_result, main = "reproduced vs published source data")
dev.off()

png(file.path(out_dir, "sourcedata_side_by_side.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
plot(source_umap, pch = 16, cex = 0.3, col = "steelblue", main = "published source data", xlab = "UMAP_1", ylab = "UMAP_2")
plot(repro_umap, pch = 16, cex = 0.3, col = "firebrick", main = "reproduced", xlab = "UMAP_1", ylab = "UMAP_2")
dev.off()

# ---- write report instead of cat to console ----

report <- c(
"==============================================",
"UMAP reproduction comparison report",
paste("generated:", format(Sys.time())),
"==============================================",
"",
"INPUT FILES",
paste("  published source data:", sourcedata_path),
paste("  reproduced data:", cobledata_path),
"",
"CELL COUNTS",
paste("  source data cells:", nrow(source_data)),
paste("  reproduced cells:", nrow(repro_data)),
paste("  common cells used in comparison:", length(common_cells)),
"",
"------------------------------------------------",
"MAIN STATISTICS  (most important number first)",
"------------------------------------------------",
"",
paste("*** protest correlation:", round(proc_test$t0, 4), "***"),
paste("    p-value:", proc_test$signif, "(999 permutations)"),
"  THIS IS THE MOST IMPORTANT NUMBER. it is a correlation-like",
"  statistic on a 0-1 scale (like Pearson's r) measuring how well",
"  the two shapes align after the best-fit rotation/scale/translation,",
"  tested against random permutations of the data. closer to 1 means",
"  stronger agreement. the p-value tells you whether this alignment",
"  is better than chance (below 0.05 means yes, reliably).",
"",
paste("Procrustes sum of squares:", round(proc_result$ss, 6)),
"  what it is: total squared distance between the two point sets",
"  after the best-fit rotation, scaling and translation. lower is",
"  more similar. not directly comparable across datasets with",
"  different scale or cell counts. mainly useful as the raw number",
"  RMSE below is calculated from.",
"",
paste("RMSE:", round(rmse, 6), "coordinate units"),
paste("  coordinate range (max - min of source umap):", round(coord_range, 6)),
paste("  RMSE as % of coordinate range:", round(rmse, 6), "/", round(coord_range, 6),
      "x 100 =", round(rmse_pct, 2), "%"),
"  what it is: root mean squared error per cell (average positional",
"  offset), expressed as a percentage of the coordinate range so it",
"  can be compared across runs regardless of the plot's absolute scale.",
"",
"------------------------------------------------",
"OUTPUT FILES",
"------------------------------------------------",
paste(" ", file.path(out_dir, "procrustes_sourcedata.png"), "- alignment/residual plot"),
paste(" ", file.path(out_dir, "sourcedata_side_by_side.png"), "- raw umap scatter comparison"),
""
)

writeLines(report, report_path)
cat("report written to", report_path, "\n")