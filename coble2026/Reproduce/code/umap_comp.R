# compare_to_sourcedata.R
# compares reproduced umap embedding against fig1g published source data
# https://arxiv.org/pdf/1802.03426

library(Seurat)
library(vegan)

sourcedata_path <- "csv/origin_fig1g_sourcedata.csv"
cobledata_path <- "csv/coble_fig1g_sourcedata.csv"
out_dir <- "csv"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- load published source data ----

source_data <- read.csv(sourcedata_path, row.names = 1)
colnames(source_data) <- c("umap_1", "umap_2", "CellType")
source_umap <- as.matrix(source_data[, c("umap_1", "umap_2")])

cat("source data cells:", nrow(source_umap), "\n")

# ---- load reproduced object ----

repro_data <- read.csv(cobledata_path, row.names = 1)
colnames(repro_data) <- c("umap_1", "umap_2", "CellType")
repro_umap <- as.matrix(repro_data[, c("umap_1", "umap_2")])

cat("reproduced cells:", nrow(repro_umap), "\n")

# ---- match cells by name ----

common_cells <- intersect(rownames(source_umap), rownames(repro_umap))
cat("common cells:", length(common_cells), "\n")

source_umap <- source_umap[common_cells, ]
repro_umap <- repro_umap[common_cells, ]

# ---- procrustes comparison ----

proc_result <- procrustes(source_umap, repro_umap)
proc_test <- protest(source_umap, repro_umap, permutations = 999)

rmse <- sqrt(proc_result$ss / length(common_cells))
rmse_pct <- 100 * rmse / diff(range(source_umap))

cat("\nsum of squares:", proc_result$ss, "\n")
cat("rmse:", rmse, " (", round(rmse_pct, 2), "% of coord range)\n")
cat("protest correlation:", round(proc_test$t0, 4), " p-value:", proc_test$signif, "\n")

png(file.path(out_dir, "procrustes_sourcedata.png"), width = 800, height = 800)
plot(proc_result, main = "reproduced vs published source data")
dev.off()

# ---- side by side raw scatter ----

png(file.path(out_dir, "sourcedata_side_by_side.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
plot(source_umap, pch = 16, cex = 0.3, col = "steelblue", main = "published source data", xlab = "UMAP_1", ylab = "UMAP_2")
plot(repro_umap, pch = 16, cex = 0.3, col = "firebrick", main = "reproduced", xlab = "UMAP_1", ylab = "UMAP_2")
dev.off()

cat("\nplots saved to", out_dir, "\n")