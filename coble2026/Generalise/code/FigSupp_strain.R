suppressPackageStartupMessages(library(Matrix))

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

BACH_DATA_DIR <- "data"
PROC_DATA_DIR <- "data_proc"
ANCHOR        <- "Csn2"
PLOT_GENES    <- c("Lpl", "Acbd7", "Cd36", "Olah")
GENES         <- c(ANCHOR, PLOT_GENES)
RESULTS_DIR   <- "results"

CELLTYPES_BACH <- c("Lp", "Avd")

# Strain from GEO sample titles, except Adult-FVB-D12 confirmed directly via
# its own GSM record. Pre-BL6's strain is name-inferred only.
GSE_SAMPLES <- list(
  "Adult-FVB"     = list(strain = "FVB/NJ",  stage = "Adult virgin"),
  "Adult-BL6"     = list(strain = "C57BL/6", stage = "Adult virgin"),
  "Adult-SW"      = list(strain = "Swiss",   stage = "Adult virgin"),
  "Adult-FVB-P7"  = list(strain = "FVB/NJ",  stage = "Adult virgin"),
  "Adult-FVB-D12" = list(strain = "FVB/NJ",  stage = "Adult virgin"),
  "Pre-BL6"       = list(strain = "C57BL/6", stage = "Pre-puberty")
)

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# File input
BACH_RDS_OBJ <- file.path(BACH_DATA_DIR, "BRCA1_SCE.rds")
PAL_RDS_OBJ <- file.path(PROC_DATA_DIR, "gse164017_prepared.rds")

# File output
OUT_CSV <- file.path(RESULTS_DIR, "virgin_gene_detection_by_strain.csv")

# ---------------------------------------------------------------------------
# GSE164017: per-sample, per-gene detection and mean expression
# ---------------------------------------------------------------------------

x <- readRDS(PAL_RDS_OBJ)

gse_rows <- do.call(rbind, lapply(names(GSE_SAMPLES), function(s) {
  info <- GSE_SAMPLES[[s]]
  i <- x$condition == info$stage & x$sample == s
  do.call(rbind, lapply(GENES, function(g) {
    v <- as.numeric(x$logcounts[g, i])
    data.frame(sample = s, stage = info$stage, strain = info$strain,
               source = "GSE164017", gene = g, n_cells = sum(i),
               pct_detected = 100 * mean(v > 0), mean_logexpr = mean(v))
  }))
}))

# ---------------------------------------------------------------------------
# Bach: Nulliparous (WTYoung), confirmed C57BL/6 throughout (ArrayExpress
# E-MTAB-10043 SDRF)
# ---------------------------------------------------------------------------

for (pkg in c("SingleCellExperiment", "SummarizedExperiment"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, " required: BiocManager::install('", pkg, "')")

sce <- readRDS(BACH_RDS_OBJ)
cd  <- SummarizedExperiment::colData(sce)
keep <- cd$CellTypesFinal %in% CELLTYPES_BACH & cd$Condition == "WTYoung" &
        cd$Experiment == "Pregnancy"
logc_bach <- SummarizedExperiment::assay(sce, "logcounts")[, keep, drop = FALSE]

bach_rows <- do.call(rbind, lapply(GENES, function(g) {
  v <- as.numeric(logc_bach[g, ])
  data.frame(sample = "Bach Nulliparous", stage = "Nulliparous", strain = "C57BL/6",
             source = "Bach", gene = g, n_cells = sum(keep),
             pct_detected = 100 * mean(v > 0), mean_logexpr = mean(v))
}))

# ---------------------------------------------------------------------------
# WRITE
# ---------------------------------------------------------------------------

d <- rbind(gse_rows, bach_rows)
d$sample <- factor(d$sample, levels = c(names(GSE_SAMPLES), "Bach Nulliparous"))
d$gene   <- factor(d$gene, levels = GENES)
d <- d[order(d$gene, d$sample), ]


write.csv(d, OUT_CSV, row.names = FALSE)
message("Wrote ", OUT_CSV)
print(d, row.names = FALSE, digits = 3)