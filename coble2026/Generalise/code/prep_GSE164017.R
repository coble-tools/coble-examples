#!/usr/bin/env Rscript
# =============================================================================
# Prepare GSE164017 as an independent gestation reference
# =============================================================================
#
# WHY
#   Bach et al. 2021 (Nat Commun 12:1502) Fig. 3i contrasts co-expression with
#   Csn2 during gestation against tumourigenesis. Their gestation arm, and the
#   2017 atlas usually used alongside it, both come from the same laboratory,
#   so neither tests whether the gestation baseline holds in anyone else's
#   hands. GSE164017 does: a different lab, different mice, 10x Chromium,
#   with resting adult, mid-pregnancy (12.5dG) and late-pregnancy (18.5dG)
#   samples and two animals at each pregnancy timepoint.
#
#     GSE164017: "10X Genomics scRNA-seq profiling of mammary gland cells from
#     mice at different life stages" (PMID 34187545, 37994351). 18 samples,
#     149,371 cells.
#
# WHAT THIS SCRIPT DOES
#   1. Reads GSE164017_RAW.tar and the shared features file (download these
#      first with the separate bash script, into RAW_DATA_DIR).
#   2. Discovers the per-sample matrices and selects the stages you ask for.
#   3. Applies QC thresholds matched to the other datasets in this comparison.
#   4. Assigns each cell to a compartment, and keeps the luminal-progenitor /
#      alveolar cells -- the population equivalent to Bach's Lp + Avd union,
#      which is what their compute_DCARS.R restricts to.
#   5. Writes a prepared object plus the evidence you need to judge whether the
#      compartment call is trustworthy.
#
# HOW COMPARTMENTS ARE CALLED, AND WHY NOT BY CLUSTERING
#   Cells are assigned by marker-panel scores rather than by clustering and
#   annotating clusters. The reason is specific to this dataset: the stages
#   differ enormously in transcriptome (a lactating cell is dominated by milk
#   transcripts), so a merged clustering separates largely BY STAGE and would
#   need batch integration before the clusters meant anything -- and the choice
#   of integration method would then sit underneath every downstream number.
#   Scoring each cell independently against fixed marker panels avoids that:
#   it is deterministic, batch-free, and inspectable. The cost is that it is
#   coarser than a good clustering, so the script writes out the per-compartment
#   marker means and a figure for you to check the call rather than trust it.
#   If you would rather use clustering, the object written at step 5 keeps all
#   epithelial cells with their scores, so you can re-label from there.
#
# USAGE
#   Rscript prepare_GSE164017.R
#
# OUTPUT
#   PROC_DATA_DIR (the prepared data, for downstream use):
#     gse164017_prepared.rds          list(logcounts, condition, sample,
#                                          compartment) for Lp/alveolar cells
#     gse164017_epithelial.rds        all epithelial cells with compartment
#                                          calls and scores, if you want to
#                                          re-label
#   QC_DATA_DIR (diagnostics, for you to check before trusting the above):
#     gse164017_qc_summary.csv        cells before and after QC, per sample
#     gse164017_compartment_calls.csv compartment x stage cell counts
#     gse164017_marker_means.csv      mean expression of every marker panel in
#                                          every compartment -- the evidence
#     gse164017_anchor_detection.csv  Csn2 and fatty-acid gene detection rates
#                                          in the kept cells, per stage
#     gse164017_anchor_variance.csv   Csn2 detection and spread, all vs kept
#     gse164017_compartment_check_a.png / _b.png  verification figures
#
# DEPENDENCIES
#   Matrix, ggplot2  (both CRAN). No Seurat or Bioconductor needed.
#   Disk: ~2.5 GB during extraction. Raw tarball downloaded separately.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

RAW_DATA_DIR  <- "data_raw"
PROC_DATA_DIR <- "data_proc"
QC_DATA_DIR   <- "data_qc"

# Which stages to prepare. Keys are labels used downstream; values are regular
# expressions matched against the sample name in the GEO filenames.
#
# RESTRICTED TO FVB/NJ MICE ONLY. Every sample below has been individually
#   confirmed FVB/NJ via its own GSM record.
#   "Adult virgin" narrowed from "^Adult" to "^Adult-FVB", dropping Adult-BL6
#   (C57BL/6) and Adult-SW (Swiss) -- keeps Adult-FVB, Adult-FVB-P7,
#   Adult-FVB-D12, all confirmed FVB/NJ.
#   "Pre-puberty" (Pre-BL6, confirmed C57BL/6) and "Early postnatal"
#   (Pre-D5-BL6, C57BL/6 by name) removed entirely.
#   12.5dG (Preg-D12, Preg-D12-G), 18.5dG (Preg-D18, Preg-D18-B4), Lactation
#   (Lac-D10, Lac-D10-D6) and Post-involution (PI-W3) are all confirmed
#   FVB/NJ.
STAGE_PATTERNS <- c(
  "Adult virgin"    = "^Adult-FVB",
  "12.5dG"          = "^Preg-?D12",
  "18.5dG"          = "^Preg-?D18",
  "Lactation"       = "^Lac",
  "Post-involution" = "^PI-"
)

# QC, matched to the thresholds used for the tumour data in this comparison
MIN_GENES <- 500
MIN_UMIS  <- 1000
MAX_MITO  <- 20      # percent

# Marker panels. First block separates epithelium from everything else; second
# block splits the epithelium into its three lineages. Genes absent from the
# annotation are dropped with a warning rather than failing.
LINEAGE_PANELS <- list(
  immune      = c("Ptprc", "Cd52", "Laptm5", "Coro1a", "Lyz2", "Cd3e", "Cd79a"),
  fibroblast  = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfrb"),
  endothelial = c("Pecam1", "Cdh5", "Eng", "Egfl7"),
  epithelial  = c("Epcam", "Krt8", "Krt18", "Cdh1", "Krt19")
)
COMPARTMENT_PANELS <- list(
  basal           = c("Krt5", "Krt14", "Krt17", "Acta2", "Myl9", "Trp63"),
  hormone_sensing = c("Esr1", "Pgr", "Prlr", "Foxa1", "Areg", "Cited1"),
  luminal_prog    = c("Elf5", "Kit", "Aldh1a3", "Cd14", "Lypd3", "Kcnk6"),
  alveolar        = c("Csn1s1", "Csn1s2b", "Wap", "Lalba", "Glycam1")
)
# NOTE: Csn2 is deliberately ABSENT from the alveolar panel above. It is the
# anchor gene of the downstream correlation, and selecting cells on the anchor
# then correlating against it within the selected cells conditions on the
# variable under test -- the first run of this script included Csn2 here and
# produced a kept set in which Csn2 was detected in 100% of gestation cells,
# which is selection, not biology. The remaining markers are still part of the
# secretory programme and so correlated with Csn2, which is a residual and
# unavoidable limitation: Bach's own Avd cluster has the same property, being
# defined by clustering that milk genes dominate. Keeping the luminal-progenitor
# compartment alongside the alveolar one is what preserves the full range of
# Csn2 in the kept cells.
# Compartments kept for the figure: the equivalent of Bach's Lp + Avd union.
KEEP_COMPARTMENTS <- c("luminal_prog", "alveolar")

# Genes the downstream figure needs to be present
ANCHOR   <- "Csn2"
FA_GENES <- c("Lpl", "Acbd7", "Cd36", "Olah")

set.seed(42)
dir.create(PROC_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QC_DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. LOCATE RAW DATA (download handled separately, see download script)
# ---------------------------------------------------------------------------

fetch <- function(path) {
  if (!file.exists(path)) stop("missing: ", path, " -- run the download script first")
  path
}

tarball  <- fetch(file.path(RAW_DATA_DIR, "GSE164017_RAW.tar"))
features <- fetch(file.path(RAW_DATA_DIR, "GSE164017_features.tsv.gz"))

raw_dir <- file.path(RAW_DATA_DIR, "raw")
if (!dir.exists(raw_dir) || !length(list.files(raw_dir, pattern = "matrix\\.mtx"))) {
  dir.create(raw_dir, showWarnings = FALSE)
  message("Extracting archive ...")
  untar(tarball, exdir = raw_dir)
}

# ---------------------------------------------------------------------------
# 2. DISCOVER SAMPLES
# ---------------------------------------------------------------------------
# Files are named GSMxxxxxxx_<sample>-matrix.mtx.gz / -barcodes.tsv.gz, and the
# feature annotation is shared across the whole series.

mtx <- list.files(raw_dir, pattern = "-matrix\\.mtx\\.gz$")
stems <- sub("-matrix\\.mtx\\.gz$", "", mtx)
names(stems) <- sub("^GSM\\d+_", "", stems)
message("\nSamples discovered in the archive (", length(stems), "):")
message(paste(" ", sort(names(stems)), collapse = "\n"))

assign_stage <- function(nm) {
  hit <- names(STAGE_PATTERNS)[vapply(STAGE_PATTERNS, grepl, logical(1), x = nm)]
  if (length(hit) == 1) hit else NA_character_
}
stage <- vapply(names(stems), assign_stage, character(1))
selected <- stems[!is.na(stage)]
sel_stage <- stage[!is.na(stage)]
if (!length(selected))
  stop("No samples matched STAGE_PATTERNS. Edit it using the listing above.")
message("\nSelected for preparation:")
message(paste(sprintf("  %-22s -> %s", names(selected), sel_stage), collapse = "\n"))

# ---------------------------------------------------------------------------
# 3. LOAD AND QC
# ---------------------------------------------------------------------------

gene_tbl <- read.delim(gzfile(features), header = FALSE, stringsAsFactors = FALSE)
symbols  <- gene_tbl[[if (ncol(gene_tbl) >= 2) 2 else 1]]

qc_rows <- list()
message("Reading and filtering samples:")
mats <- lapply(seq_along(selected), function(k) {
  stem <- selected[[k]]; nm <- names(selected)[k]
  m <- as(Matrix::readMM(gzfile(file.path(raw_dir, paste0(stem, "-matrix.mtx.gz")))),
          "CsparseMatrix")
  bcs <- readLines(gzfile(file.path(raw_dir, paste0(stem, "-barcodes.tsv.gz"))))
  if (nrow(m) != length(symbols))
    stop(sprintf("%s: matrix has %d rows but the shared features file has %d",
                 nm, nrow(m), length(symbols)))
  rownames(m) <- symbols
  colnames(m) <- paste0(nm, "_", bcs)
  m <- m[!duplicated(rownames(m)), , drop = FALSE]

  total   <- Matrix::colSums(m)
  n_genes <- Matrix::colSums(m > 0)
  mito    <- Matrix::colSums(m[grepl("^mt-", rownames(m), ignore.case = TRUE), ,
                               drop = FALSE]) / pmax(total, 1) * 100
  keep <- n_genes >= MIN_GENES & total >= MIN_UMIS & mito < MAX_MITO
  qc_rows[[nm]] <<- data.frame(sample = nm, stage = sel_stage[k],
                               cells_before = ncol(m), cells_after = sum(keep),
                               median_umis = median(total[keep]),
                               median_genes = median(n_genes[keep]))
  message(sprintf("  %-22s %6d -> %6d cells", nm, ncol(m), sum(keep)))
  m[, keep, drop = FALSE]
})

genes  <- Reduce(intersect, lapply(mats, rownames))
counts <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
rm(mats); invisible(gc(verbose = FALSE))

qc <- do.call(rbind, qc_rows)
meta <- data.frame(
  cell   = colnames(counts),
  sample = rep(qc$sample, qc$cells_after),
  stage  = rep(qc$stage,  qc$cells_after),
  stringsAsFactors = FALSE
)
stopifnot(nrow(meta) == ncol(counts), !any(is.na(meta$stage)))
write.csv(qc, file.path(QC_DATA_DIR, "gse164017_qc_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. NORMALISE AND SCORE
# ---------------------------------------------------------------------------
# log1p(counts / total per cell * 1e4) -- identical to Seurat's LogNormalize,
# and to the normalisation used for the tumour data in this comparison.

sf <- Matrix::colSums(counts); sf[sf == 0] <- 1
logc <- counts %*% Matrix::Diagonal(x = 1e4 / sf)
dimnames(logc) <- dimnames(counts)
logc <- log1p(logc)

# Panel score = mean z-scored expression of the panel's genes. Z-scoring across
# cells puts panels of different absolute expression on a comparable footing,
# so the argmax across panels is meaningful.
row_z <- function(g) {
  v <- as.numeric(logc[g, ])
  s <- stats::sd(v)
  if (s == 0) rep(0, length(v)) else (v - mean(v)) / s
}
panel_score <- function(panel, label) {
  present <- intersect(panel, rownames(logc))
  if (length(present) < length(panel))
    warning(sprintf("%s: %d/%d marker genes found (missing: %s)", label,
                    length(present), length(panel),
                    paste(setdiff(panel, present), collapse = ", ")),
            call. = FALSE)
  if (!length(present)) stop("No marker genes found for panel ", label)
  Reduce(`+`, lapply(present, row_z)) / length(present)
}

# Is there any non-epithelial population to gate out? Several mammary series
# (this one included) are FACS-sorted Lin- epithelium, in which case the
# immune/endothelial/fibroblast panels are uniformly zero. Taking an argmax
# across panels that are all zero assigns cells arbitrarily: the first run of
# this script called 46,934 of 84,000 cells "non-epithelial" in a sorted
# dataset, and inspection showed they were myoepithelial cells (Acta2 3.8,
# Krt14 2.1) losing the argmax to an all-zero fibroblast panel. So gate only
# when there is something to gate.
nonepi_markers <- intersect(c("Ptprc", "Cd52", "Pecam1", "Col1a1", "Dcn"),
                            rownames(logc))
nonepi_pct <- 100 * Matrix::rowMeans(logc[nonepi_markers, , drop = FALSE] > 0)
message("\nNon-epithelial marker detection across all cells (%):")
print(round(nonepi_pct, 2))
# Judge on immune and endothelial markers only. Col1a1 is expressed by some
# myoepithelial cells, so including it put this sorted dataset over a 2%
# threshold on the previous run and left the broken gating switched on.
# Ptprc and Pecam1 are unambiguous and are abundant in any unsorted gland.
PRESORTED <- max(nonepi_pct[intersect(c("Ptprc", "Pecam1"),
                                      names(nonepi_pct))]) < 5
message(if (PRESORTED)
  "  -> treating as pre-sorted epithelium; lineage gating SKIPPED"
  else "  -> non-epithelial cells present; lineage gating applied")

compart <- vapply(names(COMPARTMENT_PANELS),
                  function(p) panel_score(COMPARTMENT_PANELS[[p]], p),
                  numeric(ncol(logc)))
meta$compartment <- colnames(compart)[apply(compart, 1, which.max)]

if (PRESORTED) {
  meta$lineage <- "epithelial"
} else {
  lineage <- vapply(names(LINEAGE_PANELS),
                    function(p) panel_score(LINEAGE_PANELS[[p]], p),
                    numeric(ncol(logc)))
  meta$lineage <- colnames(lineage)[apply(lineage, 1, which.max)]
  meta$compartment[meta$lineage != "epithelial"] <- "non-epithelial"
}
meta <- cbind(meta, as.data.frame(compart))

message("\nLineage assignment (all cells):")
print(table(meta$lineage, meta$stage))
message("\nCompartment assignment:")
print(table(meta$compartment, meta$stage))
write.csv(as.data.frame.matrix(table(meta$compartment, meta$stage)),
          file.path(QC_DATA_DIR, "gse164017_compartment_calls.csv"))

# ---- the evidence: marker means per compartment ---------------------------

all_markers <- unique(c(ANCHOR, unlist(c(LINEAGE_PANELS, COMPARTMENT_PANELS))))
all_markers <- intersect(all_markers, rownames(logc))
marker_means <- do.call(rbind, lapply(sort(unique(meta$compartment)), function(cc) {
  i <- meta$compartment == cc
  data.frame(compartment = cc, n_cells = sum(i), gene = all_markers,
             mean_logexpr = Matrix::rowMeans(logc[all_markers, i, drop = FALSE]),
             pct_detected = 100 * Matrix::rowMeans(logc[all_markers, i, drop = FALSE] > 0),
             row.names = NULL, stringsAsFactors = FALSE)
}))
write.csv(marker_means, file.path(QC_DATA_DIR, "gse164017_marker_means.csv"),
          row.names = FALSE)

message("\nMean log-expression of key markers by compartment:")
key <- intersect(c("Epcam", "Ptprc", "Col1a1", "Krt14", "Acta2", "Esr1", "Prlr",
                   "Elf5", "Kit", "Aldh1a3", "Csn2", "Wap"), all_markers)
wide <- reshape(marker_means[marker_means$gene %in% key,
                             c("compartment", "gene", "mean_logexpr")],
                idvar = "gene", timevar = "compartment", direction = "wide")
names(wide) <- sub("^mean_logexpr\\.", "", names(wide))
print(wide[match(key, wide$gene), ], row.names = FALSE, digits = 2)

# ---------------------------------------------------------------------------
# 5. SUBSET AND WRITE
# ---------------------------------------------------------------------------

keep_i <- meta$compartment %in% KEEP_COMPARTMENTS
if (!sum(keep_i)) stop("No cells assigned to ", paste(KEEP_COMPARTMENTS, collapse = "/"))
message(sprintf("\nKeeping %d cells in %s (of %d total)",
                sum(keep_i), paste(KEEP_COMPARTMENTS, collapse = " + "), nrow(meta)))

need <- c(ANCHOR, FA_GENES)
absent <- setdiff(need, rownames(logc))
if (length(absent)) warning("genes needed downstream are absent: ",
                            paste(absent, collapse = ", "), call. = FALSE)

# Is the anchor gene actually usable in each stage? A correlation cannot be
# estimated where Csn2 is detected in ~0% of cells, which is the expected
# situation in the virgin gland.
det <- do.call(rbind, lapply(sort(unique(meta$stage)), function(st) {
  i <- keep_i & meta$stage == st
  data.frame(stage = st, n_cells = sum(i),
             gene = intersect(need, rownames(logc)),
             pct_detected = 100 * Matrix::rowMeans(
               logc[intersect(need, rownames(logc)), i, drop = FALSE] > 0),
             row.names = NULL)
}))
write.csv(det, file.path(QC_DATA_DIR, "gse164017_anchor_detection.csv"),
          row.names = FALSE)
message("\nDetection in the kept cells (%, by stage):")
dw <- reshape(det[, c("stage", "gene", "pct_detected")], idvar = "gene",
              timevar = "stage", direction = "wide")
names(dw) <- sub("^pct_detected\\.", "", names(dw))
print(dw[match(intersect(need, dw$gene), dw$gene), ], row.names = FALSE, digits = 2)

anchor_stats <- do.call(rbind, lapply(sort(unique(meta$stage)), function(st) {
  a <- as.numeric(logc[ANCHOR, meta$stage == st])
  k <- as.numeric(logc[ANCHOR, keep_i & meta$stage == st])
  data.frame(stage = st,
             n_all = length(a), pct_det_all = 100 * mean(a > 0), sd_all = sd(a),
             n_kept = length(k), pct_det_kept = 100 * mean(k > 0), sd_kept = sd(k))
}))
write.csv(anchor_stats, file.path(QC_DATA_DIR, "gse164017_anchor_variance.csv"),
          row.names = FALSE)
message("\n", ANCHOR, ": detection and spread, all cells vs kept cells")
message("(high detection is not a problem for a rank correlation; zero SPREAD is)")
print(anchor_stats, row.names = FALSE, digits = 3)

saveRDS(list(
  logcounts   = logc[, keep_i, drop = FALSE],
  condition   = meta$stage[keep_i],
  sample      = meta$sample[keep_i],
  compartment = meta$compartment[keep_i],
  provenance  = list(accession = "GSE164017", prepared = Sys.time(),
                     qc = c(min_genes = MIN_GENES, min_umis = MIN_UMIS,
                            max_mito_pct = MAX_MITO),
                     compartments_kept = KEEP_COMPARTMENTS,
                     normalisation = "log1p(counts/colSums * 1e4)",
                     compartment_method = "marker-panel argmax on z-scored expression")
), file.path(PROC_DATA_DIR, "gse164017_prepared.rds"))

epi_i <- meta$lineage == "epithelial"
saveRDS(list(logcounts = logc[, epi_i, drop = FALSE], meta = meta[epi_i, ]),
        file.path(PROC_DATA_DIR, "gse164017_epithelial.rds"))
message("\nWrote gse164017_prepared.rds and gse164017_epithelial.rds to ", PROC_DATA_DIR)

# ---------------------------------------------------------------------------
# 6. VERIFICATION FIGURE
# ---------------------------------------------------------------------------
# Panel a: marker means by compartment -- does each compartment express the
#          markers it is named for, and not the others?
# Panel b: Csn2 detection by stage in the kept cells -- is the anchor usable?

pa_dat <- marker_means[marker_means$gene %in% key, ]
pa_dat$gene <- factor(pa_dat$gene, levels = key)
pa_dat$compartment <- factor(pa_dat$compartment)
pa <- ggplot(pa_dat, aes(compartment, gene)) +
  geom_point(aes(size = pct_detected, colour = mean_logexpr)) +
  scale_colour_gradient(low = "grey90", high = "#B02418", name = "mean\nlog-expr") +
  scale_size_continuous(range = c(0.5, 6), name = "% detected") +
  labs(title = "a  Compartment calls against their marker panels",
       subtitle = "each compartment should be brightest for the genes it is named for",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        panel.grid.minor = element_blank())

pb <- ggplot(anchor_stats, aes(stage, sd_kept)) +
  geom_col(fill = "#1F5FA8", width = 0.6) +
  geom_text(aes(label = sprintf("SD %.2f\n%.0f%% detected\nn=%s", sd_kept,
                                pct_det_kept, format(n_kept, big.mark = ","))),
            vjust = -0.25, size = 2.6, lineheight = 0.95) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(title = sprintf("b  Spread of %s in the kept cells", ANCHOR),
       subtitle = paste("a rank correlation needs spread, not absence:",
                        "\nnear-100% detection is fine, zero SD is not"),
       x = NULL, y = sprintf("SD of %s (log-normalised)", ANCHOR)) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

ggsave(file.path(QC_DATA_DIR, "gse164017_compartment_check_a.png"), pa,
       width = 5.2, height = 4.0, dpi = 300)
ggsave(file.path(QC_DATA_DIR, "gse164017_compartment_check_b.png"), pb,
       width = 5.4, height = 3.6, dpi = 300)
message("Wrote gse164017_compartment_check_a.png and _b.png to ", QC_DATA_DIR)

message("\n--------------------------------------------------------------------")
message("BEFORE USING THIS DOWNSTREAM, CHECK:")
message(" 1. Panel a: luminal_prog is Elf5/Kit/Aldh1a3-high and Krt14/Esr1-low;")
message("    alveolar is Csn2/Wap-high. If the compartments are smeared, the")
message("    marker-score call has failed and clustering is needed instead.")
message(" 2. Panel b: the anchor must have non-zero SPREAD in each stage it is")
message("    used for. High detection is not a defect; zero SD is.")
message(" 3. gse164017_compartment_calls.csv: the alveolar compartment should")
message("    appear at 12.5/18.5dG and be near-absent in the virgin gland.")
message("--------------------------------------------------------------------")