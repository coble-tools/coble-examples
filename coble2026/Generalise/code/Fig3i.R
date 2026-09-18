suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

BACH_DATA_DIR <- "data"
PROC_DATA_DIR <- "data_proc"
ANCHOR        <- "Csn2"
RESULTS_DIR   <- "results"

PLOT_GENES <- c("Lpl", "Acbd7", "Cd36", "Olah")

CELLTYPES_BACH <- c("Lp", "Avd")
EXCLUDE_COND   <- "WTOld"

# GSE164017 stages pooled into "Pregnancy": all non-lactation stages,
# for closer parity with Bach's own pool (Nulliparous+4.5dG+9.5dG+14.5dG).
# Pre-puberty removed: no longer prepared by prepare_GSE164017.R (FVB/NJ-only
# restriction dropped it entirely).
GSE_PREGNANCY_STAGES <- c("Adult virgin", "12.5dG", "18.5dG")
GSE_LACTATION_STAGE  <- "Lactation"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# File input
BACH_RDS_OBJ <- file.path(BACH_DATA_DIR, "BRCA1_SCE.rds")
PAL_RDS_OBJ <- file.path(PROC_DATA_DIR, "gse164017_prepared.rds")

# File output
FIG_BASE  <- file.path(RESULTS_DIR, "Fig2i_bach")
FIG3_BASE <- file.path(RESULTS_DIR, "Supp_bach_timeseries_by_gene")
STATS_CSV <- file.path(RESULTS_DIR, "corr_2i_bach.csv")
STATS3_CSV <- file.path(RESULTS_DIR, "corr_supp_bach_timeseries.csv")

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading ", basename(dest), " ...")
    download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  }
  dest
}

# ---------------------------------------------------------------------------
# LOAD -- Bach's own gestation arm, pooled to a single "Pregnancy" row
# ---------------------------------------------------------------------------

load_bach_pregnancy <- function() {
  for (pkg in c("SingleCellExperiment", "SummarizedExperiment"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " required: BiocManager::install('", pkg, "')")

  rds <- fetch("https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/BRCA1_SCE.rds",BACH_RDS_OBJ)
  sce <- readRDS(rds)
  cd  <- SummarizedExperiment::colData(sce)

  keep <- cd$CellTypesFinal %in% CELLTYPES_BACH & cd$Condition != EXCLUDE_COND &
          cd$Experiment == "Pregnancy"
  message(sprintf("  Bach object: %d cells -> %d Lp/Avd pregnancy cells", ncol(sce), sum(keep)))

  ids <- if (!is.null(colnames(sce))) colnames(sce)[keep]
         else paste0(cd$SampleID[keep], "_", cd$barcode[keep])

  logc <- SummarizedExperiment::assay(sce, "logcounts")[, keep, drop = FALSE]
  colnames(logc) <- make.unique(as.character(ids))

  meta <- data.frame(
    row    = "Pregnancy (Bach)",
    sample = as.character(if ("SampleID" %in% colnames(cd)) cd$SampleID[keep] else "Bach"),
    stringsAsFactors = FALSE
  )
  list(logcounts = logc, meta = meta)
}

# ---------------------------------------------------------------------------
# LOAD -- GSE164017, pooled to "Pregnancy" and "Lactation" rows
# ---------------------------------------------------------------------------

load_gse164017_rows <- function() {
  if (!file.exists(PAL_RDS_OBJ)) stop("missing: ", PAL_RDS_OBJ, " -- run prepare_GSE164017.R first")
  x <- readRDS(PAL_RDS_OBJ)

  row_label <- rep(NA_character_, length(x$condition))
  row_label[x$condition %in% GSE_PREGNANCY_STAGES] <- "Pregnancy (GSE164017)"
  row_label[x$condition == GSE_LACTATION_STAGE]    <- "Lactation (GSE164017)"
  keep <- !is.na(row_label)

  message(sprintf("  GSE164017: %d cells -> %d pregnancy (%s) + lactation",
                  length(x$condition), sum(keep),
                  paste(GSE_PREGNANCY_STAGES, collapse = "/")))

  meta <- data.frame(row = row_label[keep], sample = x$sample[keep], stringsAsFactors = FALSE)
  list(logcounts = x$logcounts[, keep, drop = FALSE], meta = meta)
}

bach <- load_bach_pregnancy()
gse  <- load_gse164017_rows()

# ---------------------------------------------------------------------------
# ASSEMBLE
# ---------------------------------------------------------------------------

genes <- intersect(rownames(bach$logcounts), rownames(gse$logcounts))
message(sprintf("Genes shared by the two datasets: %d", length(genes)))
if (!ANCHOR %in% genes) stop("anchor gene not shared between datasets: ", ANCHOR)

drop_absent <- function(v, nm) {
  keep <- intersect(v, genes)
  if (length(keep) < length(v))
    warning(sprintf("%s: dropping %s (absent from one of the two datasets)",
                    nm, paste(setdiff(v, keep), collapse = ", ")), call. = FALSE)
  keep
}
PLOT_GENES <- drop_absent(PLOT_GENES, "PLOT_GENES")

expr <- cbind(bach$logcounts[genes, , drop = FALSE],
              gse$logcounts[genes, , drop = FALSE])
meta <- rbind(bach$meta, gse$meta)
stopifnot(ncol(expr) == nrow(meta))

ROWS <- c("Pregnancy (Bach)", "Pregnancy (GSE164017)", "Lactation (GSE164017)")



# ---------------------------------------------------------------------------
# CORRELATIONS
# ---------------------------------------------------------------------------

rho_table <- function(group, levels, gene_set) {
  do.call(rbind, lapply(gene_set, function(g) {
    do.call(rbind, lapply(levels, function(rr) {
      i <- which(!is.na(group) & group == rr)
      if (!length(i)) return(NULL)
      x <- as.numeric(expr[ANCHOR, i]); y <- as.numeric(expr[g, i])
      ok <- stats::sd(x) > 0 && stats::sd(y) > 0
      ct <- if (ok) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
            else list(estimate = NA_real_, p.value = NA_real_)
      data.frame(gene = g, row = rr,
                 rho = unname(ct$estimate), p = ct$p.value,
                 pct_detected = 100 * mean(y > 0),
                 pct_detected_anchor = 100 * mean(x > 0),
                 n_cells = length(i), stringsAsFactors = FALSE)
    }))
  }))
}

stats <- rho_table(meta$row, ROWS, PLOT_GENES)

show <- function(s, what, label, gene_set) {
  w <- reshape(s[, c("gene", "row", what)], idvar = "gene",
               timevar = "row", direction = "wide")
  names(w) <- sub(paste0("^", what, "\\."), "", names(w))
  message("\n-- ", label, " --")
  print(w[match(gene_set, w$gene), ], row.names = FALSE, digits = 2)
}
show(stats, "rho", "Bach vs GSE164017 (pregnancy/lactation): Spearman rho with Csn2", PLOT_GENES)

write.csv(stats, STATS_CSV, row.names = FALSE)

# ---------------------------------------------------------------------------
# PLOTTING -- genes as facet columns (labels along the top), rows are the
# three pooled groups. PDF carries the stat annotations (rho/detected/anchor/
# n_cells); PNG is the bare scatter with none of that text.
# ---------------------------------------------------------------------------

GG_34 <- utils::packageVersion("ggplot2") >= "3.4.0"
GG_35 <- utils::packageVersion("ggplot2") >= "3.5.0"
message(sprintf("\nggplot2 %s", as.character(utils::packageVersion("ggplot2"))))
geom_lm_fit <- function(...) {
  a <- list(...)
  if (!GG_34 && "linewidth" %in% names(a)) { a$size <- a$linewidth; a$linewidth <- NULL }
  do.call(geom_smooth, a)
}
el_line <- function(...) {
  a <- list(...)
  if (!GG_34 && "linewidth" %in% names(a)) { a$size <- a$linewidth; a$linewidth <- NULL }
  do.call(element_line, a)
}
geom_lab <- function(...) {
  a <- list(...)
  if (GG_35 && "label.size" %in% names(a)) { a$linewidth <- a$label.size; a$label.size <- NULL }
  do.call(geom_label, a)
}

write_legend <- function(file, title, subtitle, caption) {
  txt <- c(title, subtitle, "", strwrap(caption, width = 100))
  writeLines(txt, paste0(file, "_legend.txt"))
  message("Wrote ", file, "_legend.txt")
}

panel_figure <- function(group, levels, fills, stats, title, subtitle, caption,
                         file, width, height, gene_set, x_breaks = NULL) {
  i <- which(!is.na(group) & group %in% levels)
  long <- do.call(rbind, lapply(gene_set, function(g)
    data.frame(Csn2Expr = as.numeric(expr[ANCHOR, i]),
               Expr = as.numeric(expr[g, i]), gene = g,
               row = group[i], stringsAsFactors = FALSE)))
  long$gene <- factor(long$gene, levels = gene_set)
  long$row  <- factor(long$row,  levels = levels)

  lab <- stats
  lab$gene <- factor(lab$gene, levels = gene_set)
  lab$row  <- factor(lab$row,  levels = levels)
  lab$Correlation <- sprintf("rho == %0.2f", lab$rho)
  lab$Detected    <- sprintf("%0.0f%% of cells", lab$pct_detected)
  lab$Anchor      <- sprintf("%0.0f%% of csn2", lab$pct_detected_anchor)
  lab$n_cells     <- sprintf("%s total cells", format(lab$n_cells, big.mark = ","))
  lab_ok <- lab[!is.na(lab$rho), ]
  lab_na <- lab[is.na(lab$rho), ]
  if (nrow(lab_na)) lab_na$Correlation <- "rho: not estimable"

  # Wash data carries only "row" -- ggplot recycles it across every gene
  # column automatically, same trick the earlier scripts used the other way.
  wash <- data.frame(row = factor(levels, levels = levels))

  # Base plot: everything shared between the labelled (PDF) and unlabelled
  # (PNG) versions -- no rho/detected/anchor/n_cells text yet.
  p_base <- ggplot(long, aes(x = Csn2Expr, y = Expr)) +
    geom_rect(data = wash, aes(fill = row), inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.35) +
    geom_point(size = 0.4, colour = "black", alpha = 0.40) +
    geom_lm_fit(method = "lm", formula = y ~ x, colour = "grey95",
                linewidth = 1.1, se = FALSE) +
    facet_grid(row ~ gene, scales = "fixed") +
    scale_fill_manual(values = fills) +
    (if (!is.null(x_breaks)) scale_x_continuous(breaks = x_breaks) else NULL) +
    labs(x = expression(italic("Csn2") * " expression (log-normalised)"),
         y = "gene expression (log-normalised)") +
    theme_bw(base_size = 11) +
    theme(legend.position   = "none",
          panel.grid.major  = el_line(colour = "grey80", linewidth = 0.1,
                                      linetype = "dashed"),
          panel.grid.minor  = element_blank(),
          strip.background  = element_rect(fill = "white", colour = NA),
          strip.text.x      = element_text(face = "italic", size = 10),
          strip.text.y      = element_text(size = 9.5, angle = -90),
          axis.line  = el_line(linewidth = 0.7),
          axis.ticks = el_line(linewidth = 0.7))

  label_layers <- list(
    geom_lab(data = lab_ok, aes(label = Correlation), parse = TRUE,
             x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15, size = 2.9,
             fontface = "bold", colour = "black", fill = "white", alpha = 0.92,
             label.size = 0, label.padding = unit(0.10, "lines"),
             inherit.aes = FALSE),
    if (nrow(lab_na))
      geom_lab(data = lab_na, aes(label = Correlation), x = -Inf, y = Inf,
               hjust = -0.05, vjust = 1.15, size = 2.5, fontface = "italic",
               colour = "grey30", fill = "white", alpha = 0.92,
               label.size = 0, label.padding = unit(0.10, "lines"),
               inherit.aes = FALSE) else NULL,
    geom_lab(data = lab, aes(label = Detected), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 2.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE),
    geom_lab(data = lab, aes(label = Anchor), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 3.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE),
    geom_lab(data = lab, aes(label = n_cells), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 4.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE)
  )

  p_labelled <- p_base
  for (ly in label_layers) if (!is.null(ly)) p_labelled <- p_labelled + ly

  ggsave(paste0(file, ".pdf"), plot = p_labelled, width = width, height = height,
         limitsize = FALSE)
  ggsave(paste0(file, ".png"), plot = p_base, width = width, height = height, dpi = 300,
         limitsize = FALSE)
  message("Wrote ", file, ".pdf (labelled) and .png (no stat labels)")
  write_legend(file, title, subtitle, caption)
  invisible(p_labelled)
}

# Three distinct washes, one per row -- not a gradient, since these three
# groups aren't a single ordered series (different datasets and conditions).
# Lactation uses a yellow tone, kept clearly distinct from the middle row's
# orange rather than the earlier pink/PuRd.
fills <- setNames(
  c(brewer.pal(5, "GnBu")[5], brewer.pal(6, "Oranges")[4], brewer.pal(9, "YlOrBr")[3]),
  ROWS)

panel_figure(
  group = meta$row, levels = ROWS, fills = fills, stats = stats,
  title = "Fatty-acid coupling to Csn2, gene by gene",
  subtitle = "Bach's own pregnancy average vs GSE164017 pregnancy and lactation",
  caption = paste0(
    "Row 1: Bach et al. 2021 (Nat Commun 12:1502), their Lp/Avd cells restricted to the ",
    "pregnancy arm (Nulliparous/4.5dG/9.5dG/14.5dG pooled), their published logcounts. Rows 2-3: ",
    "GSE164017, 10X Genomics scRNA-seq of mouse mammary gland, cells restricted to ",
    "luminal-progenitor and alveolar compartments (marker-panel argmax, see ",
    "prepare_GSE164017.R); row 2 pools Adult virgin, 12.5dG and 18.5dG, matching the way Bach's ",
    "own row pools its baseline in with the rest; row 3 is Lactation alone. The two datasets come ",
    "from different labs, mice and pipelines and are not jointly normalised. Spearman's rho ",
    "throughout; grey text (PDF only) gives the percentage of cells expressing the target gene ",
    "and Csn2, and the total cells in each panel."),
  file = FIG_BASE, width = 2.3 * length(PLOT_GENES) + 1.5, height = 2.3 * length(ROWS) + 1.8,
  gene_set = PLOT_GENES)

# ---------------------------------------------------------------------------
# SECOND PLOT -- Bach's 4 timepoints, then all 5 GSE164017 stages
# (previously the third of three plots; the Nulliparous-vs-Pre-puberty plot
# that used to sit here has been removed, since Pre-puberty is no longer
# prepared by prepare_GSE164017.R)
# ---------------------------------------------------------------------------

BACH_TIMEPOINTS3 <- c("Nulliparous", "4.5dG", "9.5dG", "14.5dG")
GSE_STAGE_ORDER3 <- c("Adult virgin", "12.5dG", "18.5dG", "Lactation", "Post-involution")

load_bach_three_timepoints <- function() {
  sce <- readRDS(BACH_RDS_OBJ)
  cd  <- SummarizedExperiment::colData(sce)

  keep <- cd$CellTypesFinal %in% CELLTYPES_BACH & cd$Condition != EXCLUDE_COND &
          cd$Experiment == "Pregnancy"
  message(sprintf("  Bach object: %d cells -> %d Lp/Avd cells across all 4 timepoints",
                  ncol(sce), sum(keep)))

  ids <- if (!is.null(colnames(sce))) colnames(sce)[keep]
         else paste0(cd$SampleID[keep], "_", cd$barcode[keep])

  logc <- SummarizedExperiment::assay(sce, "logcounts")[, keep, drop = FALSE]
  colnames(logc) <- make.unique(as.character(ids))

  # WTYoung is Bach's internal label for Nulliparous -- same mapping used
  # elsewhere in this project.
  cond  <- as.character(cd$Condition[keep])
  stage <- ifelse(cond == "WTYoung", "Nulliparous", cond)

  meta <- data.frame(
    row    = paste0(stage, " (Bach)"),
    sample = as.character(if ("SampleID" %in% colnames(cd)) cd$SampleID[keep] else "Bach"),
    stringsAsFactors = FALSE
  )
  list(logcounts = logc, meta = meta)
}

load_gse164017_all_stages <- function() {
  if (!file.exists(PAL_RDS_OBJ)) stop("missing: ", PAL_RDS_OBJ, " -- run prepare_GSE164017.R first")
  x <- readRDS(PAL_RDS_OBJ)
  meta <- data.frame(
    row    = paste0(x$condition, " (GSE164017)"),
    sample = x$sample,
    stringsAsFactors = FALSE
  )
  list(logcounts = x$logcounts, meta = meta)
}

bach3 <- load_bach_three_timepoints()
gse3  <- load_gse164017_all_stages()

genes3 <- intersect(rownames(bach3$logcounts), rownames(gse3$logcounts))
message(sprintf("Genes shared by the two datasets: %d", length(genes3)))
if (!ANCHOR %in% genes3) stop("anchor gene not shared between datasets: ", ANCHOR)

PLOT_GENES3 <- intersect(PLOT_GENES, genes3)
if (length(PLOT_GENES3) < length(PLOT_GENES))
  warning("PLOT_GENES3: dropping genes absent from one of the two datasets", call. = FALSE)

expr <- cbind(bach3$logcounts[genes3, , drop = FALSE],
              gse3$logcounts[genes3, , drop = FALSE])
meta <- rbind(bach3$meta, gse3$meta)
stopifnot(ncol(expr) == nrow(meta))

BACH_ROWS3 <- paste0(BACH_TIMEPOINTS3, " (Bach)")
GSE_ROWS3  <- paste0(GSE_STAGE_ORDER3, " (GSE164017)")
ROWS3      <- c(BACH_ROWS3, GSE_ROWS3)
message("Rows: ", paste(intersect(ROWS3, unique(meta$row)), collapse = ", "))

stats3 <- rho_table(meta$row, ROWS3, PLOT_GENES3)
show(stats3, "rho", "Bach's 4 timepoints + all 5 GSE164017 stages: Spearman rho with Csn2",
     PLOT_GENES3)


write.csv(stats3, STATS3_CSV, row.names = FALSE)

# Blues for Bach's 4 timepoints, orange-to-yellow for GSE164017's 5 stages.
col_bach3 <- colorRampPalette(brewer.pal(5, "GnBu")[2:5])(length(BACH_ROWS3))
col_gse3  <- colorRampPalette(c("#D2691E", "#FFD700"))(length(GSE_ROWS3))
fills3    <- setNames(c(col_bach3, col_gse3), ROWS3)

panel_figure(
  group = meta$row, levels = ROWS3, fills = fills3, stats = stats3,
  title = "Fatty-acid coupling to Csn2, full time series, gene by gene",
  subtitle = "All four Bach timepoints then all five GSE164017 stages",
  caption = paste0(
    "Rows 1-4: Bach et al. 2021 (Nat Commun 12:1502), their Lp/Avd cells, ",
    "Nulliparous/4.5dG/9.5dG/14.5dG, their published logcounts. Rows 5-9: GSE164017, 10X ",
    "Genomics scRNA-seq of mouse mammary gland across life stages, cells restricted to ",
    "luminal-progenitor and alveolar compartments (marker-panel argmax, see ",
    "prepare_GSE164017.R); Post-involution (row 9) is a single mouse. The two datasets come ",
    "from different labs, mice and pipelines and are not jointly normalised, so treat this as two ",
    "independent series shown side by side rather than one merged trajectory. Spearman's rho ",
    "throughout; grey text (PDF only) gives the percentage of cells expressing the target gene ",
    "and Csn2, and the total cells in each panel."),
  file = FIG3_BASE, width = 2.3 * length(PLOT_GENES3) + 1.5, height = 2.3 * length(ROWS3) + 1.8,
  gene_set = PLOT_GENES3)