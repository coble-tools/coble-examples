suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

DATA_DIR    <- "data"
ANCHOR      <- "Csn2"
RESULTS_DIR <- "results"

# The four fatty-acid genes carrying the Fig. 3i claim.
PLOT_GENES <- c("Lpl", "Acbd7", "Cd36", "Olah")
GENES_FIG1 <- PLOT_GENES

# Cell selection, exactly as in their compute_DCARS.R
CELLTYPES_BACH <- c("Lp", "Avd")
EXCLUDE_COND   <- "WTOld"

TIMEPOINTS <- c("Nulliparous", "4.5dG", "9.5dG", "14.5dG")
COLS_FIG1  <- c("Pregnancy (Bach)", "Tumour (Bach)")

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# File input
BACH_RDS_OBJ <- file.path(DATA_DIR, "BRCA1_SCE.rds")

# File output
FIG1_BASE   <- file.path(RESULTS_DIR, "Supp_bach_replication.png")
STATS1_CSV <- file.path(RESULTS_DIR, "corr_supp_bach_replication.csv")

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
# LOAD -- Bach's own processed object
# ---------------------------------------------------------------------------

load_bach <- function() {
  for (pkg in c("SingleCellExperiment", "SummarizedExperiment"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " required: BiocManager::install('", pkg, "')")

  rds <- fetch("https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/BRCA1_SCE.rds", BACH_RDS_OBJ)   # ~2.2 GB
  sce <- readRDS(rds)
  cd  <- SummarizedExperiment::colData(sce)

  need <- c("CellTypesFinal", "Experiment", "Condition")
  if (!all(need %in% colnames(cd)))
    stop("colData columns not found: ", paste(setdiff(need, colnames(cd)), collapse = ", "),
         "\n  available: ", paste(colnames(cd), collapse = ", "))

  keep <- cd$CellTypesFinal %in% CELLTYPES_BACH & cd$Condition != EXCLUDE_COND
  message(sprintf("  Bach object: %d cells -> %d Lp/Avd", ncol(sce), sum(keep)))

  # Their distributed object carries no colnames on the assays; cell ids are in
  # colData. Without this the matrix and the metadata cannot be kept aligned.
  ids <- if (!is.null(colnames(sce))) colnames(sce)[keep]
         else paste0(cd$SampleID[keep], "_", cd$barcode[keep])

  logc <- SummarizedExperiment::assay(sce, "logcounts")[, keep, drop = FALSE]
  colnames(logc) <- make.unique(as.character(ids))

  cond <- as.character(cd$Condition[keep])
  expt <- as.character(cd$Experiment[keep])
  # WTYoung is the nulliparous group of the gestation timecourse.
  timepoint <- ifelse(cond == "WTYoung", "Nulliparous", cond)

  meta <- data.frame(
    fig1 = ifelse(expt == "Pregnancy", "Pregnancy (Bach)", "Tumour (Bach)"),
    sample = as.character(if ("SampleID" %in% colnames(cd)) cd$SampleID[keep] else expt),
    stringsAsFactors = FALSE
  )
  message("  pregnancy timepoints found: ",
          paste(intersect(TIMEPOINTS, unique(timepoint)), collapse = ", "))
  list(logcounts = logc, meta = meta)
}

bach <- load_bach()

# ---------------------------------------------------------------------------
# ASSEMBLE
# ---------------------------------------------------------------------------

expr <- bach$logcounts
meta <- bach$meta
if (!ANCHOR %in% rownames(expr)) stop("anchor gene not in Bach's object: ", ANCHOR)

drop_absent <- function(v, nm) {
  keep <- intersect(v, rownames(expr))
  if (length(keep) < length(v))
    warning(sprintf("%s: dropping %s (absent from Bach's object)",
                    nm, paste(setdiff(v, keep), collapse = ", ")), call. = FALSE)
  keep
}
GENES_FIG1 <- drop_absent(GENES_FIG1, "GENES_FIG1")



# ---------------------------------------------------------------------------
# CORRELATIONS
# ---------------------------------------------------------------------------

rho_table <- function(group, levels, gene_set = PLOT_GENES) {
  do.call(rbind, lapply(gene_set, function(g) {
    do.call(rbind, lapply(levels, function(cc) {
      i <- which(!is.na(group) & group == cc)
      if (!length(i)) return(NULL)
      x <- as.numeric(expr[ANCHOR, i]); y <- as.numeric(expr[g, i])
      # Spearman is undefined if either vector is constant (e.g. Csn2 is
      # detected in 0% of nulliparous cells).
      ok <- stats::sd(x) > 0 && stats::sd(y) > 0
      ct <- if (ok) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
            else list(estimate = NA_real_, p.value = NA_real_)
      data.frame(gene = g, column = cc,
                 rho = unname(ct$estimate), p = ct$p.value,
                 pct_detected = 100 * mean(y > 0),
                 pct_detected_anchor = 100 * mean(x > 0),
                 n_cells = length(i), stringsAsFactors = FALSE)
    }))
  }))
}

stats1 <- rho_table(meta$fig1, COLS_FIG1, GENES_FIG1)

show <- function(s, what, label, gene_set = PLOT_GENES) {
  w <- reshape(s[, c("gene", "column", what)], idvar = "gene",
               timevar = "column", direction = "wide")
  names(w) <- sub(paste0("^", what, "\\."), "", names(w))
  message("\n-- ", label, " --")
  print(w[match(gene_set, w$gene), ], row.names = FALSE, digits = 2)
}
show(stats1, "rho", "FIGURE 1: Spearman rho with Csn2 (must match published values)",
     GENES_FIG1)

write.csv(stats1, STATS1_CSV, row.names = FALSE)

published <- c(Lpl = 0.70, Acbd7 = 0.60, Cd36 = 0.58, Olah = 0.53)
got <- setNames(stats1$rho[stats1$column == "Pregnancy (Bach)"],
                stats1$gene[stats1$column == "Pregnancy (Bach)"])[names(published)]
message("\nPositive control, pregnancy column vs published Fig. 3i:")
message(paste(sprintf("  %-6s got %.2f, published %.2f%s", names(published), got, published,
                      ifelse(abs(got - published) < 0.015, "", "   <-- MISMATCH")),
              collapse = "\n"))

# ---------------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------------
# Bach's idiom: points, an lm fit, a bold parsed "rho == x.xx" label on a
# colour-washed panel. Their washes are GnBu (gestation) and PuRd
# (tumourigenesis).
#
# DEVIATION: their panel also carries a dashed 2-D density contour at
# h = bandwidth.nrd * 6. That setting does not transfer -- the zero-inflation
# and discreteness of log-normalised counts turn it into a mesh that hides the
# points it is meant to summarise. Contours are omitted throughout.

# ggplot2 renamed the line-width aesthetic from `size` to `linewidth` in 3.4.0,
# and geom_label's border argument from `label.size` to `linewidth` in 3.5.0.
# These shims let the script run unchanged on any of those versions.
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
                         file, width, gene_set = PLOT_GENES) {
  i <- which(!is.na(group) & group %in% levels)
  long <- do.call(rbind, lapply(gene_set, function(g)
    data.frame(Csn2Expr = as.numeric(expr[ANCHOR, i]),
               Expr = as.numeric(expr[g, i]), gene = g,
               column = group[i], stringsAsFactors = FALSE)))
  long$gene   <- factor(long$gene, levels = gene_set)
  long$column <- factor(long$column, levels = levels)

  lab <- stats
  lab$gene   <- factor(lab$gene, levels = gene_set)
  lab$column <- factor(lab$column, levels = levels)
  lab$Correlation <- sprintf("rho == %0.2f", lab$rho)
  lab$Detected    <- sprintf("%0.0f%% of cells", lab$pct_detected)
  lab$Anchor    <- sprintf("%0.0f%% of csn2", lab$pct_detected_anchor)
  lab$n_cells <- sprintf("%s total cells", format(lab$n_cells, big.mark = ","))
  lab_ok <- lab[!is.na(lab$rho), ]
  lab_na <- lab[is.na(lab$rho), ]
  if (nrow(lab_na)) lab_na$Correlation <- "rho: not estimable"

  wash <- data.frame(column = factor(levels, levels = levels))

  p <- ggplot(long, aes(x = Csn2Expr, y = Expr)) +
    geom_rect(data = wash, aes(fill = column), inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.35) +
    geom_point(size = 0.4, colour = "black", alpha = 0.40) +
    geom_lm_fit(method = "lm", formula = y ~ x, colour = "grey95",
                linewidth = 1.1, se = FALSE) +
    # Opaque label backgrounds: several panels have dense point clouds in the
    # top-left corner and plain text is illegible against them.
    geom_lab(data = lab_ok, aes(label = Correlation), parse = TRUE,
             x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15, size = 2.9,
             fontface = "bold", colour = "black", fill = "white", alpha = 0.92,
             label.size = 0, label.padding = unit(0.10, "lines"),
             inherit.aes = FALSE) +
    (if (nrow(lab_na))
       geom_lab(data = lab_na, aes(label = Correlation), x = -Inf, y = Inf,
                hjust = -0.05, vjust = 1.15, size = 2.5, fontface = "italic",
                colour = "grey30", fill = "white", alpha = 0.92,
                label.size = 0, label.padding = unit(0.10, "lines"),
                inherit.aes = FALSE) else NULL) +
    geom_lab(data = lab, aes(label = Detected), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 2.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE) +
    geom_lab(data = lab, aes(label = Anchor), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 3.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE) +
    geom_lab(data = lab, aes(label = n_cells), x = -Inf, y = Inf,
             hjust = -0.06, vjust = 4.75, size = 2.5, colour = "grey20",
             fill = "white", alpha = 0.92, label.size = 0,
             label.padding = unit(0.10, "lines"), inherit.aes = FALSE) +
    facet_grid(gene ~ column, scales = "free_y", switch = "y") +
    scale_fill_manual(values = fills) +
    labs(x = expression(italic("Csn2") * " expression (log-normalised)"),
         y = "gene expression (log-normalised)") +
    theme_bw(base_size = 11) +
    theme(legend.position   = "none",
          panel.grid.major  = el_line(colour = "grey80", linewidth = 0.1,
                                      linetype = "dashed"),
          panel.grid.minor  = element_blank(),
          strip.background  = element_rect(fill = "white", colour = NA),
          strip.text.y.left = element_text(face = "italic", size = 10, angle = 0),
          strip.text.x      = element_text(size = 9.5),
          strip.placement   = "outside",
          axis.line  = el_line(linewidth = 0.7),
          axis.ticks = el_line(linewidth = 0.7))
    ggsave(paste0(file, ".pdf"), plot = p, width = width, 
          height = 1.9 * length(gene_set) + 1.8, limitsize = FALSE)
    ggsave(paste0(file, ".png"), plot = p, width = width, 
          height = 1.9 * length(gene_set) + 1.8, dpi = 300, limitsize = FALSE)
    message("Wrote ", file, ".pdf and .png")
    write_legend(file, title, subtitle, caption)
    invisible(p)
}

col_gest <- brewer.pal(n = 5, "GnBu")[5]
col_btum <- brewer.pal(n = 6, "PuRd")[5]

# ---- Figure 1: Bach reproduction --------------------------------------

panel_figure(
  group = meta$fig1, levels = COLS_FIG1,
  fills = setNames(c(col_gest, col_btum), COLS_FIG1),
  stats = stats1,
  title = "Bach et al. Fig. 3i, reproduced on their own data",
  subtitle = "Fatty-acid genes lose their coupling to Csn2 in tumour",
  caption = paste0("Bach et al. 2021 (Nat Commun 12:1502), their processed object ",
    "BRCA1_SCE.rds, restricted as in their compute_DCARS.R to Lp and Avd cells excluding the ",
    "WTOld condition; expression is their published logcounts. Gene list, cell selection and ",
    "panel design from their Figure3_DCARS.Rmd. Values are Spearman's rho; grey text gives the ",
    "percentage of cells in which the gene is detected. Lines are linear least-squares fits, ",
    "shown for fidelity to the original panel; note that a fit can slope where the rank ",
    "correlation is flat."),
  file = FIG1_BASE, width = 6.4, gene_set = GENES_FIG1)