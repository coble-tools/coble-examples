#!/usr/bin/env Rscript
# =============================================================================
# Driver for Figure 1 of Bach et al. 2021 (Nat Commun 12:1502).
#
# Runs code/Figure1_chunks.R -- the authors' own Figure1.Rmd chunks, verbatim --
# under R-notebook evaluation semantics, inside the COBLE environment built
# from cbl/.coble.cbl.
#
#   SCE_PATH   the authors' distributed SingleCellExperiment
#              (BRCA1_SCE.rds, https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/)
#   OUT_DIR    where the figure files are written
# =============================================================================

SCE_PATH    <- Sys.getenv("SCE_PATH",  "data/BRCA1_SCE.rds")
OUT_DIR     <- Sys.getenv("OUT_DIR",   "results/Figure1/")
CODE_DIR    <- Sys.getenv("CODE_DIR",  "code")
FUNCTIONS_R <- file.path(CODE_DIR, "functions.R")
CHUNKS_R    <- file.path(CODE_DIR, "Figure1_chunks.R")

stopifnot(file.exists(SCE_PATH), file.exists(FUNCTIONS_R), file.exists(CHUNKS_R))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!grepl("/$", OUT_DIR)) OUT_DIR <- paste0(OUT_DIR, "/")

cat("[Figure1] SCE_PATH =", SCE_PATH, "\n")
cat("[Figure1] OUT_DIR  =", OUT_DIR, "\n")

# The notebook auto-printed top-level plot objects; several ggsave() calls in
# the authors' code take no plot= argument and therefore rely on last_plot().
# print.eval=TRUE reproduces that, and the echo device catches the renders.
pdf(file.path(OUT_DIR, "notebook_echo.pdf"), onefile = TRUE, width = 9, height = 7)
t0 <- Sys.time()
source(CHUNKS_R, print.eval = TRUE, echo = FALSE, local = FALSE)
invisible(dev.off())

cat("[Figure1] elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
print(data.frame(file = list.files(OUT_DIR),
                 bytes = file.size(file.path(OUT_DIR, list.files(OUT_DIR)))))
