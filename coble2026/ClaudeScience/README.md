# Claude Science Reproduction of Bach et al

Given this prompt:
```
I want you to build an environment so that I can replicate a computational result from a publication, using the COBLE tool. I want you to use the environment internally to create the results and then give me the coble file and the script so I can verify them. Here is the doumentation for COBLE and the repo: 

https://github.com/coble-tools/coble
https://coble-tools.github.io/coble/

And here is an example:
https://github.com/coble-tools/coble-examples/blob/main/coble2026/Reproduce/cbl/.coble.cblI would like you to make an environment to replicate this paper's code.
https://www.nature.com/articles/s41467-021-21783-3
```
The covble scrit and the code scripts were created by Claude Science and run by me on my local machine - fully runnable environment from the coble definition and the Claude Science code scripts ran without error.

Everything below is written by Claude Science and all other files in this folder are created by Claude Science.


# Reproducing Figure 1 of Bach et al. 2021 with COBLE

*Bach, K., Pensa, S., Zarocsinceva, M. et al. Time-resolved single-cell analysis of
Brca1 associated mammary tumourigenesis reveals aberrant differentiation of luminal
progenitors. **Nat Commun** 12, 1502 (2021).*
<https://doi.org/10.1038/s41467-021-21783-3>

Publication: <https://www.nature.com/articles/s41467-021-21783-3>
Authors' code: <https://github.com/MarioniLab/Tumorigenesis2018>
Authors' processed data: <https://content.cruk.cam.ac.uk/jmlab/BRCA1Tumourigenesis/>

`same inputs, different team` — the authors' own analysis code and their own
distributed object, run in a COBLE environment rebuilt from the software
generation they used.

**Target figure:** Figure 1 — the integrated single-cell atlas of the Brca1/p53
tumourigenesis time course and the wild-type gestation time course.
Produced by `src/DownstreamAnalysis/Figure1.Rmd`: the cell-type UMAP of the
integrated atlas, the same UMAP split by tumour-time bin and by condition, and
the stacked barplot of compartment composition per condition.

---

## Layout

```
bach2021/
  cbl/
    .coble.cbl          the recipe -- the intent
    .coble.sh           the bash recipe COBLE generated from it
    .coble_export.cbl   the frozen environment COBLE captured after building
    .coble_summary.txt  per-command timing log of the build
    validate.sh         version assertions, installed as `coble validate`
  code/
    Figure1.R           driver: runs the chunks under notebook semantics
    Figure1_chunks.R    the authors' Figure1.Rmd chunks, verbatim
    functions.R         vendored unchanged from MarioniLab/Tumorigenesis2018
  data/                 BRCA1_SCE.rds, downloaded by run_figure1.sh
  results/Figure1/      output
  run_figure1.sh        download + run, in conda, docker or singularity
```

## Build and run

```bash
conda install rachelsa::coble        # or clone coble-tools/coble and add code/ to PATH

coble build --recipe cbl/.coble.cbl --env bach-rep --rebuild --validate cbl/validate.sh
# containers instead of / as well as conda:
# coble build --recipe cbl/.coble.cbl --env bach-rep --containers docker,singularity --validate cbl/validate.sh

conda activate bach-rep
coble me                             # who built it, when, from which recipe
coble validate                       # assert the versions the figure depends on

bash run_figure1.sh                  # downloads BRCA1_SCE.rds (2.0 GB) then runs
# bash run_figure1.sh singularity
```

To rebuild the environment exactly rather than by intent:

```bash
coble build --recipe cbl/.coble_export.cbl --env bach-rep-mirror
```

---

## The environment, and why these versions

The authors' analysis container (`Container/Tumorigenesis.def`) builds R 3.6.3
from source on Ubuntu and then installs the packages in
`Container/ListOfPackages.txt` with `BiocManager::install()`. Nothing is pinned:
the versions you get depend on the day you build. For R 3.6.3 that resolves to
**Bioconductor 3.10**, so the recipe states that release pair explicitly rather
than leaving it to the calendar:

| | recipe | why |
|---|---|---|
| R | 3.6.3 | the version built in `Tumorigenesis.def` |
| SingleCellExperiment | 1.8.0 | Bioconductor 3.10 |
| scater | 1.14.0 | Bioconductor 3.10 |
| scran | 1.14.1 | Bioconductor 3.10 |
| edgeR | 3.28.0 | Bioconductor 3.10 |
| ggplot2, cowplot, dplyr, plyr, ggrepel, viridis, RColorBrewer, ggthemes, Matrix, svglite | CRAN as of March 2021 | `Container/README.md` states the downstream plotting scripts were run with a local R, not the container; March 2021 is publication |
| ggrastr | 0.2.3 | release current at publication; no conda build exists for R 3.6, so it is fetched from the CRAN archive. It has no compiled code, and its Cairo, ggbeeswarm, png and ragg dependencies come from conda, so nothing is compiled during the build |

The paper's data-availability statement says R 3.4.1, which does not match the
container in the authors' own repository (R 3.6.3). The container is the
better evidence and is what the recipe follows.

### Two things the recipe has to work around

Both are 2019-era packaging that has since rotted; both are fixed in the recipe
rather than by hand, so the build is repeatable.

1. **`libopenblas`.** The bioconda r36 binaries are linked against
   `libopenblasp-r0.3.7.dylib` by soname. Left free, the solver takes the
   current libopenblas, and `edgeR.dylib` — and therefore `scran` — fails to
   `dyn.load`. The recipe pins `libopenblas=0.3.7` before anything Bioconductor
   lands.
2. **`bioconductor-genomeinfodbdata`.** Its conda post-link script fetches its
   data tarball with `curl` and no `-L`. bioconductor.org has since moved to 302
   redirects, so the fetch quietly stores the redirect page, the md5 check
   fails, and the whole conda transaction aborts. The recipe's `bash:` section
   writes a `.curlrc` containing `location` into the environment and points
   `CURL_HOME` at it.

### Deliberately not included

`gganimate` / `gifski`. `Figure1.Rmd` ends with a block that renders the
decorative repository animation `fancyUMAP.gif`. It is not a panel of Figure 1,
and it is the only part of the script omitted.

---

## What was changed in the authors' code

`code/Figure1_chunks.R` is the five R chunks of `Figure1.Rmd`, concatenated in
document order. The changes, and nothing else:

1. the two hard-coded paths in the setup chunk read `SCE_PATH` and `OUT_DIR`,
   set by the driver, instead of the author's working-tree paths;
2. `source("../functions.R")` reads the vendored copy in `code/`;
3. one added `ggsave` for `p0`, the faceted UMAP, which the notebook displayed
   inline and never wrote to disk;
4. the `gganimate` block is omitted, as above.

`code/Figure1.R` is a driver, not analysis. It exists because several of the
authors' `ggsave()` calls take no `plot=` argument and so depend on
`last_plot()`, which an R notebook sets by auto-printing top-level expressions
and `Rscript` does not. The driver sources the chunk file with
`print.eval = TRUE` and an open device, which restores exactly that behaviour
without editing their plotting code. The renders it catches are kept in
`results/Figure1/notebook_echo.pdf`.

## Outputs

`results/Figure1/`

| file | chunk |
|---|---|
| `CellTypesUMAP.png` | cell-type UMAP, tumourigenesis + wild type |
| `CellTypesUMAP_FullData.png` | the same over the full integrated atlas |
| `CellTypesUMAP_WithLabels.pdf` | cell-type UMAP with cell-type centroid labels |
| `Barplot_Compartments.svg` | compartment composition per condition |
| `UMAP_ptimeBin_Groups.pdf` | UMAP faceted by tumour-time bin |
| `UMAP_AllTumorigenesis.png`, `UMAP_Bin_{1..5}.png`, `UMAP_Bin_WT_{Young,Old}.png` | per-condition split UMAPs, 12,000 cells sampled per bin under `set.seed(42)` |
| `UMAP_TumorAndPregnancyCellGroups.png`, `UMAP_TumorAndPregnancyExperiment.png` | the full atlas by compartment and by experiment |
| `sessionInfo.txt` | R session of the run |
| `notebook_echo.pdf` | the notebook auto-prints |

Runtime is about 4.5 minutes, nearly all of it reading the 2.0 GB object and
rasterising the PNGs; peak memory is a few GB.
