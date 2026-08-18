# ****************************************************************************#
# Installation: Integration-specific packages not available via pixi/conda
# ****************************************************************************#

# --- THE INFRASTRUCTURE REQUIREMENTS ---
# As with Part 2, every integration package with a conda-forge or bioconda
# build is declared in `pixi.toml` (`part3-feature`) and installed via
# `pixi install` — Harmony, batchelor, scran, ggrepel, RColorBrewer, viridis,
# remotes, future, FNN, cluster, reshape2. Do NOT add `install.packages()`
# calls here for anything already declared there.
#
# The ONE exception is `SeuratWrappers`, which is distributed only via GitHub
# (satijalab/seurat-wrappers) with no CRAN or conda release, so it must be
# installed at the R level using `remotes` (itself pixi-managed).

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    remotes::install_github("satijalab/seurat-wrappers")
}

cat("\n✓ SeuratWrappers verified. All other packages are managed by pixi.toml.\n")
