# ****************************************************************************#
# Installation: Integration-specific packages not available via pixi/conda
# ****************************************************************************#

# --- THE INFRASTRUCTURE REQUIREMENTS ---
# Every integration package with a conda-forge or bioconda build is declared
# in `pixi.toml` — shared tooling (Seurat, dplyr, ggplot2, ggrepel,
# RColorBrewer, viridis, reshape2, remotes, future) in `core-feature`, and
# integration-specific packages (Harmony, batchelor, scran, FNN, cluster) in
# `part3-feature`. Both are installed via `pixi install -e part3`. Do NOT add
# `install.packages()` calls here for anything already declared there.
#
# What's left are two packages with no conda package at all, so pixi's
# channels can't reach them, and they must be installed at the R level:
#   1. `colorout`, hosted on the independent R-specific repository
#      community.r-multiverse.org (same package/exception as Part 2).
#   2. `SeuratWrappers`, distributed only via GitHub
#      (satijalab/seurat-wrappers), installed using `remotes` (itself
#      pixi-managed).

# REPRODUCIBILITY & ECOSYSTEM BENEFIT:
# By wrapping each install inside an `if (!requireNamespace(...))` conditional
# structure, the script validates your local environment first. It only
# downloads a package if it is missing, preventing your script from wasting
# compute hours and bandwidth re-downloading it during repetitive runs.

# --- 1. colorout ---
if (!requireNamespace("colorout", quietly = TRUE)) {
  install.packages("colorout", repos = "https://community.r-multiverse.org")
}

# --- 2. SeuratWrappers ---
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  # `ref`: pin to a specific commit SHA (not @HEAD) so re-running this script
  #   in a year installs the exact same code, not whatever's newest on the
  #   default branch. Find/update the SHA at:
  #   https://github.com/satijalab/seurat-wrappers/commits/master
  # `upgrade = "never"`: stops remotes from checking SeuratWrappers'
  #   dependencies against CRAN and offering to upgrade them. Without this,
  #   it will happily upgrade Seurat, future, rlang, etc. to whatever's
  #   newest on CRAN — silently overwriting the exact versions pixi.toml
  #   pinned and defeating the whole point of using pixi.
  remotes::install_github(
    "satijalab/seurat-wrappers",
    ref = "8df8343", # replace with the commit SHA you've validated against
    upgrade = "never"
  )
}

cat("\n✓ colorout and SeuratWrappers verified. All other packages are managed by pixi.toml.\n")
