# ****************************************************************************#
# Installation: Required packages for single-sample Quality Control
# ****************************************************************************#

# --- THE INFRASTRUCTURE REQUIREMENTS ---
# Single-cell QC depends on tools spanning different software ecosystems. As of
# this pipeline's current setup, those pools are handled by two different
# managers, on purpose:
#   1. pixi (pixi.toml): The single source of truth for every package with a
#      conda-forge or bioconda build — e.g. Seurat, SeuratObject, SoupX,
#      DropletUtils, scDblFinder, SingleCellExperiment, ggplot2, patchwork,
#      dplyr, scales, etc. These are version-pinned in `pixi.toml` and
#      installed via `pixi install`. Do NOT add `install.packages()` calls
#      here for anything already declared there — that would let this script's
#      unpinned, "whatever's current on CRAN today" install silently drift
#      away from the version pixi resolved, defeating the point of pinning it.
#   2. External/Community (this script): Packages with no conda package at
#      all, so pixi's channels can't reach them, and they must be installed at
#      the R level instead. Currently just `colorout`, hosted on the
#      independent R-specific repository community.r-multiverse.org.
#
# REPRODUCIBILITY & ECOSYSTEM BENEFIT:
# By wrapping the install inside an `if (!requireNamespace(...))` conditional
# structure, the script validates your local environment first. It only
# downloads the package if it is missing, preventing your script from wasting
# compute hours and bandwidth re-downloading it during repetitive runs.

if (!requireNamespace("colorout", quietly = TRUE)) {
  install.packages("colorout", repos = "https://community.r-multiverse.org")
}

cat("\n✓ colorout verified. All other packages are managed by pixi.toml.\n")
