# ****************************************************************************#
# STEP 3.0 :Installation of R-level packages (not available via pixi/conda)
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   None — this is the first step of the pipeline. It only installs R-level
#   packages (colorout, SeuratWrappers) and is safe to re-run. The dependency
#   check below confirms nothing else is required before continuing.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-3.0-install-packages.R")

# --- INFRASTRUCTURE BOUNDARY ---
# Packages on conda-forge/bioconda (Seurat, Harmony, ggplot2, etc.) are managed 
# by `pixi.toml` (`core-feature` + `part3-feature`). Do NOT add install.packages() 
# for them here. Only packages missing from Conda are installed at the R level:
#   1. colorout       - R terminal colorizer (from community.r-multiverse.org)
#   2. SeuratWrappers - Community tools (GitHub-only via remotes)

# --- 1. colorout ---
if (!requireNamespace("colorout", quietly = TRUE)) {
  install.packages("colorout", repos = "https://community.r-multiverse.org")
}

# --- 2. SeuratWrappers ---
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  # `ref`: Commit SHA locks the code version for future reproducibility.
  # `upgrade = "never"`: Prevents remotes from checking CRAN and overriding
  # versions pinned by pixi.toml (e.g., Seurat, rlang, future).
  remotes::install_github(
    "satijalab/seurat-wrappers",
    ref = "8df8343",
    upgrade = "never"
  )
}

cat("\n✓ colorout and SeuratWrappers verified. All other packages are managed by pixi.toml.\n")