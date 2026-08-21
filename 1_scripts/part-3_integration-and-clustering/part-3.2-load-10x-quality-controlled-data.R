# ****************************************************************************#
# STEP 2: Load QC-filtered samples from Part 2
# ****************************************************************************#


# --- LOAD SAMPLE METADATA (SINGLE SOURCE OF TRUTH) ---
# ****************************************************************************#
# Rather than re-declaring sample identities in a separate hardcoded
# data.frame (which drifts out of sync the moment Part 2's TSV changes), we
# read from the exact same `sample_names.tsv` that Step 2 of Part 2 uses to
# select each sample's Cell Ranger input and QC thresholds. This guarantees
# the condition/patient labels attached here always match what was actually
# used to generate the `.rds` files we're about to load.
#
# The TSV contains all 12 cohort samples (Healthy, Periodontitis_Pre, and
# Periodontitis_Post). For this integration run we deliberately exclude the
# Pre-Treatment timepoint, keeping only Healthy vs Post-Treatment — adjust
# the `filter()` below if your comparison changes.

# Read the sample metadata
sample_metadata <- read.delim("2_input/sample-metadata/sample_names.tsv",
    stringsAsFactors = FALSE) %>%
    select(sample_id, condition, patient_id) %>%
    filter(condition != "Periodontitis_Pre_Treatment")

# Sanity check that the metadata looks right
cat("Sample metadata loaded:\n")
cat("  Samples:   ", nrow(sample_metadata), "\n")
cat("  Conditions:", paste(unique(sample_metadata$condition), collapse = ", "), "\n\n")


# --- LOCATE PART 2 OUTPUT ---
# ****************************************************************************#
# Points at the same date-stamped `RUN_ID` job folder that Part 2's
# `10-save-clean-dataset.R` wrote each sample's `_qc_filtered.rds` and
# `cell_metadata.csv` into (see the `filtered_data/<sample_id>/` tree).
qc_data_path <- "3_output/2026_06_09_brown_job_3058993/qc_and_filtering/filtered_data"


# --- PARALLEL LOAD OF QC-FILTERED SEURAT OBJECTS ---
# ****************************************************************************#
# Reading 8+ compressed `.rds` files sequentially is I/O-bound rather than
# CPU-bound, so we parallelize across samples with `future`/`furrr` instead
# of a plain `lapply()`. `.options = furrr_options(seed = TRUE)` ensures any
# incidental RNG use inside Seurat's deserialization (e.g. tempfile naming)
# gets a proper, parallel-safe RNG stream per worker rather than triggering
# an "UNRELIABLE VALUE" warning.
#
# Match `workers` to the cores requested in your srun/sbatch job (see the
# `srun` task in pixi.toml — bump `--cpus-per-task` there if you want more).

# Set up parallel workers for loading the Seurat objects
plan(multisession, workers = min(nrow(sample_metadata), availableCores()))

# Load the Seurat objects in parallel
seurat_list <- LOG_STEP("Loading QC-filtered Seurat objects in parallel...", {
  future_map(
    sample_metadata$sample_id,
    function(sample_id) {
      obj <- readRDS(file.path(qc_data_path, sample_id, paste0(sample_id, "_qc_filtered.rds")))
      # Re-stamp metadata from the TSV, rather than trusting whatever each
      # sample's saved object already carries — keeps this script the single
      # authority on which condition/patient labels feed into integration.
      obj$sample_id <- sample_id
      obj$condition <- sample_metadata$condition[sample_metadata$sample_id == sample_id]
      obj$patient_id <- sample_metadata$patient_id[sample_metadata$sample_id == sample_id]
      obj
    },
    .options = furrr_options(seed = TRUE) # equivalent to future.seed = TRUE
  ) %>%
    set_names(sample_metadata$sample_id)
})

# Release workers once done, good practice before heavy in-process steps
plan(sequential)


# --- ENFORCE A SHARED GENE SPACE ACROSS SAMPLES ---
# ****************************************************************************#
# WHY THIS MATTERS FOR INTEGRATION:
# Each sample's Part 2 gene-level QC filter (Step 8) can retain a slightly
# different gene set — a gene detected in ≥0.1% of cells in one sample might
# fall just below that bar in another. FastMNN (via Seurat's IntegrateLayers)
# and Harmony both require every sample to share the exact same feature
# (gene) rows before batch correction; a mismatched gene space either errors
# outright or silently misaligns the shared PCA space integration depends on.
#
# We take the intersection of every sample's gene list, then subset all
# samples down to that common set.

# Get the common genes across all samples
common_genes <- seurat_list %>%
  map(rownames) %>%
  reduce(intersect)

# Subset all samples to the common gene set
seurat_list <- LOG_STEP("Subsetting all samples to common gene set...", {
  seurat_list %>%
    map(~ .x[common_genes, ])
})


# --- DATASET DIMENSIONS QUICK CHECK ---
# ****************************************************************************#
cat("\nLoaded", length(seurat_list), "QC-filtered samples:\n")
iwalk(seurat_list, ~ cat(sprintf(
  "  %s: %d cells × %d genes\n", .y, ncol(.x), nrow(.x)
)))

cat("\nTotal cells across all samples:", sum(map_int(seurat_list, ncol)), "\n")
cat("Common genes retained:", length(common_genes), "\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Before entering this step, Part 2 produced 12 independently QC-filtered,
#   per-sample Seurat objects — each cleaned of empty droplets, ambient RNA,
#   doublets, and dead cells, but living in isolation on disk with no shared
#   structure between them.
#
# WHAT WE HAVE ACCOMPLISHED:
#   In this step, we re-established `sample_names.tsv` as the single source
#   of truth for sample identity, filtered our cohort down to the
#   Healthy vs Post-Treatment comparison, and loaded the corresponding
#   `.rds` files in parallel via `future`/`furrr`. We then resolved the one
#   structural obstacle standing between "a folder of separate objects" and
#   "an integratable dataset": each sample's independently-filtered gene set.
#   By intersecting all sample gene lists into a single common feature space,
#   every object entering Step 3 now shares identical rows.
#
# WHERE WE ARE HEADING (STEP 3):
#   With a common gene space established, our sample list is now
#   structurally ready for batch correction. In Step 3, we will merge these
#   objects and run integration methods (Harmony/FastMNN) to correct for
#   technical, patient-driven variation between samples, allowing biological
#   differences between Healthy and Periodontitis states to emerge in
#   downstream clustering rather than being masked by batch effects.
# ****************************************************************************#