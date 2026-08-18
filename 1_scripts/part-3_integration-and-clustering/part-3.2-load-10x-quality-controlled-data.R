#-----------------------------------------------
# STEP 2: Load QC-filtered samples from Part 2
#-----------------------------------------------

# Define sample metadata
sample_metadata <- data.frame(
  sample_id = c(
    "Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4",
    "Periodontitis_Post_1", "Periodontitis_Post_2", "Periodontitis_Post_3", "Periodontitis_Post_4"
  ),
  condition = c(
    rep("Healthy", 4),
    rep("Periodontitis_Post_Treatment", 4)
  ),
  patient_id = c(
    "Donor_1", "Donor_2", "Donor_3", "Donor_4",
    "Patient_1", "Patient_2", "Patient_3", "Patient_4"
  ),
  stringsAsFactors = FALSE
)

# Load QC-filtered Seurat objects
qc_data_path <- "3_output/2026_06_09_brown_job_3058993/qc_and_filtering/filtered_data"


# Match this to the cores you requested in your srun/sbatch job (see the
# `srun` task in pixi.toml — bump --cpus-per-task there if you want more).
plan(multisession, workers = min(length(sample_metadata$sample_id), availableCores()))

seurat_list <- future_map(
  sample_metadata$sample_id,
  function(sample_id) {
    obj <- readRDS(file.path(qc_data_path, sample_id, paste0(sample_id, "_qc_filtered.rds")))
    obj$sample_id <- sample_id
    obj$condition <- sample_metadata$condition[sample_metadata$sample_id == sample_id]
    obj$patient_id <- sample_metadata$patient_id[sample_metadata$sample_id == sample_id]
    obj
  },
  .options = furrr_options(seed = TRUE)  # equivalent to future.seed = TRUE
)

plan(sequential)  # release workers once done, good practice before heavy in-process steps

names(seurat_list) <- sample_metadata$sample_id

# Ensure all samples have the same genes (critical for FastMNN integration)
all_genes <- lapply(seurat_list, rownames)
common_genes <- Reduce(intersect, all_genes)

# Subset all samples to common genes
seurat_list <- lapply(seurat_list, function(obj) {
  obj[common_genes, ]
})

# Report dimensions
cat("\nLoaded", length(seurat_list), "QC-filtered samples:\n")
for (sample_id in names(seurat_list)) {
  cat(sprintf("  %s: %d cells × %d genes\n",
              sample_id,
              ncol(seurat_list[[sample_id]]),
              nrow(seurat_list[[sample_id]])))
}

cat("\nTotal cells across all samples:",
    sum(sapply(seurat_list, ncol)), "\n")