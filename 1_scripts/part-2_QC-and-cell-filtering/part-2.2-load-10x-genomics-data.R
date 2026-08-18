# ****************************************************************************#
# STEP 2: Load 10x Genomics data
# ****************************************************************************#


# --- DEFINE PATHS/VARIABLES ---
# ****************************************************************************#

# sample metadata definition
# ----------------------------------------------------------------------------#
# method 1 - manual definition (change these variables accordingly)
# (META_SAMPLE_NAME <- "Healthy_1")
# (META_SRA_ID <- "SRR14575500")
# (META_CONDITION <- "Healthy")
# (META_PATIENT_ID <- "Donor_1")


# method 2 - read from file (change this line accordingly)
SAMPLE_ROW <- subset(
  read.delim("2_input/sample-metadata/sample_names.tsv",
    stringsAsFactors = FALSE
  ),
  # choose from sample_name column in "2_input/sample-metadata/sample_names.tsv"
  sample_id == "Healthy_1" # e.g. "Participant_1"
)

META_SAMPLE_NAME <- SAMPLE_ROW$sample_id
META_SRA_ID <- SAMPLE_ROW$sra_id
META_CONDITION <- SAMPLE_ROW$condition
META_PATIENT_ID <- SAMPLE_ROW$patient_id

cat("Sample metadata loaded:\n")
cat("  Sample:   ", META_SAMPLE_NAME, "\n")
cat("  SRA ID:   ", META_SRA_ID, "\n")
cat("  Condition:", META_CONDITION, "\n")
cat("  Patient:  ", META_PATIENT_ID, "\n\n")


# The directory containing the 10x files (from Cell Ranger)
# ----------------------------------------------------------------------------#
CELLRANGER_INPUT <- paste0(
  "2_input/cellranger-matrix-counts/GSE174609_All_Participants/",
  META_SAMPLE_NAME
)


# The date-stamped job run ID (from HPC scheduler) for tracking
# ----------------------------------------------------------------------------#
RUN_ID <- "2026_06_09_brown_job_3058993"


# --- DYNAMIC OUTPUT PATH REGISTRATION ---
# ****************************************************************************#
# Dynamically construct paths inside a date-stamped job run ID.
PLOTS_OUT_DIR <- file.path(
  "3_output", RUN_ID, "qc_and_filtering", "plots", META_SAMPLE_NAME
)
DATA_OUT_DIR <- file.path(
  "3_output", RUN_ID, "qc_and_filtering", "filtered_data", META_SAMPLE_NAME
)
METRICS_OUT_DIR <- file.path(
  "3_output", RUN_ID, "qc_and_filtering", "metrics", META_SAMPLE_NAME
)


# --- AUTOMATED DIRECTORY PROVISIONING ---
# ****************************************************************************#
# Construct the folder trees on disk:
#  - `recursive = TRUE`: Instructs R to build missing parent paths on the fly
#  - `showWarnings = FALSE`: Prevents the script from halting or throwing noise
#    if the directories were already initialized by an earlier script step

dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METRICS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("✓ Output directories securely synchronized to disk:\n")
cat("  • Data Master Archive:", DATA_OUT_DIR, "\n")
cat("  • Pipeline Metrics Ledger:", METRICS_OUT_DIR, "\n")
cat("  • Visual Diagnostic Canvas:", PLOTS_OUT_DIR, "\n\n")


# --- THE CRITICAL MATRIX DECISION: RAW VS FILTERED ---
# ****************************************************************************#
# Cell Ranger (10x software) outputs two entirely distinct count matrices:
#
# 1. RAW (`raw_feature_bc_matrix`):
#    Includes EVERY single droplet that passed through the microfluidic chip
#    (often 100,000+ barcodes). This includes true cells, empty background
#    droplets, ambient "soup" RNA, and doublets.
#
# 2. FILTERED (`filtered_feature_bc_matrix`):
#    Includes only the barcodes that Cell Ranger's internal, conservative
#    algorithm guessed were actual cells (usually ~3,000 - 10,000 cells).
#
# WHY THIS TUTORIAL USES OPTION A (RAW):
#   Cell Ranger's built-in filtering often accidentally throws away real, ultra-
#   small cells (like resting lymphocytes) or retains large dead cell debris.
#   Loading the RAW matrix allows you to run modern, superior downstream tools
#   like EmptyDrops (Step 4), SoupX (Step 5), and scDblFinder (Step 6) to perform
#   customized, high-precision purification that outperforms Cell Ranger.

# OPTION A: Load RAW matrix (what this tutorial demonstrates)
#   - LOG_STEP(<message>, { <code block> }) — prints the start message, runs the
#     block, times it, prints "✓ Done in Xs", and returns the block's result.
#   - Assign it like a normal expression: RESULT <- LOG_STEP("msg...", { ... })
#   - The LOG_STEP() function is loaded in Step 1

COUNTS <- LOG_STEP("Loading 10x matrix (this may take some time)...", {
  Read10X(data.dir = file.path(CELLRANGER_INPUT, "raw_feature_bc_matrix"))
})

# # OPTION B: Load FILTERED matrix (alternative)
# COUNTS <- LOG_STEP("Loading 10x matrix...", {
#   Read10X(data.dir = file.path(CELLRANGER_INPUT, "filtered_feature_bc_matrix"))
# })


# --- CREATE SEURAT OBJECT ---
# ****************************************************************************#
# We initialize the object with ZERO pre-filtering constraints.
# Setting `min.cells = 0` and `min.features = 0` is required because we want to
# keep the ambient background noise intact for now. If we filtered genes or
# cells here, tools like SoupX and EmptyDrops would lack the background
# baseline data they need to estimate contamination accurately.

SEURAT_OBJ <- LOG_STEP("Creating Seurat object (this may take some time)...", {
  CreateSeuratObject(
    counts = COUNTS,
    project = META_SAMPLE_NAME,
    min.cells = 0, # Explicitly do NOT filter genes yet
    min.features = 0 # Explicitly do NOT filter cells yet
  )
})


# --- ADD EXPERIMENTAL METADATA ---
# ****************************************************************************#
# Annotates the dataset with project-specific structural metadata. This is
# essential for downstream analyses when merging multiple patient batches or
# performing multi-condition differential expression (e.g., Healthy vs Disease)

SEURAT_OBJ@meta.data <- SEURAT_OBJ@meta.data %>%
  mutate(
    sample_id  = META_SAMPLE_NAME,
    sra_id     = META_SRA_ID,
    condition  = META_CONDITION,
    patient_id = META_PATIENT_ID,
    time_point = NA_character_ # Explicit typed NA for column consistency
  )


# --- DATASET DIMENSIONS QUICK CHECK ---
# ****************************************************************************#
cat("Loaded:", ncol(SEURAT_OBJ), "droplets ×", nrow(SEURAT_OBJ), "genes\n")
cat("(Most are empty - EmptyDrops will filter in Step 4)\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Before entering this step, we initialized our computational workspace (Step
#   1) by loading our specialized genomic libraries, pinning a global reproduc-
#   ible random seed, and aligning our high-performance Pixi environment
#   packages. At that stage, our project existed only as raw, unparsed Cell
#   Ranger output matrices stored blindly on the cluster disk.
#
# WHAT WE HAVE ACCOMPLISHED:
#   In this step, we successfully completed our data ingestion and environment
#   isolation milestone. We built a dynamic path provisioning framework rooted
#   in our active HPC scheduler `RUN_ID`, creating parallel, sample-isolated
#   folder tracks on disk. We made a deliberate architectural choice to load the
#   raw feature matrix rather than Cell Ranger's pre-filtered matrix,
#   retaining every captured microfluidic barcode with absolute fidelity. By
#   initializing a master Seurat object with zero filtering constraints
#   (`min.cells = 0`, `min.features = 0`) and mutating sample annotations into
#   a single database
#   memory block, we preserved the raw data structure alongside crucial
#   experimental metadata.
#
# WHERE WE ARE HEADING (STEP 3):
#   Our dataset is now successfully instantiated in R, housing over a million
#   captured barcodes alongside our global clinical metadata tags. However, we
#   are currently holding a massive, raw, and un-vetted sequencing dump.
#
#   Because the microfluidic channel captures everything it touches, the
#   overwhelming majority of these million barcodes represent fluid-only, empty
#   droplets or free-floating background molecules. Before we can deploy
#   advanced filtering algorithms, we must understand the baseline noise profile
#   of this specific run. In Step 3, we will perform an initial diagnostic
#   exploration. We will measure the global sparsity of the expression matrix
#   and audit the UMI count distributions to map out the technical footprints of
#   our sample before initiating data purification.
# ****************************************************************************#
