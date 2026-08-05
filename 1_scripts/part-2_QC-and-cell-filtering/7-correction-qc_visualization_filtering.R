# ****************************************************************************#
# STEP 7: Calculate QC metrics, visualize, and filter cells
# ****************************************************************************#


# ****************************************************************************#
# --- THE PURPOSE OF CELL-LEVEL QC ---
# ****************************************************************************#
#   Up until this point, we have performed technical background clean-up
#   (Steps 4, 5, and 6). Step 7 is our final safety net: cell-level Quality
#   Control. Here, we separate viable, healthy cells from technical artifacts
#   and cellular fragments.


# WHY MANUAL INSPECTION BEATS BLIND STATISTICAL AUTOMATION:
#   1. Visibility: You see exactly what you are cutting and why.
#   2. Biological Context: Different tissues possess wildly different baselines
#      (e.g., highly metabolic tumor biopsies or cardiac cells naturally retain
#      elevated mitochondrial fractions that would be auto-deleted by a rigid
#      statistical pipeline).
#   3. Skew Prevention: Outliers can heavily drag down mathematical medians,
#      leading automated methods to accidentally over-filter healthy populations
#
# Always look at your data distribution profiles before committing to your
# cutoffs!


# ****************************************************************************#
# --- CALCULATE CORE DIAGNOSTIC QC METRICS ---
# ****************************************************************************#
#   This block calculates three vital biological sensors used to evaluate
#   whether a cell barcode represents a healthy living cell or dying debris
#
# THE BIOLOGICAL "WHY" BEHIND EACH METRIC:
#   1. percent.mt (Mitochondrial Load): Measures the proportion of reads
#      mapping to mitochondrial genes (^MT-). When a cell's outer membrane
#      ruptures during tissue preparation, its small cytoplasmic RNA leaks
#      out, but the heavy mitochondria get trapped inside. Therefore, a high
#      percentage indicates a dead, burst, or severely stressed cell.
#   2. percent.ribo (Ribosomal Abundance): Scans for structural ribosomal
#      proteins (^RP[SL]). Highly metabolic cells (like active immune cells
#      or dividing tumor cells) naturally translate massive amounts of protein
#      and carry a high ribosomal signal. Sudden drops can signal cell death.
#   3. log10GenesPerUMI (Transcriptional Complexity): Measures how many unique
#      genes are discovered relative to total sequencing depth. A healthy cell
#      is complex and expresses thousands of diverse genes. Debris or empty
#      droplets contain a high count of the same repetitive transcripts,
#      causing this ratio to plunge.

# Calculate per-cell quality control statistics
SEURAT_OBJ@meta.data <- SEURAT_OBJ@meta.data %>%
  mutate(
    percent.mt        = PercentageFeatureSet(SEURAT_OBJ, pattern = "^MT-"),
    percent.ribo      = PercentageFeatureSet(SEURAT_OBJ, pattern = "^RP[SL]"),
    log10GenesPerUMI  = log10(nFeature_RNA) / log10(nCount_RNA)
  )


# ****************************************************************************#
# --- PRINT SUMMARY STATISTICS & DIAGNOSTIC INTERPRETATION ---
# ****************************************************************************#
#   These statistical quantiles and medians serve as initial numeric guardrails
#   to guide your manual threshold settings before plotting.
cat("\nQC Metric Distributions:\n")

# ==============================================================================
# SEQUENCING DEPTH & CAPTURE (nCount_RNA / nFeature_RNA)
#   - High-Quality Target: Medians >3,000 UMIs and >1,200 unique genes indicate
#     robust sequencing depth and exceptional transcript diversity per cell.
#   - Technical Failures: Medians below 1,000 UMIs or 500 genes indicate severe
#     under-sequencing, loading failures, or pervasive RNA degradation.
# ==============================================================================
cat("nCount_RNA (UMI):\n")
cat(
  "  Median:", median(SEURAT_OBJ$nCount_RNA),
  "| Q1-Q3:", quantile(SEURAT_OBJ$nCount_RNA, 0.25), "-",
  quantile(SEURAT_OBJ$nCount_RNA, 0.75), "\n"
)

cat("nFeature_RNA (genes):\n")
cat(
  "  Median:", median(SEURAT_OBJ$nFeature_RNA),
  "| Q1-Q3:", quantile(SEURAT_OBJ$nFeature_RNA, 0.25), "-",
  quantile(SEURAT_OBJ$nFeature_RNA, 0.75), "\n"
)


# ==============================================================================
# 2. MITOCHONDRIAL FRACTIONAL LOAD (percent.mt)
#   - Live Baseline: In a healthy cell, most RNA belongs to the cytoplasm and
#     only a tiny fraction (typically 1% to 5%) belongs to mitochondria. A low
#     median confirms that the vast majority of your data is highly viable.
#   - Membrane Tear Inflection Point: The 95th percentile acts as a technical
#     sensor marking exactly where membranes ruptured and leaked cytoplasmic
#     RNA during tissue processing, splitting the data into two groups:
#       * The 95% Healthy Majority: 95% of all cells kept their membranes
#         intact, retaining a low mitochondrial load between 0% and this score.
#       * The 5% Dead/Dying Tail: The remaining 5% suffered extreme stress;
#         their membranes tore, normal transcripts escaped, and they morphed
#         into "ghost cells" whose mitochondrial signals skyrocketed past
#         this score (potentially up to 50% or 80%).
#   - Rationale: Use this 95th percentile value to anchor your upcoming hard
#     filtering caps (typically 15% to 20%) to cleanly slice away only the
#     necrotic 5% tail.
# ==============================================================================
cat("percent.mt:\n")
cat(
  "  Median:", round(median(SEURAT_OBJ$percent.mt), 2), "%",
  "| 95th percentile:", round(quantile(SEURAT_OBJ$percent.mt, 0.95), 2), "%\n"
)


# ==============================================================================
# 3. RIBOSOMAL TRANSCRIPT ABUNDANCE (percent.ribo)
#   - Translational Baseline: Ribosomes act as a cell's internal protein
#     factories. A normal, high-quality sample typically prints a median score
#     falling safely within the expected 10% to 30% healthy benchmark range.
#   - Metabolic Health Signatures: A robust baseline in this window confirms
#     that the overwhelming majority of your captured singlets represent
#     translationally active, metabolically vibrant cells (a highly positive
#     quality signature classic for immune and expanding tissue lineages).
#   - The Apoptotic/Necrotic Shift: Tracking this score acts as a vital cross-
#     validation tool when evaluating cell death states:
#     * The Co-Dependent Fall: Healthy, quiet cells can naturally have lower
#       ribosomal counts, but their mitochondrial fractions will stay low.
#     * The Necrotic Signpost: If you look at your upcoming scatter plots and
#       locate individual barcodes where percent.ribo plunges toward 0%
#       while percent.mt simultaneously spikes into the extreme tail, you
#       have uncovered a completely dead, lysed cellular carcass ("ghost
#       cell") that must be filtered out.
cat("percent.ribo:\n")
cat("  Median:", round(median(SEURAT_OBJ$percent.ribo), 2), "%\n")


# ****************************************************************************#
# --- GENERATE DISTRIBUTION DIAGNOSTICS ---
# ****************************************************************************#
# THE VISUAL EXAMINATION CRITERIA:
# Open your violin and scatter plots and actively look for these flags:
#   - Bimodal Violins: Double curves indicate distinct populations. Set your
#     filtering line in the "valley" between those two peaks.
#   - High Mitochondrial Tails: Cells pointing up to high MT percentages with
#     low UMI counts are dying cells whose membranes burst. Cut them out!
#   - Outlier Scatter Clusters: Points breaking away from the standard upward
#     diagonal trend line represent aberrant droplets or hidden doublets.

# --- VIOLIN PLOT: QC METRICS DISTRIBUTION ---
P4 <- VlnPlot(
  SEURAT_OBJ,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  ncol     = 4,
  pt.size  = 0.1
)

# Panel 1: nCount_RNA — apply K-scale labels (e.g., 20K, 40K)
P4[[1]] <- P4[[1]] +
  scale_y_continuous(
    breaks = seq(0, max(SEURAT_OBJ$nCount_RNA), by = 20000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  theme(plot.title = element_text(face = "bold"))

# Panel 2: nFeature_RNA — apply K-scale labels (e.g., 2K, 4K)
P4[[2]] <- P4[[2]] +
  scale_y_continuous(
    breaks = seq(0, max(SEURAT_OBJ$nFeature_RNA), by = 2000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  theme(plot.title = element_text(face = "bold"))

# Apply bold title only to panels 3 and 4
P4[[3]] <- P4[[3]] + theme(plot.title = element_text(face = "bold"))
P4[[4]] <- P4[[4]] + theme(plot.title = element_text(face = "bold"))

# Save Visual Diagnostics to Disk
ggsave(
  file.path(PLOTS_OUT_DIR, "03_qc_violins.png"),
  plot = P4, width = 16, height = 4, dpi = 300
)


# --- Scatter Plot Construction ---
# Plot UMI counts against gene counts.
P5 <- FeatureScatter(
  SEURAT_OBJ,
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
) +
  labs(
    title = "UMI vs Genes Detected",
    x = "nCount_RNA (K = 1,000)", y = "nFeature_RNA (K = 1,000)"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(SEURAT_OBJ$nCount_RNA), by = 5000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  scale_y_continuous(
    breaks = seq(0, max(SEURAT_OBJ$nFeature_RNA), by = 1000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )


# Plot UMI counts against mitochondrial percentage
P6 <- FeatureScatter(
  SEURAT_OBJ,
  feature1 = "nCount_RNA", feature2 = "percent.mt"
) +
  labs(
    title = "UMI vs Mitochondrial %",
    x = "nCount_RNA (K = 1,000)", y = "percent.mt (%)"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(SEURAT_OBJ$nCount_RNA), by = 5000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  scale_y_continuous(
    breaks = seq(0, max(SEURAT_OBJ$percent.mt), by = 10)
  ) +
  theme(
    axis.text.x     = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y     = element_text(size = 10),
    legend.position = "none"
  )


# Plot mitochondrial percentage against ribosomal percentage
P7 <- FeatureScatter(
  SEURAT_OBJ,
  feature1 = "percent.mt", feature2 = "percent.ribo"
) +
  labs(
    title = "Mitochondrial % vs Ribosomal %",
    x = "percent.mt (%)", y = "percent.ribo (%)"
  ) +
  scale_x_continuous(breaks = seq(0, max(SEURAT_OBJ$percent.mt), by = 10)) +
  scale_y_continuous(breaks = seq(0, max(SEURAT_OBJ$percent.ribo), by = 5)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
  )


# Assemble all three scatter plots into one figure.
P_SCATTER <- P5 + P6 + P7

# Save Visual Diagnostics to Disk
ggsave(
  file.path(PLOTS_OUT_DIR, "04_qc_scatter.png"),
  plot = P_SCATTER, width = 15, height = 5, dpi = 300
)

cat(
  "\n→ EXAMINE plots 03_qc_violins.png and 04_qc_scatter.png in",
  file.path(PLOTS_OUT_DIR), "\n"
)
cat("→ Look for:\n")
cat("   - Bimodal distributions (good vs bad cells)\n")
cat("   - Outliers in scatter plots\n")
cat("   - Relationship between metrics\n\n")


# ****************************************************************************#
# --- DEFINE THRESHOLDS (MANUAL ADJUSTMENT PHASE) ---
# ****************************************************************************#
# HOW TO SET THESE VALUES FOR ANY DATASET (VISUAL DECODING):
#   1. NFEATURE_MIN / NCOUNT_MIN (The Debris Floor): Look at the UMI vs Genes
#      scatter. Find the dense upward curve. Set your minimums right above the
#      flat, bottom-left "shelf" or "hook" of low-complexity background debris.
#   2. MT_THRESH (The Death Line): Look at the UMI vs Mitochondrial scatter.
#      Locate where the dense horizontal block ends and the vertical wall of
#      dots begins shooting up toward 100%. Set your cap right at that elbow
#      (typically 10-15% for fluid preps, up to 20% for tough solid tumors).
#   3. NFEATURE_MAX / NCOUNT_MAX (The Multiplet Ceiling): Look at the top-right
#      of the UMI vs Genes scatter. Find where the thick, primary cell cloud
#      begins to thin out into sparse, high-flying outlier dots. Place your
#      caps there to prune remaining unflagged physical doublets.

# VERIFICATION CHECK: Cross-validate with the Mitochondrial vs Ribosomal
# scatter. Your chosen MT_THRESH should align with where the ribosomal
# signal collapses.
cat("=== Setting Filtering Thresholds ===\n")
cat("Based on visual inspection of the plots above:\n\n")

# --- Load per-sample filtering thresholds from metadata TSV ---
# NFEATURE_MIN - Min genes required to keep a cell (below this = debris)
# NFEATURE_MAX - Max genes to keep a cell (above this = likely doublets)
# NCOUNT_MIN   - Min UMIs required to keep a cell (below this = debris)
# NCOUNT_MAX   - Max UMIs to keep a cell (above this = likely doublets)
# MT_THRESH    - Cells kept only if mt% is strictly below this value (drops
#                ruptured, dead, or dying cells), adjust based on tissue
#
# Values are set per-sample in sample_names.tsv (edit that file to tune
# thresholds for a given sample based on YOUR visual inspection above).

SAMPLE_ROW <- subset(
  read.delim("2_input/sample-metadata/sample_names.tsv",
    stringsAsFactors = FALSE
  ),
  sample_name == META_SAMPLE_NAME # reuses the sample selected in Step 2
)

if (nrow(SAMPLE_ROW) != 1) {
  stop(
    "Expected exactly 1 row for sample '", META_SAMPLE_NAME,
    "' in sample_names.tsv, found ", nrow(SAMPLE_ROW)
  )
}

NFEATURE_MIN <- SAMPLE_ROW$nfeature_min
NFEATURE_MAX <- SAMPLE_ROW$nfeature_max
NCOUNT_MIN <- SAMPLE_ROW$ncount_min
NCOUNT_MAX <- SAMPLE_ROW$ncount_max
MT_THRESH <- SAMPLE_ROW$mt_thresh

THRESHOLD_VALUES <- c(
  NFEATURE_MIN = NFEATURE_MIN,
  NFEATURE_MAX = NFEATURE_MAX,
  NCOUNT_MIN = NCOUNT_MIN,
  NCOUNT_MAX = NCOUNT_MAX,
  MT_THRESH = MT_THRESH
)

if (any(is.na(THRESHOLD_VALUES))) {
  stop(
    "Missing threshold value(s) for sample '", META_SAMPLE_NAME,
    "' in sample_names.tsv: ",
    paste(names(THRESHOLD_VALUES)[is.na(THRESHOLD_VALUES)], collapse = ", ")
  )
}

cat("Set thresholds (adjust in sample_names.tsv based on YOUR data):\n")
cat("  nFeature_RNA: [", NFEATURE_MIN, ",", NFEATURE_MAX, "]\n")
cat("  nCount_RNA: [", NCOUNT_MIN, ",", NCOUNT_MAX, "]\n")
cat("  percent.mt: <", MT_THRESH, "%\n\n")


# ****************************************************************************#
# --- GENERATE CUTOFF THRESHOLD PLOTS ---
# ****************************************************************************#
# This block executes our data-driven quality gates, tagging cells as TRUE
# (viable) or FALSE (debris/dead) based on the criteria set in Step 4.
#
# TECHNICAL BREAKDOWN OF THE EVALUATION LOGIC:
#   1. Data Extraction: Rips out the core numeric metrics from Seurat and
#      places them into a fast, isolated temporary dataframe (QC_DF).
#   2. Multi-Gate Classification (pass_qc): Every cell must pass all 5 gates
#      simultaneously. Failing even one gate drops it into the discard pool.
#   3. Log-Scale Synchronization: Converts our raw counts into log10 space
#      to match standard genomic density plots, overlaying our hard choices
#      as an explicit red box boundary.
#
# VISUAL NOTE (THE HIDDEN DIMENSION):
#   - If you see RED dots trapped inside the green bounding box in your saved
#     plot, that is correct behavior! Those cells pass the gene/UMI count caps,
#     but failed the hidden 3rd dimension: their mitochondrial load exceeded
#     10%.


# Create a temporary dataframe to store the QC metrics and filtering results
QC_DF <- SEURAT_OBJ@meta.data %>%
  # Extract the QC metrics from the Seurat object
  select(nCount_RNA, nFeature_RNA, percent.mt) %>%
  # Generate the pass_qc column based on the 5-point logical check
  mutate(
    pass_qc = nCount_RNA >= NCOUNT_MIN & # keep cells with at least $NCOUNT_MIN UMIs
      nCount_RNA <= NCOUNT_MAX & # keep cells with at most $NCOUNT_MAX UMIs
      nFeature_RNA >= NFEATURE_MIN & # keep cells with at least $NFEATURE_MIN genes
      nFeature_RNA <= NFEATURE_MAX & # keep cells with at most $NFEATURE_MAX genes
      percent.mt < MT_THRESH # keep cells with less than $MT_THRESH% mito. reads
  )


# Plot the 2D distribution with boundary overlays
# --- Option A: Log-Scale Bounding Box Plot ---
P8_LOG <- ggplot(QC_DF, aes(
  x = log10(nCount_RNA + 1), y = log10(nFeature_RNA + 1), color = pass_qc
)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = log10(c(NCOUNT_MIN, NCOUNT_MAX)),
    linetype = "dashed", color = "red"
  ) +
  geom_hline(
    yintercept = log10(c(NFEATURE_MIN, NFEATURE_MAX)),
    linetype = "dashed", color = "red"
  ) +
  scale_color_manual(values = c("TRUE" = "#06D6A0", "FALSE" = "#EF476F")) +
  labs(
    title = "Cell Filtering Thresholds (Log Scale)",
    subtitle = paste0(
      sum(QC_DF$pass_qc), " cells pass QC (",
      round(sum(QC_DF$pass_qc) / nrow(QC_DF) * 100, 1), "%)"
    ),
    x = "log10(UMI + 1)",
    y = "log10(Genes + 1)",
    color = "Pass QC"
  ) +
  theme_classic()


# --- Option B: Linear-Scale Bounding Box Plot ---
P8_LINEAR <- ggplot(QC_DF, aes(
  x = nCount_RNA, y = nFeature_RNA, color = pass_qc
)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = c(NCOUNT_MIN, NCOUNT_MAX),
    linetype = "dashed", color = "red"
  ) +
  geom_hline(
    yintercept = c(NFEATURE_MIN, NFEATURE_MAX),
    linetype = "dashed", color = "red"
  ) +
  scale_color_manual(values = c("TRUE" = "#06D6A0", "FALSE" = "#EF476F")) +
  scale_x_continuous(
    breaks = seq(0, max(QC_DF$nCount_RNA), by = 5000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  scale_y_continuous(
    breaks = seq(0, max(QC_DF$nFeature_RNA), by = 1000),
    labels = scales::label_number(scale_cut = cut_short_scale())
  ) +
  labs(
    title = "Cell Filtering Thresholds (Linear Scale)",
    subtitle = paste0(
      sum(QC_DF$pass_qc), " cells pass QC (",
      round(sum(QC_DF$pass_qc) / nrow(QC_DF) * 100, 1), "%)"
    ),
    x = "\nTotal Transcripts (UMI counts, K = 1,000)",
    y = "Unique Genes Detected (K = 1,000)\n",
    color = "Pass QC"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# --- Save Visual Diagnostics to Disk ---
ggsave(
  file.path(PLOTS_OUT_DIR, "05_filtering_thresholds_log.png"),
  plot = P8_LOG, width = 8, height = 7, dpi = 300
)

# Save Visual Diagnostics to Disk
ggsave(
  file.path(PLOTS_OUT_DIR, "05_filtering_thresholds_linear.png"),
  plot = P8_LINEAR, width = 8, height = 7, dpi = 300
)

# ****************************************************************************#
# --- REPORT FILTERING METRICS ---
# ****************************************************************************#
# This block audits your threshold configurations, generating a granular report
# of exactly how many cells fell victim to each individual technical filter.
#
# INTERPRETATION FRAMEWORK (THE DIAGNOSTIC AUDIT):
#   - Total Attrition: Tracks your overarching data loss. If you lose most of
#     your data here, downstream clusters will lack statistical power.
#   - Overlapping Casualties: Note that the granular sum can exceed total cells
#     removed. Damaged cells frequently fail multiple gates simultaneously
#     (e.g., a ruptured cell often has low genes AND high mitochondrial load).
#   - Tissue Baseline Check: Look at which gate caught the most cells. A massive
#     spike in the "High MT%" bucket indicates bad tissue dissociation stress,
#     whereas high gene counts indicate an abundance of physical doublets.

REMOVAL_PCT <- sum(!QC_DF$pass_qc) / nrow(QC_DF) * 100
cat("Filtering impact:\n")
cat("  Cells before:", nrow(QC_DF), "\n")
cat("  Cells passing QC:", sum(QC_DF$pass_qc), "\n")
cat(
  "  Cells removed:", sum(!QC_DF$pass_qc),
  "(", round(REMOVAL_PCT, 1), "%)\n\n"
)

cat("Removal breakdown:\n")
cat(
  "  Low genes (<", NFEATURE_MIN, "):",
  sum(QC_DF$nFeature_RNA < NFEATURE_MIN), "\n"
)
cat(
  "  High genes (>", NFEATURE_MAX, "):",
  sum(QC_DF$nFeature_RNA > NFEATURE_MAX), "\n"
)
cat(
  "  High MT% (>", MT_THRESH, "%):",
  sum(QC_DF$percent.mt > MT_THRESH), "\n\n"
)


# ****************************************************************************#
# --- AUTOMATED SAFETY GUARD CHECK ---
# ****************************************************************************#
# This programmatic guardrail prevents manual data manipulation errors. It
# cross-references your total attrition percentage against universal genomic
# standards.
#
# BIOLOGICAL RATIONALE FOR THE GUARDRAILS:
#   - The Over-Filtering Danger (>30%): Slicing away nearly a third of your data
#     suggests your thresholds are too aggressive, or your tissue preparation
#     suffered catastrophic survival rates. It risks erasing whole cell types.
#   - The Under-Filtering Danger (<5%): Removing almost nothing means you are
#     retaining low-quality cell fragments, empty ambient soup droplets, and
#     doublets, which will corrupt downstream clustering with technical noise.
#   - The Healthy Corridor (10-25%): The sweet spot for standard single-cell
#     runs that ensures clean cell states without sacrificing biological
#     diversity.

if (REMOVAL_PCT > 30) {
  cat("WARNING: Removing", round(REMOVAL_PCT, 1), "% of cells is high!\n")
  cat("   Consider relaxing thresholds (especially MT%)\n\n")
} else if (REMOVAL_PCT < 5) {
  cat("WARNING: Only removing", round(REMOVAL_PCT, 1), "% of cells is low.\n")
  cat("   Check if you're retaining low-quality cells.\n\n")
} else {
  cat("✓ Removal rate is reasonable (typical: 10-25%)\n\n")
}

# Cache final diagnostic tracker for historical workflow tracking
CELLS_BEFORE_CELL_QC <- nrow(QC_DF)


# ****************************************************************************#
# --- APPLY SUBSETTING FILTERS ---
# ****************************************************************************#
# Execution Phase: This is the point of no return. Up until this step, all
# numbers were just passive visual predictions inside a temporary dataframe.
#
# TECHNICAL IMPACT ON SEURAT STORAGE:
#   - Hard Memory Erasure: The subset() function permanently deletes the low-
#     quality barcode columns from your active Seurat object's gene expression
#     matrices, metadata, and data slots.
#   - Downstream Readiness: Purging these dead cells compresses the object size
#     in your computer's RAM, preparing a clean, pristine environment for
#     accurate normalisation, variable gene discovery, and dimensional
#     reduction.

SEURAT_OBJ <- subset(
  SEURAT_OBJ,
  subset = nCount_RNA >= NCOUNT_MIN &
    nCount_RNA <= NCOUNT_MAX &
    nFeature_RNA >= NFEATURE_MIN &
    nFeature_RNA <= NFEATURE_MAX &
    percent.mt < MT_THRESH
)

cat("After filtering:", ncol(SEURAT_OBJ), "cells remaining\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
# Before entering this step, we had processed our dataset through automated
# background clean-up filters (Steps 4, 5, and 6), which successfully left us
# with a refined matrix of singlet droplets. However, these droplets still
# contained an unknown mixture of highly viable living cells, dead or dying
# cellular debris ("ghost cells"), and low-complexity technical fragments.
#
# WHAT WE HAVE ACCOMPLISHED:
# In this step, we successfully established our final and most critical cell-
# level safety net through data-driven Quality Control (QC). By mapping out
# core biological sensors—such as mitochondrial load (percent.mt) and
# transcriptional complexity (log10GenesPerUMI)—we performed a detailed visual
# audit of our cell distributions. Instead of blindly trusting automated caps,
# we drew precise, custom threshold lines directly through the cell clouds
# on our diagnostic log and linear scatter plots. This multi-gate filter allowed
# us to confidently slice away the necrotic, high-mitochondrial 5% tail and
# low-complexity debris, ensuring our active cellular matrix represents strictly
# healthy, metabolically viable singlet systems.
#
# WHERE WE ARE HEADING (STEP 8):
# Our dataset has reached its definitive cell count, compressed in memory and
# purified of low-quality barcodes. However, while our cells are clean, our
# gene dimensions are still completely unvetted and cluttered with noise.
#
# Because thousands of genes are completely inactive in our targeted tissue or
# were only captured in a tiny handful of cells, they act as massive, un-
# informative computational weight that distorts statistical averages. In Step
# 8, we will transition to Gene-Level Quality Control. We will calculate the
# global prevalence profile of every gene across our verified cell pool and
# apply a clean percentage threshold (keeping only genes detected in >= 0.1%
# of cells) to purge rare artifacts and zero-count columns from our active
# feature space.
# ****************************************************************************#
