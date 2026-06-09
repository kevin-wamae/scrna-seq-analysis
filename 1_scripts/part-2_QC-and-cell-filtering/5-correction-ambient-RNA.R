# -----------------------------------------------------------------------------#
# STEP 5: Ambient RNA correction with SoupX
# -----------------------------------------------------------------------------#

# --- THE PROBLEM OF AMBIENT RNA ("THE SOUP") ---
# Even droplets that successfully capture a healthy, living cell are not pure.
# They are submerged in an extracellular liquid suspension saturated with free-
# floating RNA molecules that leaked out of fragile cells that lysed (burst)
# during tissue dissociation.
#
# WHY AMBIENT RNA MATTERS:
# This background "soup" is pulled into EVERY droplet alongside your true cells.
# It acts as a universal contaminant, creating false-positive gene expression.
# It severely inflates marker genes and heavily skews tissue-specific profiles
# (e.g., free-floating Hemoglobin flooding immune cells in blood, or Myelin
# transcripts contaminating non-neuronal cells in brain tissue). If left
# uncorrected, it misleads downstream clustering and compromises the
# statistical specificity of your differential expression analysis.
#
# HOW SOUPX RESOLVES IT:
# SoupX calculates the exact expression fingerprint of the empty droplets,
# profiles your cell-containing droplets, and uses advanced regression modeling
# to figure out exactly how much ambient background noise needs to be subtracted
# out of your true expression counts.

# --- WORKFLOW CONDITIONAL PREREQUISITE ---
# NOTE: If you chose to load a Cell Ranger FILTERED matrix instead of a RAW
# matrix in Step 2, you MUST skip this step. Cell Ranger has already applied
# its own basic ambient correction. Skip directly to Step 6 (Doublet Detection).


# --- 1. Prepare SoupX Channel ---
# We feed SoupX two distinct components:
# - `tod` (Table of Droplets): The completely raw, unfiltered count matrix
#   containing ALL barcodes (saved in Step 4) used to establish the soup profile.
# - `toc` (Table of Cells): The purified count matrix containing only the
#   validated cell barcodes that passed Step 4's EmptyDrops filter.
tod <- raw_counts_all_droplets
toc <- LayerData(seurat_obj, layer = "counts")
sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)


# --- 2. Rapid Clustering for Background Calibration ---
# SoupX needs an approximate understanding of cell identity to figure out which
# genes do not biologically belong in specific cells. For example, if an
# entire cluster of T-cells expresses 3% Hemoglobin, SoupX knows that is
# biological nonsense and tags Hemoglobin as ambient contamination.
# We utilize standard Seurat pipelines linked via dplyr piping syntax (`%>%`)
# to rapidly cluster the cells and feed the assignments directly into SoupX.

# PIPELINE STEP BREAKDOWN:
#   1. NormalizeData: Scales counts per cell to 10k and log-transforms them
#   2. FindVariableFeatures: Selects top 2,000 variance-driving genes
#   3. ScaleData: Shifts gene expression means to 0 and scales variance to 1
#   4. RunPCA: Compresses 2,000 variable genes down to 30 principal components
#   5. FindNeighbors: Builds an SNN graph connecting cells with similar PCs
#   6. FindClusters: Uses the Louvain algorithm to partition the SNN graph

temp_obj <- seurat_obj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)

# Map assignments: Link cluster numbers with barcodes and pass to SoupX
sc <- setClusters(
  sc,
  setNames(as.character(temp_obj$seurat_clusters), colnames(temp_obj))
)


# --- 3. Estimate Contamination Levels ---
# The `autoEstCont` function looks for highly cell-type-specific genes that
# should be completely silent in most other cells, monitoring how far they
# spread across the whole dataset to automatically compute a global
# contamination fraction (Rho).
# Wrapping this in a `tryCatch` block ensures that if your sample is so clean
# that the math cannot isolate a background signal, the code won't crash.
sc <- tryCatch(
  {
    autoEstCont(sc, verbose = FALSE)
  },
  error = function(e) {
    cat("Note: autoEstCont failed to find contamination signatures\n")
    sc$fit$rho <- NULL
    return(sc)
  }
)


# Save the estimated fraction of ambient RNA contamination
contamination_fraction <- sc$fit$rho


# --- 4. Evaluating Running Metrics & Ranges ---
# EXPECTED RESULTS & THRESHOLDS:
# - Contamination < 5% (or NULL): Extremely clean prep (common in PBMCs or
#   well-preserved suspensions). No correction is needed; skip to avoid noise.
# - Contamination 5% to 15%: Standard single-cell run. Worth correcting to
#   sharpen your downstream cluster markers.
# - Contamination > 20%: High cell lysis sample (frequent in fibrous solid
#   tumors or frozen brain biopsies). Correction is absolutely vital here.

cat(
  "Contamination fraction:",
  ifelse(is.null(contamination_fraction), "NULL (very clean sample)",
    paste0(round(contamination_fraction * 100, 2), "%")
  ), "\n"
)

# --- DECODING A "NULL" CONTAMINATION FRACTION ---
# WHAT A "NULL" RESULT MEANS:
# If `contamination_fraction` returns `NULL`, it indicates that the automated
# estimation engine (`autoEstCont`) could not locate a statistically viable
# ambient background signature across your cell clusters.
#
# THIS HAPPENS FOR TWO MAIN BIOLOGICAL/TECHNICAL REASONS:
#   1. High Sample Integrity: The tissue dissociation was exceptionally gentle,
#      resulting in negligible cell lysis (bursting). There is simply no ambient
#      "soup" floating in the suspension to be detected.
#   2. Lack of Baseline Markers: The dataset may lack clear, hyper-specific
#      marker genes (like Hemoglobin in blood or Myelin in brain tissue) that
#      the algorithm relies on to anchor and measure ambient leakage.
#
# HOW TO HANDLE IT IN ANY DATASET:
# A `NULL` result is an excellent sign. It means your sample is cleanly preserved
# and does not require background subtraction. Attempting to force an arbitrary
# correction value on a `NULL` sample risks introducing artificial noise and
# eroding genuine biological transcripts. Safely skip correction and proceed.


# --- 5. Data-Driven Correction Decision ---
# WHY WE USE A THRESHOLD-BASED DECISION GATE:
# Running SoupX is not a mandatory box-checking exercise; it is an active
# intervention. If your contamination fraction returns as `NULL` (meaning the
# algorithm could not find enough background noise to model a signature) or
# falls below 5% (<0.05), your sample is exceptionally clean.
#
# THE RISK OF OVER-CORRECTION:
# Attempting to force background subtraction on a pristine sample is highly
# dangerous. Without a clear ambient signal to target, `adjustCounts` can
# begin eroding genuine, low-abundance biological transcripts, introducing
# artificial processing noise into your clean data.
#
# THE DECISION MATRIX:
#   - Fraction < 5% or NULL -> Skip correction to preserve true raw biology.
#   - Fraction >= 5%        -> Apply correction to scrub ambient noise.

if (is.null(contamination_fraction) || contamination_fraction < 0.05) {
  cat("→ Contamination is NULL or <5% - skipping SoupX correction\n")
  cat("   Your sample is clean! Proceeding with original counts.\n")
} else {
  cat("→ Contamination is", round(contamination_fraction * 100, 2), "% (≥5%)\n")
  cat("   Applying SoupX correction...\n")

  # `adjustCounts` safely subtracts the background ambient counts from the
  # expression matrix and outputs a cleaned, integer-rounded matrix.
  suppressWarnings({
    corrected_counts <- adjustCounts(sc)
  })

  # Inject the corrected matrix right back into the Seurat object's count layer
  seurat_obj <- SetAssayData(
    seurat_obj,
    layer = "counts",
    new.data = corrected_counts
  )
  cat("   ✓ SoupX correction applied\n")
}


# -----------------------------------------------------------------------------#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# -----------------------------------------------------------------------------#
# WHERE WE STARTED:
# Before entering this step, we used EmptyDrops (Step 4) to filter our massive
# raw matrix down to ~9,600 validated cell barcodes. However, those cells were
# still flagged for potential extracellular contamination from cell-free
# background ambient RNA.
#
# WHAT WE HAVE ACCOMPLISHED:
# We calculated the expression profile of the ambient background and generated
# a rapid clustering reference for SoupX. By evaluating the global contamination
# fraction, our data-driven decision loop determined that our sample is
# completely clean (returning a fraction of NULL/<5%). Consequently, we safely
# skipped background subtraction to avoid introducing processing artifacts,
# successfully preserving our pristine, unmanipulated raw expression counts.
#
# WHERE WE ARE HEADING (STEP 6):
# While our transcript expressions are completely free of ambient background
# soup noise, we still face a significant physical hazard: multi-cell
# overlapping captures within individual droplets.
#
# In Step 6, we will pass our preserved counts to `scDblFinder`. The algorithm
# will generate synthetic doublet profiles to train a machine learning model,
# allowing us to scan our droplets, identify structural multiplets, and
# systematically strip out physical doublets to ensure our downstream results
# represent individual cells.
# -----------------------------------------------------------------------------#
