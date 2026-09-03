# ****************************************************************************#
# STEP 5.1A: Generate integration comparison UMAPs
# ****************************************************************************#
# TODO: Define color palettes in 3.1 to standardize colors across all plots
#       across all plots to ensure that the same sample/condition is always
#       the same color throughout the analysis pipeline.
# TODO: Consider creating a separate script for loading the Seurat objects
#       from the checkpoint directory because subsequent steps in the
#       need them and a user may not be able to debug the cause of missing
#       data

# ****************************************************************************#
# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Steps 3.3A, 4.2A-D produced five checkpointed objects on disk: the naive,
#   uncorrected baseline plus four independently-run integration methods
#   (CCA, RPCA, Harmony, FastMNN). Each method reported its own cluster count
#   in isolation, but a number alone can't tell you which correction actually
#   worked best — that requires looking at all five side-by-side.
#
#   This step reloads all five checkpoints fresh from disk (rather than
#   relying on in-memory objects from previous steps in this script) and
#   renders two comparison grids:
#     - By SAMPLE: does each method actually mix samples together, or do
#       patient-specific islands persist?
#     - By CONDITION: does biological signal (Healthy vs Post-Treatment)
#       survive the correction, or does integration over-correct and blend
#       away real disease-driven differences along with the batch noise?
#
#   Step 3.3C's per-cluster mixing statistic is the quantitative version of
#   this comparison; these UMAPs are the visual companion piece.


# --- 1. Reload Sample Metadata (Single Source of Truth) ---
# ****************************************************************************#
#   Re-read from the same `sample_names.tsv` used throughout the pipeline
#   rather than trusting condition/patient labels already baked into the
#   checkpointed objects. This script may run as a standalone job well after
#   previous steps finished, so metadata is reloaded explicitly instead of
#   assumed to still be in memory.
sample_metadata <- read.delim("2_input/sample-metadata/sample_names.tsv",
    stringsAsFactors = FALSE
) %>%
    select(sample_id, condition, patient_id) %>%
    filter(condition != "Periodontitis_Pre_Treatment")

cat("Sample metadata loaded:\n")
cat("  Samples:   ", nrow(sample_metadata), "\n")
cat("  Conditions:", paste(unique(sample_metadata$condition), collapse = ", "), "\n\n")


# --- 2. Define Color Palettes ---
# ****************************************************************************#
#   Matches the same two palettes used in Step 3.3B, for visual continuity
#   between the "before" (naive merge) figure and these "after" (integrated)
#   figures.
#     - `sample_colors`: one distinct color per patient/sample (8 total),
#       used to check whether each method actually mixes samples together.
#     - `condition_colors`: one color per biological condition, used to
#       check whether real disease-state separation survives correction.
#
#   NAMING CAVEAT (see Step 3.3B): the names below must exactly match the
#   values in `merged_naive$condition` / `integrated_*$condition` — check
#   with `unique(sample_metadata$condition)`. A mismatched name doesn't
#   error, DimPlot just silently falls back to its default palette for that
#   condition.
sample_colors <- c(
    # Healthy
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    # Post_Patient
    "#FF7F00", "#A65628", "#F781BF", "#999999"Gray
)
names(sample_colors) <- sample_metadata$sample_id

condition_colors <- c(
    "Healthy" = "#2E86AB",                      # Blue
    "Periodontitis_Post_Treatment" = "#F18F01"  # Orange
)


# --- 3. Helper: Standardized Comparison UMAP Builder ---
# ****************************************************************************#
#   Wraps `DimPlot()` with the consistent styling (title, legend size) we
#   want across all ten panels (5 methods × 2 grouping variables), so each
#   individual plot call below only has to specify what differs: which
#   object, which reduction, which grouping variable, which title.
make_comparison_plot <- function(seurat_obj, reduction_name, title,
                                group_by = "sample_id", colors = NULL) {
    p <- DimPlot(seurat_obj,
        reduction = reduction_name, group.by = group_by, pt.size = 0.05
    ) +
        ggtitle(title) +
        theme(
            plot.title = element_text(face = "bold", size = 12),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.3, "cm")
        )

    # Apply custom colors if provided
    if (!is.null(colors)) {
        p <- p + scale_color_manual(values = colors)
    }

    return(p)
}


# --- 4. Load All Five Checkpointed Objects ---
# ****************************************************************************#
#   Loads the naive baseline (Step 3.3A) plus all four integration methods
#   (Steps 4.2A-D) back from the checkpoint vault built in Step 3.1
#   (`DATA_CHECKPOINT_DIR`).
#
#   FORMAT HANDLING: `CHECKPOINT_FORMAT` (set once in Step 3.1: 1 = qs2,
#   2 = rds) governs every checkpoint read/write in this pipeline. This
#   loader respects that same switch rather than hardcoding `qs2::qs_read()`,
#   so a run configured for `rds` doesn't silently fail here.

# Create a named vector of checkpoint names to load
CHECKPOINT_NAMES <- c(
    naive   = "01_merged_naive",
    cca     = "02_integrated_cca",
    rpca    = "03_integrated_rpca",
    harmony = "04_integrated_harmony",
    fastmnn = "05_integrated_fastmnn"
)

# Get the checkpoint extension based on the checkpoint format
CHECKPOINT_EXT <- if (CHECKPOINT_FORMAT == 1) {
    "qs2"
} else if (CHECKPOINT_FORMAT == 2) {
    "rds"
} else {
    stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}

# Helper function to load a checkpoint
load_checkpoint <- function(name) {
    path <- file.path(DATA_CHECKPOINT_DIR, paste0(name, ".", CHECKPOINT_EXT))
    if (CHECKPOINT_FORMAT == 1) {
        qs2::qs_read(path)
    } else {
        readRDS(path)
    }
}

# Load all checkpoints
CHECKPOINTS <- LOG_STEP("Loading all 5 checkpointed objects for comparison...", {
    setNames(
        lapply(CHECKPOINT_NAMES, load_checkpoint),
        names(CHECKPOINT_NAMES)
    )
})

# Extract the Seurat objects from the named vector
merged_naive        <- CHECKPOINTS$naive
integrated_cca      <- CHECKPOINTS$cca
integrated_rpca     <- CHECKPOINTS$rpca
integrated_harmony  <- CHECKPOINTS$harmony
integrated_fastmnn  <- CHECKPOINTS$fastmnn

# Print summary of loaded checkpoints
cat("✓ All 5 checkpoints loaded\n\n")


# --- 5. Build Comparison Panels: Colored by Sample (Mixing Check) ---
# ****************************************************************************#
#   Five panels, one per method, all colored the same way (by sample_id).
#   A method that successfully corrected batch effects should show samples
#   thoroughly interleaved within shared clusters, in visible contrast to
#   the naive panel where patient-specific islands are expected.
p_sample_naive <- make_comparison_plot(
    merged_naive, "umap", "Naive Merge",
    colors = sample_colors
)
p_sample_cca <- make_comparison_plot(
    integrated_cca, "umap.cca", "CCA",
    colors = sample_colors
)
p_sample_rpca <- make_comparison_plot(
    integrated_rpca, "umap.rpca", "RPCA",
    colors = sample_colors
)
p_sample_harmony <- make_comparison_plot(
    integrated_harmony, "umap.harmony", "Harmony",
    colors = sample_colors
)
p_sample_mnn <- make_comparison_plot(
    integrated_fastmnn, "umap.mnn", "FastMNN",
    colors = sample_colors
)

# Assemble into a 2x3 grid (5 plots + 1 empty spacer to fill the grid)
combined_samples <- (p_sample_naive | p_sample_cca | p_sample_rpca) /
    (p_sample_harmony | p_sample_mnn | plot_spacer())

# Save comparison by sample image
ggsave(
    file.path(PLOTS_COMPARISON_DIR, "02_integration_by_sample.png"),
    plot = combined_samples, width = 16, height = 12, dpi = 300
)

# Print output message
cat("→ Saved: ", file.path(PLOTS_COMPARISON_DIR, "02_integration_by_sample.png"), "\n")
cat("   Look for samples interleaving within shared clusters (good) vs\n")
cat("   forming isolated single-sample islands (batch effect persists)\n\n")


# --- 6. Build Comparison Panels: Colored by Condition (Biology Check) ---
# ****************************************************************************#
#   The counterpart check to Step 5: correcting batch effects is only a
#   win if real biological signal (Healthy vs Post-Treatment) is preserved
#   alongside it. A method that blends conditions together as thoroughly as
#   it blends samples has likely over-corrected, erasing genuine disease
#   signal along with the technical noise.
p_cond_naive <- make_comparison_plot(
    merged_naive, "umap", "Naive Merge", "condition",
    colors = condition_colors
)
p_cond_cca <- make_comparison_plot(
    integrated_cca, "umap.cca", "CCA", "condition",
    colors = condition_colors
)
p_cond_rpca <- make_comparison_plot(
    integrated_rpca, "umap.rpca", "RPCA", "condition",
    colors = condition_colors
)
p_cond_harmony <- make_comparison_plot(
    integrated_harmony, "umap.harmony", "Harmony", "condition",
    colors = condition_colors
)
p_cond_mnn <- make_comparison_plot(
    integrated_fastmnn, "umap.mnn", "FastMNN", "condition",
    colors = condition_colors
)

combined_conditions <- (p_cond_naive | p_cond_cca | p_cond_rpca) /
    (p_cond_harmony | p_cond_mnn | plot_spacer())

# Save comparison by condition image
ggsave(
    file.path(PLOTS_COMPARISON_DIR, "03_integration_by_condition.png"),
    plot = combined_conditions, width = 16, height = 12, dpi = 300
)

# Print output message
cat("→ Saved: ", file.path(PLOTS_COMPARISON_DIR, "03_integration_by_condition.png"), "\n")
cat("   Look for condition-driven separation surviving correction (good) vs\n")
cat("   conditions blending together indiscriminately (over-correction)\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Previous Steps left us with five independently checkpointed objects on
#   disk — the naive baseline plus four integration methods — each only
#   evaluated in isolation via its own cluster count.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We reloaded all five checkpoints together and rendered two 2x3
#   comparison grids: one colored by sample (the batch-mixing check) and
#   one colored by condition (the biological-preservation check). These two
#   figures (02_integration_by_sample.png, 03_integration_by_condition.png)
#   are now the direct visual evidence for judging which method, if any,
#   best balances removing technical batch effects against preserving real
#   Healthy-vs-Post-Treatment biological signal.
#
# WHERE WE ARE HEADING (STEP 5.1B — QUANTIFYING INTEGRATION SUCCESS):
#   A visual grid is persuasive but, exactly as Step 3.3B cautioned for the
#   naive baseline, not something you can report as a number. Next we
#   re-apply Step 3.3C's per-cluster sample-composition statistic to all
#   four integrated objects, checking whether n_batch_clusters actually
#   drops relative to the naive merge's 1/19 — turning this visual
#   impression into the same kind of citable, quantitative evidence used
#   to justify integration in the first place.
# ****************************************************************************#