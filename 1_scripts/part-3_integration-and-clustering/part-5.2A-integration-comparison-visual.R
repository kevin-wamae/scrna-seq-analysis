# ****************************************************************************#
# STEP 5.2A: Generate integration comparison UMAPs
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP,
#     output-directory variables, and the shared `sample_colors` /
#     `condition_colors` palettes.
#   • part-5.1-load-seurat-object-checkpoints.R — supplies the five objects
#     plotted here: `merged_naive`, `integrated_cca`, `integrated_rpca`,
#     `integrated_harmony`, `integrated_fastmnn`.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source the prerequisites (interactive) or stops with a clear
#   message (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-5.2A-integration-comparison-visual.R")

# NOTE: Requires Step 5.1 (load-seurat-object-checkpoints) to have already
#       run in this session — that step is what populates `merged_naive`,
#       `integrated_cca`, `integrated_rpca`, `integrated_harmony`, and
#       `integrated_fastmnn`.

# ****************************************************************************#
# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Steps 3.3A, 4.2A-D produced five checkpointed objects on disk: the naive,
#   uncorrected baseline plus four independently-run integration methods
#   (CCA, RPCA, Harmony, FastMNN). Each method reported its own cluster count
#   in isolation, but a number alone can't tell you which correction actually
#   worked best — that requires looking at all five side-by-side.
#
#   This step takes the five checkpointed objects loaded by Step 5.1 and
#   renders two comparison grids:
#     - By SAMPLE: does each method actually mix samples together, or do
#       patient-specific islands persist?
#     - By CONDITION: does biological signal (Healthy vs Post-Treatment)
#       survive the correction, or does integration over-correct and blend
#       away real disease-driven differences along with the batch noise?
#
#   Step 3.3C's per-cluster mixing statistic is the quantitative version of
#   this comparison; these UMAPs are the visual companion piece.


# --- 1. Shared Metadata & Color Palettes ---
# ****************************************************************************#
#   `sample_metadata`, `sample_colors`, and `condition_colors` are defined
#   once in Step 3.1 and reused here. Because this step can run as a
#   standalone job well after 3.1, they arrive via `ensure_dependencies()`
#   like every other shared object — nothing is re-read from the TSV or
#   re-defined locally, so this comparison colours each sample and condition
#   exactly as the naive "before" figure (Step 3.3B) and every later figure.


# --- 2. Helper: Standardized Comparison UMAP Builder ---
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


# --- 3. Build Comparison Panels: Colored by Sample (Mixing Check) ---
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


# --- 4. Build Comparison Panels: Colored by Condition (Biology Check) ---
# ****************************************************************************#
#   The counterpart check to Section 4 above: correcting batch effects is
#   only a win if real biological signal (Healthy vs Post-Treatment) is
#   preserved alongside it. A method that blends conditions together as
#   thoroughly as it blends samples has likely over-corrected, erasing
#   genuine disease signal along with the technical noise.
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
#   Step 5.1 left us with five checkpointed objects loaded into memory — the
#   naive baseline plus four integration methods — each only evaluated in
#   isolation via its own cluster count.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We rendered two 2x3 comparison grids from those five objects: one
#   colored by sample (the batch-mixing check) and one colored by condition
#   (the biological-preservation check). These two figures
#   (02_integration_by_sample.png, 03_integration_by_condition.png) are now
#   the direct visual evidence for judging which method, if any, best
#   balances removing technical batch effects against preserving real
#   Healthy-vs-Post-Treatment biological signal.
#
# WHERE WE ARE HEADING (STEP 5.2B — QUANTIFYING INTEGRATION SUCCESS):
#   A visual grid is persuasive but, exactly as Step 3.3B cautioned for the
#   naive baseline, not something you can report as a number. Next we
#   re-apply Step 3.3C's per-cluster sample-composition statistic to all
#   four integrated objects, checking whether n_batch_clusters actually
#   drops relative to the naive merge's 1/19 — turning this visual
#   impression into the same kind of citable, quantitative evidence used
#   to justify integration in the first place.
# ****************************************************************************#
