# ****************************************************************************#
# STEP 5.3A: Assess biological signal preservation
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and
#     output-directory variables.
#   • part-5.1-load-seurat-object-checkpoints.R — supplies the five objects
#     scored here (via `methods_list`).
#   • part-5.2B-integration-comparison-quantitatively.R — supplies
#     `methods_list` and `mixing_results`, which this step extends with a
#     `condition_separation` column.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if `mixing_results` is already in the environment, otherwise it
#   offers to source the prerequisites (interactive) or stops with a clear
#   message (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-5.3A-assess-biological-preservation.R")

# NOTE: Requires Step 5.1 (load-seurat-object-checkpoints) and Step 5.2B
#       (integration-comparison-quantitatively) to have already run in this
#       session — this step extends `mixing_results` and reuses
#       `methods_list`, both built in 5.2B.

# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   A good integration method should remove technical variation (samples
#   mix) while preserving biological variation (conditions remain
#   distinguishable). Step 5.2B's mixing score only checks the first half
#   of that — it says nothing about whether real biology survived.
#
#   HOW THIS METRIC WORKS (FOR BEGINNERS):
#   Integration is a balancing act. We want to remove batch effects
#   (technical variation) but keep biological differences (biological
#   variation). This metric helps us check if we've struck the right
#   balance.
#
#   THE MEASUREMENT PROCESS:
#     1. Find condition centers: calculate the average UMAP position for
#        all Healthy cells and all Post_Treatment cells.
#     2. Measure separation: calculate the distance between these two
#        centers.
#     3. Interpret the distance:
#          - High distance = conditions are well-separated (biological
#            signal preserved) ✓
#          - Low distance  = conditions overlap (might have
#            over-integrated, losing biology)
#
#   WHAT WE'RE LOOKING FOR:
#     - High mixing score (Step 5.2B)  → samples within each condition mix well
#     - Moderate-to-high separation (this step) → conditions remain distinguishable
#
#   EXAMPLES:
#     - Good integration:  mixing high, separation high  (samples blend,
#       but Healthy vs Post_Treatment stay in distinct regions)
#     - Over-integration:  mixing high, separation low   (everything
#       mixes together, biological differences lost)
#     - Under-integration: mixing low,  separation high  (samples stay
#       separated, batch effects remain uncorrected)


# --- 1. Helper: Condition Separation for One Seurat Object ---
# ****************************************************************************#
#   Computes each condition's mean position ("center") in 2D UMAP space,
#   then returns the mean pairwise distance between those centers. With
#   two conditions (Healthy vs Post_Treatment) this is just a single
#   center-to-center distance; written generally so a third condition
#   wouldn't silently break it.
#
#   EFFICIENCY NOTE: pulls the embedding matrix and condition vector once
#   up front (avoids repeated `seurat_obj[[...]]` slot access), and calls
#   `stats::dist()` explicitly — this pipeline also loads packages (e.g.
#   `spam`, `BiocGenerics`) that register an unrelated S4 class of the
#   same name, so the namespace is spelled out to keep it unambiguous.
assess_condition_separation <- function(seurat_obj, reduction) {
    # Pull embedding and condition labels once, up front
    umap_coords <- Embeddings(seurat_obj, reduction = reduction)[, 1:2]
    conditions <- seurat_obj$condition

    # Step 1: Calculate mean UMAP position for each condition
    condition_centers <- aggregate(umap_coords, by = list(conditions), FUN = mean)

    # Step 2: Calculate distance between condition centers
    # (with only 2 conditions this is a single number; dist() generalizes
    # to the mean pairwise distance if more conditions are ever added)
    condition_centers[, -1] %>%
        stats::dist() %>%
        mean()
}


# --- 2. Calculate Separation for All Methods ---
# ****************************************************************************#
#   `map_dbl()` iterates the helper over the same `methods_list` built in
#   Step 5.2B, returning one separation score per method in a single pass —
#   the direct counterpart to that step's mixing-score calculation.
separation_results <- methods_list %>%
    map_dbl(~ assess_condition_separation(.x$obj, .x$reduction))

# Attach to the existing mixing_results table, matched explicitly by method
# name rather than assumed row order
mixing_results$condition_separation <- separation_results[mixing_results$method]

# Display results
cat("\nBiological Signal Preservation:\n")
cat("(Higher separation = conditions remain distinguishable)\n\n")
print(mixing_results[, c("method", "mixing_score", "condition_separation", "n_clusters")])




# --- 3. Visualize the Mixing vs. Separation Trade-off ---
# ****************************************************************************#
#   Ideal integration methods land in the upper right: high mixing (samples
#   blend within conditions) AND high separation (conditions remain
#   distinguishable). `geom_text_repel` keeps the five method labels legible
#   even when points sit close together.
p_tradeoff <- mixing_results %>%
    ggplot(aes(x = mixing_score, y = condition_separation, label = method)) +
    geom_point(size = 4, aes(color = method)) +
    geom_text_repel(size = 4, box.padding = 0.5) +
    labs(
        title = "Integration Quality: Mixing vs Biological Preservation",
        subtitle = "Ideal: Upper right (high mixing + preserved biology)",
        x = "Sample Mixing Score (higher = better integration)",
        y = "Condition Separation (higher = preserved biology)"
    ) +
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set2")

# Save mixing vs preservation image
ggsave(file.path(PLOTS_COMPARISON_DIR, "05_mixing_vs_preservation.png"),
    p_tradeoff, width = 10, height = 7, dpi = 300
)


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 5.2B gave us a mixing score per method but no way to tell whether
#   high mixing reflected genuine batch correction or an over-integrated
#   method that blended away real biological signal along with it.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We computed a condition-separation score for all five methods — the
#   distance between each condition's mean UMAP position — and combined it
#   with the existing mixing scores into one table and one trade-off plot
#   (05_mixing_vs_preservation.png). Methods can now be judged on both axes
#   at once: which one, if any, sits in the "high mixing + preserved
#   biology" upper-right region rather than trading one for the other.
#
# WHERE WE ARE HEADING (SELECTING THE BEST INTEGRATION METHOD):
#   With mixing score, condition separation, and cluster count now sitting
#   side by side in `mixing_results` for all five methods, the next step
#   turns this comparison table into a decision: weighing the trade-offs
#   captured in 05_mixing_vs_preservation.png to pick the one integration
#   method that will carry forward into downstream clustering and
#   differential expression analysis.
# ****************************************************************************#