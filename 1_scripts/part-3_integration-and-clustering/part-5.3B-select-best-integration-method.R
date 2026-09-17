# ****************************************************************************#
# STEP 5.3: Select optimal integration method
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and
#     output-directory variables.
#   • part-5.1-load-seurat-object-checkpoints.R — supplies the five objects.
#   • part-5.2B-integration-comparison-quantitatively.R — supplies
#     `methods_list` and `mixing_results`.
#   • part-5.3A-assess-biological-preservation.R — adds the
#     `condition_separation` column to `mixing_results`, which the ranking
#     below prints alongside each method's mixing score.
#   If `mixing_results` (with `condition_separation`) is already in the
#   environment, this block does nothing; otherwise it offers to source the
#   prerequisites (interactive) or stops with a clear message
#   (non-interactive/batch).
if (!exists("ensure_dependencies", inherits = TRUE)) {
    source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
}
ensure_dependencies(step = "part-5.3B-select-best-integration-method.R")

# NOTE: Requires Steps 5.1, 5.2B, and 5.3A to have already run in this
#       session — this step ranks the `mixing_results` table (mixing score +
#       condition separation) built across those steps, and selects an
#       object out of `methods_list`.

# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Steps 5.2B and 5.3A gave us two numbers per method — mixing score and
#   condition separation — plus a trade-off plot to eyeball them together.
#   This step turns that comparison into an actual decision: rank the five
#   methods, print the ranking for the record, and carry the winning
#   method's object and UMAP reduction forward as `integrated_final` /
#   `reduction_final` for all downstream clustering and differential
#   expression work.
#
#   PRIMARY CRITERION: sample mixing score. Good integration requires
#   samples to mix within shared clusters; a method that fails to do this
#   hasn't corrected the batch effect regardless of how well it preserved
#   biology. Condition separation is reported alongside the ranking as a
#   sanity check — see Step 5.3A — so an obviously over-integrated top
#   pick (high mixing, near-zero separation) is visible before committing
#   to it, even though it isn't currently used to re-rank.
#
#   INTERPRETATION GUIDE:
#     - High mixing score (>3.0): good integration, samples mix well
#       within cell types
#     - Moderate separation: biological conditions remain distinguishable
#     - Balanced scores: ideal integration removes technical variation
#       while preserving biology
#
#   WARNING SIGNS TO WATCH FOR:
#     - Very high mixing + very low separation: over-integration
#       (lost biological signal)
#     - Very low mixing + very high separation: under-integration
#       (batches not corrected)


# --- 1. Rank Methods by Mixing Score ---
# ****************************************************************************#
#   `dplyr::arrange(desc(...))` replaces the base `order()` re-indexing —
#   same result, easier to read as "sort mixing_results by mixing_score,
#   descending" rather than reasoning through a negated column index.
mixing_results <- mixing_results %>%
    arrange(desc(mixing_score))

# Print the ranking, most-mixed method first
cat("\nMethod ranking (by sample mixing):\n")
mixing_results %>%
    mutate(rank = row_number()) %>%
    pwalk(function(rank, method, mixing_score, condition_separation, ...) {
        cat(sprintf(
            "  %d. %s (mixing: %.2f, separation: %.2f)\n",
            rank, method, mixing_score, condition_separation
        ))
    })


# --- 2. Select the Top-Ranked Method ---
# ****************************************************************************#
best_method <- mixing_results$method[1]
cat("\nRecommended method:", best_method, "\n")


# --- 3. Retrieve the Winning Object and Reduction ---
# ****************************************************************************#
#   `methods_list` (built in Step 5.2B) already maps each method name to
#   its Seurat object and UMAP reduction, so the winning pair can be looked
#   up directly instead of re-declaring the same five method-to-object
#   mappings a second time in an if/else chain. This also means adding a
#   sixth integration method later only requires updating `methods_list`
#   once, in 5.2B — not here as well.
if (!best_method %in% names(methods_list)) {
    # Should be unreachable, since best_method is drawn from mixing_results,
    # which is itself derived from methods_list — kept as a safety net in
    # case that assumption is ever broken upstream.
    stop("Selected method '", best_method, "' not found in methods_list.")
}

integrated_final <- methods_list[[best_method]]$obj
reduction_final <- methods_list[[best_method]]$reduction

cat("\nUsing", best_method, "for downstream clustering analysis\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Five methods, each with a mixing score and a condition-separation score
#   sitting side by side in `mixing_results`, but no single object carried
#   forward for downstream analysis.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We ranked all five methods by mixing score, printed that ranking
#   alongside each method's separation score for the record, and selected
#   `best_method`'s object and UMAP reduction into `integrated_final` /
#   `reduction_final` — the single integration result the rest of the
#   pipeline will build on.
#
# WHERE WE ARE HEADING (DOWNSTREAM CLUSTERING):
#   With one integration method now committed to, the next step re-runs
#   clustering on `integrated_final` using `reduction_final`, and downstream
#   differential expression and cell-type annotation proceed from there.
# ****************************************************************************#