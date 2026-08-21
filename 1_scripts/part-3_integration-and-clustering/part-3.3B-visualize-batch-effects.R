# ****************************************************************************#
# STEP 4: Visualize batch effects in naive merge
# ****************************************************************************#


# --- THE PURPOSE OF THIS DIAGNOSTIC ---
# ****************************************************************************#
#   Step 3 gave us cluster counts and cell totals, but numbers alone can't
#   show you WHERE batch effects live in the data. A UMAP colored by sample
#   makes technical separation immediately visible: if patients form their
#   own isolated islands rather than mixing within shared biological cell-type
#   clusters, that's the batch effect Harmony/FastMNN are meant to correct in
#   Step 5. This plot is the "before" picture we'll compare every integration
#   method against.


# --- 1. Define Color Palettes ---
# ****************************************************************************#
#   Two separate palettes, matched to two separate questions:
#     - `sample_colors`: one distinct color per patient/sample. Lets you see
#       whether any single sample clusters apart from the rest (a batch
#       effect) versus blending into shared clusters (no batch effect).
#     - `condition_colors`: one color per biological condition (Healthy vs
#       Post-Treatment). Lets you see whether disease state drives real
#       separation, independent of which specific patient a cell came from.
#
#   8 samples need 8 distinct, visible colors (Set1-style categorical palette,
#   manually assigned rather than auto-generated so colors stay stable across
#   re-runs regardless of factor level ordering).
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", # Healthy 1-4: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#999999"  # Post_Patient 1-4: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# NOTE: keys here must exactly match the values in `merged_naive$condition`
# (which come from the TSV's `condition` column — currently
# "Periodontitis_Post_Treatment", not "Post_Treatment"). A mismatched key
# means `DimPlot` silently falls back to its default palette for that
# condition instead of erroring, so double-check this against
# `unique(merged_naive$condition)` before trusting the plot's colors.
condition_colors <- c(
  "Healthy" = "#2E86AB",          # Blue
  "Post_Treatment" = "#F18F01"    # Orange
)


# --- 2. Build Diagnostic UMAP Panels ---
# ****************************************************************************#
#   Four complementary views of the same embedding, each answering a
#   different question:
#     - p1: Are samples separating out on their own? (batch effect signature)
#     - p2: Is biological condition driving separation instead? (the signal
#       we actually want to see, once batch effects are corrected)
#     - p3: What do the naive (uncorrected) cluster boundaries look like?
#     - p4: Faceted per-sample — makes it easy to spot ONE specific sample
#       behaving differently from the rest, which an overlaid plot can hide.

# UMAP colored by sample - shows batch effects
p1_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "sample_id",
                     pt.size = 0.05, cols = sample_colors) +
  ggtitle("Naive Merge: Colored by Sample") +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# UMAP colored by condition
p2_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                     pt.size = 0.05, cols = condition_colors) +
  ggtitle("Naive Merge: Colored by Condition")

# UMAP colored by clusters
p3_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "seurat_clusters",
                     pt.size = 0.05, label = TRUE, label.size = 5) +
  ggtitle("Naive Merge: Clusters") +
  NoLegend()

# Split by sample to see separation
p4_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                     split.by = "sample_id", pt.size = 0.05, ncol = 4,
                     cols = condition_colors) +
  ggtitle("Naive Merge: Split by Sample") +
  theme(strip.text = element_text(size = 9, face = "bold"))


# --- 3. Assemble & Save Diagnostic Panel ---
# ****************************************************************************#
#   `patchwork`'s `|` (side-by-side) and `/` (stacked) operators combine all
#   four panels into one 2x2 figure, so sample-level, condition-level, and
#   cluster-level views can be compared at a glance in a single saved file.
combined_naive <- (p1_naive | p2_naive) / (p3_naive | p4_naive)

ggsave(
  file.path(PLOTS_OUT_DIR, "01_naive_merge_batch_effects.png"),
  plot = combined_naive, width = 16, height = 12, dpi = 300
)

cat("→ EXAMINE file:", file.path(PLOTS_OUT_DIR, "01_naive_merge_batch_effects.png"), "\n")
cat("   Look for samples forming isolated islands (batch effect) vs\n")
cat("   blending into shared clusters (no strong batch effect)\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 3.3A gave us a merged-but-uncorrected object and a bare number: 19
#   clusters. That number alone can't distinguish real biology from batch
#   effects — a cluster count says nothing about WHY cells grouped the way
#   they did.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We rendered the naive merge as a four-panel diagnostic UMAP — colored by
#   sample, by condition, by cluster identity, and split per-sample — and
#   saved it to 01_naive_merge_batch_effects.png. This turns "batch effects
#   might be present" into something directly visible: wherever a specific
#   sample or patient forms an isolated island rather than blending into a
#   shared cluster, that's technical variation masquerading as biological
#   signal. This figure is now our fixed "before" reference — every
#   integration method's UMAP in later steps gets compared back against it.
#
# WHERE WE ARE HEADING (STEP 5):
#   A visual impression of batch effects is persuasive but not quantitative
#   — "that cluster looks mostly red" isn't a number you can report or set a
#   threshold against. In Step 5, we will calculate the actual per-cluster
#   sample composition and flag any cluster dominated by a single sample,
#   turning what Step 4 showed us into a concrete, citable statistic.
# ****************************************************************************#