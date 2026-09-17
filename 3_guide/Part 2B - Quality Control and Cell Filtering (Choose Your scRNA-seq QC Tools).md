---
title: "How to Choose Your scRNA-seq QC Tools (Part2-2): SoupX vs DecontX and DoubletFinder vs scDblFinder"
source: "https://ngs101.com/how-to-choose-your-scrna-seq-qc-tools-part2-2-soupx-vs-decontx-and-doubletfinder-vs-scdblfinder/"
author:
  - "[[Lei]]"
published: 2026-06-27
created: 2026-07-12
description: "A decision-focused quick tip for picking the right ambient RNA correction and doublet detection tools — run side by side on real 10x PBMC data Introduction: Same Job, Different Tools In Part 2 of this series, we built a full quality control (QC) workflow that corrected ambient RNA contamination and removed doublets before any downstream […]"
tags:
  - "clippings"
---
*A decision-focused quick tip for picking the right ambient RNA correction and doublet detection tools — run side by side on real 10x PBMC data*

## Introduction: Same Job, Different Tools

In [Part 2 of this series](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/), we built a full quality control (QC) workflow that corrected ambient RNA contamination and removed doublets before any downstream analysis. We picked one tool for each job and moved on. But if you have searched for scRNA-seq QC advice, you have certainly run into the two questions this quick tip answers:

- For ambient RNA correction, should I use **SoupX** or **DecontX**?
- For doublet detection, should I use **DoubletFinder** or **scDblFinder**?

These are the two most widely used tools in each category, and beginners constantly ask which one to choose. The honest answer is that they solve the same problem with different assumptions, different inputs, and different practical trade-offs. This tutorial runs all four on the same sample so you can see how they compare and — more importantly — gives you a clear decision rule for your own data.

This is a quick tip, not a full pipeline tutorial. We assume you already understand what ambient RNA and doublets are. If you need that background, read [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/) first.

### What You Will Learn

By the end of this quick tip, you will be able to:

- ✓ Run SoupX and DecontX on the same sample and compare their contamination estimates
- ✓ Run DoubletFinder and scDblFinder on the same sample and compare their doublet calls
- ✓ Understand the conceptual difference between each pair of tools
- ✓ Apply a clear decision rule for which tool to use on your own dataset
- ✓ Combine tools sensibly (for example, taking consensus doublet calls)

> **Quick orientation:** Ambient RNA correction and doublet detection are two separate QC steps. SoupX and DecontX both address ambient RNA. DoubletFinder and scDblFinder both address doublets. You do not choose between an ambient tool and a doublet tool — you run one of each.

## The Two Problems in One Minute

Before comparing tools, a one-line refresher on what each pair is correcting.

**Ambient RNA contamination.** When cells are dissociated, some lyse and release mRNA into the suspension. This free-floating “soup” gets co-captured in every droplet, so each cell’s count matrix contains a fraction of molecules that did not originate from that cell. Left uncorrected, ambient RNA inflates marker genes everywhere and creates false positive expression. SoupX and DecontX both estimate the contamination fraction and subtract the soup.

**Doublets.** When two cells are captured in the same droplet, they are sequenced under one barcode and appear as a single “cell” whose transcriptome is a blend of two real cells. Doublets create artificial intermediate cell types and inflate your cell count. DoubletFinder and scDblFinder both flag these droplets so you can remove them.

The key insight that motivates this comparison: each pair attacks the same problem from a different statistical angle, so they will not produce identical results. Understanding where they agree and disagree tells you how much to trust the output.

## Setting Up Your Environment

We need both pairs of tools plus the standard single-cell stack. SoupX and DoubletFinder are Seurat-friendly; DecontX and scDblFinder are Bioconductor tools that operate on the `SingleCellExperiment` (SCE) object covered in [Part 6](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-6-understanding-seurat-and-singlecellexperiment-objects/).

```bash
#-----------------------------------------------
# STEP 0: Install the four QC tools (run once)
#-----------------------------------------------
 
options(repos = c(CRAN = "https://cloud.r-project.org"))
 
# CRAN / GitHub tools (Seurat ecosystem)
install.packages("SoupX")
 
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
# DoubletFinder has been Seurat 5 compatible since Nov 2023; the '_v3'
# function suffix was removed in that update. Install the current version.
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
 
# Bioconductor tools (SingleCellExperiment ecosystem)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# DecontX now ships as its own package (split out from celda)
BiocManager::install(c("decontX", "scDblFinder", "SingleCellExperiment"))
```


```bash
#-----------------------------------------------
# STEP 1: Load libraries and configure the session
#-----------------------------------------------
 
library(Seurat)
library(SoupX)
library(decontX)
library(scDblFinder)
library(DoubletFinder)
library(SingleCellExperiment)
 
library(ggplot2)
library(dplyr)
library(patchwork)
 
set.seed(12345)
theme_set(theme_classic())
 
dir.create("qc_compare/plots", recursive = TRUE, showWarnings = FALSE)
```

> **Object format note:** SoupX returns a corrected sparse matrix you can drop straight back into a Seurat object. DecontX and scDblFinder expect a `SingleCellExperiment`. DoubletFinder annotates a Seurat object in place. We convert between the two formats as needed, so keep track of which object each tool writes to.

## Loading the Data

We use the same GSE174609 PBMC periodontitis dataset as the rest of the series, processed through Cell Ranger in [Part 1](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-1-from-fastq-to-count-matrix/). To keep it simple, we demonstrate on a single representative sample (`Healthy_1`). The same code applies to any sample — wrap it in `lapply` over your sample list to process all twelve.

Ambient RNA tools need different inputs, so we load both matrices up front:

- The **filtered** matrix (`filtered_feature_bc_matrix`) contains called cells only. All four tools use this.
- The **raw** matrix (`raw_feature_bc_matrix`) contains every droplet, including empty ones. SoupX uses this to profile the soup; DecontX does not need it.
```bash
#-----------------------------------------------
# STEP 2: Load Cell Ranger output for one sample
#-----------------------------------------------
 
# Adjust this path to your Part 1 Cell Ranger output
cellranger_path <- "~/GSE174609_scRNA/cellranger_output/Healthy_1/outs"
 
raw_counts  <- Read10X(file.path(cellranger_path, "raw_feature_bc_matrix"))
filt_counts <- Read10X(file.path(cellranger_path, "filtered_feature_bc_matrix"))
 
# Build a Seurat object from the called cells
seurat_obj <- CreateSeuratObject(
  counts  = filt_counts,
  project = "Healthy_1"
)
```

### A Shared Clustering for Fair Comparison

Every one of these tools benefits from cell cluster labels: SoupX uses clusters to estimate the soup, DecontX uses them to model contamination, and both doublet tools use them to generate realistic artificial doublets. To compare the tools fairly, we compute clusters **once** and feed the same labels into all four.

```R
#-----------------------------------------------
# STEP 3: One shared Seurat preprocessing + clustering
#-----------------------------------------------
 
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
 
# Shared cluster labels, named by barcode so every tool aligns the same cells
shared_clusters <- setNames(
  as.character(seurat_obj$seurat_clusters),
  colnames(seurat_obj)
)
```

> **Why share the clustering:** If each tool computed its own clusters with its own defaults, differences in output could come from the clustering rather than the decontamination or doublet model. Fixing the clustering isolates the part we actually want to compare.

## Round 1: Ambient RNA - SoupX vs DecontX

### How Each Tool Thinks About the Soup

**SoupX** (Young and Behjati, 2020) builds an explicit profile of the ambient RNA by looking at the empty droplets in the raw matrix — droplets that contain only soup, no cell. With the default `autoEstCont`, it then estimates a **single global contamination fraction** (`rho`) for the whole sample and applies it to every cell when it subtracts the soup. Because it reads the soup profile directly from empty droplets, SoupX needs the raw matrix.

**DecontX** (Yang et al., 2020) takes a different route. It models each cell’s expression as a Bayesian mixture of two distributions: counts that belong to the cell’s own population, and counts that look like they leaked in from other populations. It estimates contamination **per cell** from the cell-by-cell data alone, so it does **not** require the raw matrix or empty droplets at all. (It can optionally take a `background` matrix of empty droplets to refine the estimate, but that is not required.) This is its single biggest practical advantage: if you only have a filtered matrix — common when you inherit a public dataset — SoupX is awkward, but DecontX just works.

> **The key conceptual difference — and it shows up in the results below:** SoupX gives you one contamination number for the whole sample; DecontX gives you a contamination number for every cell. A global estimate is simpler and more conservative; a per-cell estimate can adapt to populations that are more or less contaminated, at the cost of being more aggressive. Neither is universally correct — but knowing which granularity each tool works at explains everything we see when we compare them.

### Running SoupX

```R
#-----------------------------------------------
# STEP 4: Ambient RNA correction with SoupX
#-----------------------------------------------
 
# tod = table of droplets (raw); toc = table of counts (filtered cells)
soup_channel <- SoupChannel(tod = raw_counts, toc = filt_counts)
 
# Give SoupX the shared clusters (aligned to the filtered cells)
soup_channel <- setClusters(
  soup_channel,
  shared_clusters[colnames(filt_counts)]
)
 
# Automatically estimate the contamination fraction, then subtract the soup
soup_channel  <- autoEstCont(soup_channel)
soupx_counts  <- adjustCounts(soup_channel, roundToInt = TRUE)
 
# autoEstCont estimates ONE global contamination fraction for the sample and
# stores that same value for every cell, so we just take the single number.
soupx_rho_global <- soup_channel$metaData$rho[1]
```

### Running DecontX

```R
#-----------------------------------------------
# STEP 5: Ambient RNA correction with DecontX
#-----------------------------------------------
 
# DecontX works on a SingleCellExperiment built from the filtered counts only
sce_ambient <- SingleCellExperiment(assays = list(counts = filt_counts))
 
# Pass the shared clusters via 'z' so the comparison is apples-to-apples.
# Note: no raw matrix is supplied -- DecontX does not need empty droplets.
sce_ambient <- decontX(
  sce_ambient,
  z = shared_clusters[colnames(filt_counts)]
)
 
# Per-cell contamination fraction estimated by DecontX
decontx_contam <- sce_ambient$decontX_contamination
names(decontx_contam) <- colnames(sce_ambient)
 
# Decontaminated counts (fractional; round if a downstream tool needs integers)
decontx_counts <- decontXcounts(sce_ambient)
```

### Comparing the Two Estimates

Here is the trap we have to avoid. It is tempting to scatter SoupX against DecontX cell by cell and compute a correlation — but that comparison is meaningless, because SoupX produces a single global number while DecontX produces one value per cell. Plotting a constant against a variable gives a vertical line and a correlation of `NA`. Instead, we show DecontX’s per-cell distribution and mark SoupX’s single global estimate as a reference line. This is the honest comparison: it shows directly how a global estimate relates to the spread of per-cell estimates.

```R
#-----------------------------------------------
# STEP 6: Compare the contamination estimates honestly
#-----------------------------------------------
 
# DecontX is per-cell; SoupX is one global number. Show the DecontX
# distribution with the SoupX global estimate as a dashed reference line.
contam_df <- data.frame(
  barcode = colnames(sce_ambient),
  DecontX = decontx_contam,
  cluster = shared_clusters[colnames(sce_ambient)]
)
 
p_contam_dist <- ggplot(contam_df, aes(x = DecontX)) +
  geom_histogram(bins = 60, fill = "#2A9D8F", alpha = 0.85) +
  geom_vline(xintercept = soupx_rho_global,
             linetype = "dashed", color = "#E76F51", linewidth = 1) +
  annotate("text", x = soupx_rho_global, y = Inf,
           label = paste0("SoupX global rho = ", round(soupx_rho_global, 3)),
           hjust = -0.05, vjust = 2, color = "#E76F51") +
  labs(
    title    = "Contamination estimates: per-cell (DecontX) vs global (SoupX)",
    subtitle = "DecontX assigns each cell its own fraction; SoupX assigns one value to all cells",
    x        = "DecontX per-cell contamination fraction",
    y        = "Number of cells"
  )
 
ggsave("qc_compare/plots/01_ambient_contamination_distribution.png",
       p_contam_dist, width = 8, height = 6, dpi = 300)
```


![](https://i0.wp.com/ngs101.com/wp-content/uploads/2026/06/01_ambient_contamination_distribution.png?resize=2048%2C1536&ssl=1)

> **What the real data shows:** On the `Healthy_1` sample, SoupX estimated a global contamination fraction of about 2.3 percent — a low, clean value typical of a well-prepared PBMC sample. DecontX placed most cells near that same low level but assigned a long tail of cells substantially higher contamination, all the way up toward 1.0. That tail is the whole story: SoupX, by design, cannot flag individual heavily contaminated cells because it only has one number to work with, while DecontX can. Whether that tail represents real per-cell contamination or DecontX being over-eager on certain populations is exactly the judgement call the next plot helps you make.

To make the correction concrete, we pull the single most abundant gene in the SoupX soup profile and show its mean expression per cluster before correction, after SoupX, and after DecontX:

```R
#-----------------------------------------------
# STEP 7: Show the correction on the top soup gene
#-----------------------------------------------
 
# Identify the top contributor to the ambient profile
top_soup_gene <- rownames(soup_channel$soupProfile)[
  which.max(soup_channel$soupProfile$est)
]
 
# Mean expression of that gene per cluster, before vs after each correction
gene_compare <- data.frame(
  cluster   = shared_clusters[colnames(filt_counts)],
  raw       = filt_counts[top_soup_gene, ],
  soupx     = soupx_counts[top_soup_gene, colnames(filt_counts)],
  decontx   = as.numeric(decontx_counts[top_soup_gene, colnames(filt_counts)])
) %>%
  group_by(cluster) %>%
  summarise(
    Raw     = mean(raw),
    SoupX   = mean(soupx),
    DecontX = mean(decontx),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(-cluster, names_to = "stage", values_to = "mean_expr")
 
p_gene_compare <- ggplot(gene_compare,
                         aes(x = cluster, y = mean_expr, fill = stage)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = paste0("Correction of top soup gene: ", top_soup_gene),
    x     = "Cluster",
    y     = "Mean counts per cell",
    fill  = "Stage"
  )
 
ggsave("qc_compare/plots/02_top_soup_gene_correction.png",
       p_gene_compare, width = 10, height = 6, dpi = 300)
```


![](https://i0.wp.com/ngs101.com/wp-content/uploads/2026/06/02_top_soup_gene_correction-scaled.png?resize=2048%2C1229&ssl=1)

> **What the real data shows — and why MALAT1 is a revealing test case:** On `Healthy_1`, the top soup gene was MALAT1, a nuclear long non-coding RNA that is one of the most abundant transcripts in almost every scRNA-seq dataset. The two tools treat it very differently. SoupX, holding to its ~2.3 percent global estimate, shaves MALAT1 down only slightly and uniformly across every cluster — the SoupX bars sit just below the raw bars everywhere. DecontX is far more aggressive and cluster-specific: it cuts MALAT1 by roughly a third in some clusters (for example cluster 20 drops from about 300 to about 100, and cluster 15 from about 475 to about 280), while nearly zeroing it out in a couple of small clusters (18 and 22). In other clusters (3, 9, 13, 19) DecontX leaves it almost untouched, matching SoupX.
> 
> Here is the catch, and it is the reason this gene is worth dwelling on: MALAT1 tops the soup profile because it is ubiquitously and enormously expressed, not necessarily because it is contamination. It is genuinely present in nearly all cells. So DecontX driving it toward zero in some clusters may be a correct read of disproportionate ambient signal — or it may be over-correction of a gene the cells really do express. There is no way to know from this plot alone. ==The practical lesson: when your top soup gene is a ubiquitous high-abundance transcript, treat aggressive per-cluster removal with caution and confirm against genes you know to be lineage-restricted before trusting the correction.==

## Round 2: Doublets - DoubletFinder vs scDblFinder

### How Each Tool Finds Doublets

Both tools share the same core idea: simulate artificial doublets by averaging pairs of real cells, then find real cells that look like those artificial doublets. The difference is in the machinery and the user burden.

**DoubletFinder** (McGinnis et al., 2019) requires you to tune a neighbourhood size parameter, `pK`, through a parameter sweep, and to supply an expected doublet count, `nExp`. It runs on a fully preprocessed Seurat object, one sample at a time. It is powerful but hands-on, and the sweep is slow.

==That `nExp` value comes from 10x Genomics’ published multiplet-rate table, which scales roughly linearly at about 0.8% per 1,000 cells recovered for the standard Next GEM 3′ chemistry (v2/v3/v3.1). This is an assumption set by how many cells were loaded in the lab, not something measured from the count matrix — which is why DoubletFinder makes you supply it. Note the rate is chemistry-specific: the newer GEM-X and HT kits roughly halve it to about 0.4% per 1,000 cells, so using 0.8% on those would over-remove cells.== We further discount `nExp` by the estimated fraction of homotypic doublets (two cells of the same type), because those are statistically invisible to this approach — they look just like a normal cell of that type.

**scDblFinder** (Germain et al., 2021) trains a gradient-boosted classifier on the artificial doublets and requires essentially no tuning. It estimates the expected doublet rate itself, runs in a fraction of the time, handles multiple samples through a single `samples` argument, and in independent benchmarks (Xi and Li, 2021) was among the top-performing methods. It operates on a `SingleCellExperiment`.

### Running DoubletFinder

```R
#-----------------------------------------------
# STEP 8: Doublet detection with DoubletFinder
#-----------------------------------------------
 
# DoubletFinder runs on the already-preprocessed Seurat object from STEP 3.
# Parameter sweep to find the optimal pK (function names dropped the '_v3'
# suffix in the Seurat 5 compatible release).
sweep_list  <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
bcmvn       <- find.pK(sweep_stats)
 
optimal_pk <- as.numeric(as.character(
  bcmvn$pK[which.max(bcmvn$BCmetric)]
))
 
# Expected doublet count from the 10x multiplet-rate table: ~0.8% per 1,000
# cells RECOVERED for standard Next GEM 3' chemistry (v2/v3/v3.1). The rate is
# set at the wet-lab loading step, not measured from the data, so it is an
# assumption you supply. IMPORTANT: newer kits are lower -- use ~0.4% per 1,000
# (i.e. 0.004 below) for GEM-X or HT chemistry, or you will over-remove cells.
n_cells        <- ncol(seurat_obj)
homotypic_prop <- modelHomotypic(seurat_obj$seurat_clusters)  # undetectable same-type doublets
nExp           <- round((0.008 * n_cells / 1000) * n_cells)
nExp_adj       <- round(nExp * (1 - homotypic_prop))
 
seurat_obj <- doubletFinder(
  seurat_obj,
  PCs        = 1:30,
  pN         = 0.25,
  pK         = optimal_pk,
  nExp       = nExp_adj,
  reuse.pANN = NULL,   # use NULL, not FALSE: FALSE triggers an
                       # "xtfrm.data.frame ... cannot xtfrm data frames"
                       # error in current DoubletFinder versions
  sct        = FALSE
)
 
# DoubletFinder writes a column whose name encodes the parameters; grab it
df_class_col <- grep("^DF.classifications", colnames(seurat_obj[[]]), value = TRUE)
seurat_obj$DoubletFinder <- seurat_obj[[df_class_col]][, 1]
```

### Running scDblFinder

```R
#-----------------------------------------------
# STEP 9: Doublet detection with scDblFinder
#-----------------------------------------------
 
# Convert the Seurat object to SCE; metadata (including seurat_clusters) carries over
sce_doublet <- as.SingleCellExperiment(seurat_obj)
 
# No parameter tuning required; pass clusters for cluster-aware simulation
sce_doublet <- scDblFinder(sce_doublet, clusters = "seurat_clusters")
 
# Transfer the calls back onto the Seurat object for a side-by-side comparison
seurat_obj$scDblFinder       <- sce_doublet$scDblFinder.class
seurat_obj$scDblFinder_score <- sce_doublet$scDblFinder.score
```

### Comparing the Two Sets of Calls

The cleanest comparison is a cross-tabulation: how many cells does each tool call a doublet, and how often do they agree?

```R
#-----------------------------------------------
# STEP 10: Cross-tabulate the two doublet callers
#-----------------------------------------------
 
# Note the case difference: DoubletFinder uses "Doublet"/"Singlet",
# scDblFinder uses "doublet"/"singlet"
doublet_crosstab <- table(
  DoubletFinder = seurat_obj$DoubletFinder,
  scDblFinder   = seurat_obj$scDblFinder
)
print(doublet_crosstab)
 
# Define a consensus category
seurat_obj$doublet_consensus <- dplyr::case_when(
  seurat_obj$DoubletFinder == "Doublet"&
    seurat_obj$scDblFinder == "doublet"~ "Both",
  seurat_obj$DoubletFinder == "Doublet"|
    seurat_obj$scDblFinder == "doublet"~ "One only",
  TRUE~ "Neither"
)
 
##              scDblFinder
## DoubletFinder singlet doublet
##       Doublet     106     588
##       Singlet    8654     378
```

Visualising the calls on the shared UMAP shows whether the two tools flag the same regions:

```R
#-----------------------------------------------
# STEP 11: Visualise and compare on the shared UMAP
#-----------------------------------------------
 
p_df <- DimPlot(seurat_obj, group.by = "DoubletFinder",
                cols = c("Singlet"= "#90E0EF", "Doublet"= "#EF233C")) +
  ggtitle("DoubletFinder")
 
p_scdbl <- DimPlot(seurat_obj, group.by = "scDblFinder",
                   cols = c("singlet"= "#90E0EF", "doublet"= "#EF233C")) +
  ggtitle("scDblFinder")
 
p_consensus <- DimPlot(seurat_obj, group.by = "doublet_consensus",
                       cols = c("Neither"= "#CED4DA",
                                "One only"= "#FFB703",
                                "Both"= "#D00000")) +
  ggtitle("Consensus")
 
p_doublet_panel <- p_df | p_scdbl | p_consensus
ggsave("qc_compare/plots/03_doublet_caller_comparison.png",
       p_doublet_panel, width = 16, height = 5, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2026/06/03_doublet_caller_comparison-scaled.png?resize=2048%2C640&ssl=1)

> **What the real data shows:** On `Healthy_1`, the two callers agree exactly where you would hope. Both light up the dense cluster at the bottom-center of the UMAP — a classic doublet hotspot where blended transcriptomes pile up — and the consensus panel shows that region saturated with dark-red “Both” calls. They also both flag the thin bridges connecting otherwise separate clusters and parts of the top-right population. The disagreement, shown as the scattered gold “One only” points, sits at the margins: scDblFinder is somewhat more liberal across the large left-hand cluster, calling borderline cells that DoubletFinder leaves as singlets. This is the expected and reassuring pattern — strong agreement on the obvious doublets, divergence on the ambiguous ones. The consensus panel makes the safe play obvious: the “Both” cells are your highest-confidence doublets, and if you want to minimise false removals, dropping only those is the conservative choice. If you want a more thorough clean-up, trust scDblFinder’s wider net.

## The Decision Guide

This is the part you came for. Here is how to choose.
### Ambient RNA: SoupX vs DecontX

| Consideration | SoupX | DecontX |
| --- | --- | --- |
| Needs raw matrix / empty droplets | Yes (reads soup from empty droplets) | No (infers from cell populations) |
| Works on a filtered matrix alone | Awkward (must estimate soup manually) | Yes, natively |
| Contamination estimate | One global fraction for the whole sample | Per-cell fraction |
| Statistical model | Non-parametric soup profile | Bayesian mixture model |
| Typical behavior | Gentle, uniform correction | More aggressive, cluster-adaptive correction |
| Native object | Matrix (drops into Seurat easily) | SingleCellExperiment |
| Output | Corrected integer counts | Decontaminated counts + per-cell contamination |
| Speed | Fast | Slower (iterative estimation) |
| Best when | You have full Cell Ranger output and want a transparent, conservative correction | You only have a filtered matrix, or want per-cell contamination scores and cluster-adaptive cleanup |

**Recommendation:** If you ran Cell Ranger yourself and have the raw matrix, **SoupX** is the pragmatic default — it is fast, transparent, conservative, and slots directly into a Seurat workflow. Its single global estimate is a feature when you want a light, predictable touch. If you inherited a public dataset with only a filtered matrix, or you want a per-cell contamination score (useful as a QC covariate) and cluster-level adaptivity, reach for **DecontX** — but remember it can correct aggressively, so sanity-check the result against lineage-restricted marker genes, as we saw with MALAT1 above. Both are well-validated; the choice is driven by your inputs, the granularity you need, and how conservative you want to be — not by one being more “accurate” than the other.

### Doublets: DoubletFinder vs scDblFinder

| Consideration | DoubletFinder | scDblFinder |
| --- | --- | --- |
| Parameter tuning | Required (pK sweep, expected nExp) | None (self-tuning) |
| Speed | Slow (parameter sweep per sample) | Fast |
| Multi-sample handling | Manual, one sample at a time | Built in via `samples` argument |
| Native object | Seurat | SingleCellExperiment |
| Benchmark performance | Good | Among top performers (Xi and Li, 2021) |
| Maintenance | Community-maintained GitHub package | Active Bioconductor package |
| Best when | You want fine manual control, or are reproducing a Seurat-native pipeline | You want a fast, tuning-free default — most cases |

**Recommendation:** For most users and most projects, **scDblFinder** is the better default — it is faster, requires no parameter tuning, handles multiple samples cleanly, and benchmarks well. Use **DoubletFinder** when you specifically want manual control over the parameter sweep, or when you are matching an existing Seurat-native pipeline. When stakes are high, run both and remove only the consensus doublets (the “Both” category above) to minimize false positives.

## Best Practices

1. **Run ambient correction before doublet detection.** Ambient RNA can make a cell look like a blend of populations, which inflates doublet scores. Cleaning the soup first gives the doublet caller a clearer signal.
2. **Always detect doublets per sample, never on integrated data.** Doublets form within a single 10x lane. Running a doublet caller across merged or integrated samples invites it to “find” cross-sample doublets that physically cannot exist. Split first, call doublets, then integrate.
3. **Cluster once and share the labels.** As we did here, a single shared clustering keeps tool comparisons fair and keeps your QC reproducible.
4. **Treat contamination fraction as a sanity check, not just a correction.** A sample with a much higher contamination estimate than its peers often had a dissociation problem. Flag it.
5. **Prefer consensus for high-stakes removal.** Removing only the doublets that two independent methods agree on minimizes the risk of discarding rare-but-real cell states that happen to sit between clusters.
6. **Keep the corrected and uncorrected counts.** Store both so you can check whether a surprising downstream result depends on the correction. Reversibility builds trust in your pipeline.

## Common Pitfalls

- **Removing too many doublets.** Aggressive thresholds delete genuine transitional cells (for example, activated lymphocytes between resting states). If a “doublet” cluster expresses a coherent biological program rather than two unrelated lineage markers, it may be a real state.
- **Forgetting the case difference in class labels.** DoubletFinder returns `"Doublet"` / `"Singlet"`; scDblFinder returns `"doublet"` / `"singlet"`. Comparing them without accounting for capitalization silently breaks your cross-tab.
- **Using the 0.8% doublet rate on newer chemistry.** The 0.8 percent per 1,000 cells figure is for standard Next GEM 3′ kits. GEM-X and HT chemistries roughly halve it to 0.4 percent per 1,000. Applying 0.8 percent to GEM-X data over-estimates `nExp` and removes real cells. Check your kit version, or sidestep the issue entirely by letting scDblFinder estimate the rate.
- **Comparing SoupX and DecontX contamination estimates cell by cell.** SoupX’s default output is one global fraction for the whole sample; DecontX’s is per-cell. Scattering one against the other produces a meaningless vertical line and a correlation of `NA`. Compare them as a distribution-versus-reference-line instead, or compare their effect on actual gene counts.
- **Over-correcting ambient RNA.** Pushing contamination estimates too high strips real lowly expressed genes. If a known marker disappears from the cell type that should express it, you have over-corrected. Be especially wary when the top soup gene is a ubiquitous high-abundance transcript like MALAT1 — its presence in the soup profile reflects abundance, not necessarily contamination.
- **Running DoubletFinder with the old `_v3` function names on Seurat 5.** The `paramSweep_v3` and `doubletFinder_v3` names were retired when the package became Seurat 5 compatible. On a current install, use `paramSweep` and `doubletFinder`.
- **Passing `reuse.pANN = FALSE` to `doubletFinder`.** A recent patch changed how this argument is handled, and `FALSE` now triggers a `cannot xtfrm data frames` error. Use `reuse.pANN = NULL` instead. This is purely a software bug, not a problem with your data.
- **Feeding DecontX a filtered matrix and expecting it to use empty droplets.** DecontX does not use empty droplets by design. That is a feature, not a bug — but do not assume it is doing the same thing SoupX does.

## Conclusions

1. **The pairs are interchangeable in purpose, not in mechanism.** SoupX and DecontX both correct ambient RNA but differ in whether they need empty droplets. DoubletFinder and scDblFinder both flag doublets but differ in tuning burden and speed.
2. **For ambient RNA, your inputs and your appetite for aggressiveness decide.** SoupX gives one conservative global correction and needs the raw matrix; DecontX gives per-cell, cluster-adaptive correction from the filtered matrix alone. On our clean PBMC sample SoupX touched the data lightly (~2.3 percent) while DecontX cut the top soup gene hard in select clusters — a difference in philosophy, not a bug. Pick SoupX for a transparent light touch, DecontX for per-cell granularity and when you lack a raw matrix.
3. **For doublets, scDblFinder is the better default for most users** — faster, tuning-free, multi-sample aware, and strong in benchmarks. DoubletFinder remains valuable when you want manual control or a Seurat-native pipeline.
4. **Consensus beats any single tool when removal is irreversible.** Running both doublet callers and removing only their agreed-upon doublets is the safest high-stakes strategy.
5. **Tool choice is rarely the limiting factor.** Correct ordering (ambient before doublets), per-sample doublet detection, and shared clustering matter more for data quality than which of these well-validated tools you pick.

## References

1. Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience. 2020;9(12):giaa151. doi:10.1093/gigascience/giaa151
2. Yang S, Corbett SE, Koga Y, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology. 2020;21(1):57. doi:10.1186/s13059-020-1950-6
3. McGinnis CS, Murrow LM, Gartner ZJ. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Systems. 2019;8(4):329-337.e4. doi:10.1016/j.cels.2019.03.003
4. Germain PL, Lun A, Garcia Meixide C, Macnair W, Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. F1000Research. 2021;10:979. doi:10.12688/f1000research.73600.2
5. Xi NM, Li JJ. Benchmarking Computational Doublet-Detection Methods for Single-Cell RNA Sequencing Data. Cell Systems. 2021;12(2):176-194.e6. doi:10.1016/j.cels.2020.11.008
6. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology. 2024;42(2):293-304. doi:10.1038/s41587-023-01767-y
7. Lee H, Joo JY, Sohn DH, et al. Single-cell RNA sequencing reveals rebalancing of immunological response in patients with periodontitis after non-surgical periodontal therapy. Journal of Translational Medicine. 2022;20(1):504. doi:10.1186/s12967-022-03686-8 \[Dataset source — GSE174609\]

---

*This quick tip extends the QC workflow from Part 2 of the comprehensive NGS101.com single-cell RNA-seq analysis series for beginners.*