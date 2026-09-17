---
title: "How to Analyze Single-Cell RNA-seq Data – Complete Beginner's Guide Part 2: Quality Control and Cell Filtering"
source: "https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/"
author:
  - "[[Lei]]"
published: 2025-12-17
created: 2026-07-12
description: "Introduction: Learning QC Through a Single-Sample Deep Dive Quality control in single-cell RNA sequencing is complex, with multiple layers of filtering and validation. Before tackling multi-sample experiments, it’s essential to understand the QC workflow thoroughly using a single sample. This focused approach allows you to: IMPORTANT: This QC workflow should be applied independently to each […]"
tags:
  - "clippings"
---
## Introduction: Learning QC Through a Single-Sample Deep Dive

Quality control in single-cell RNA sequencing is complex, with multiple layers of filtering and validation. Before tackling multi-sample experiments, it’s essential to understand the QC workflow thoroughly using a single sample. This focused approach allows you to:

- **Understand each QC step deeply** without the complexity of batch effects
- **Master Seurat 5’s built-in plotting functions** designed specifically for scRNA-seq
- **Build intuition** about what good vs. problematic data looks like
- **Troubleshoot effectively** by seeing the impact of each filtering decision
- **Establish a template** you can apply to any sample

> **IMPORTANT:** This QC workflow should be **applied independently to each sample** in your dataset. Different samples may have different quality profiles, and thresholds should be adjusted based on visual inspection of each sample’s QC plots.

### What Makes Single-Cell Data Quality Control Unique

Unlike bulk RNA-seq where you average across millions of cells, single-cell RNA-seq requires cell-by-cell quality assessment. Each of the thousands of cells in your sample may have:

- Different RNA capture efficiency
- Variable sequencing depth
- Potential technical artifacts (doublets, empty droplets, dying cells)
- Biological heterogeneity that must be preserved

**The QC challenge**: Remove technical noise while preserving biological signal.

### The QC Workflow We’ll Follow

```xml
Raw 10x Data (per sample)
    ↓
1. Load Cell Ranger Output
    ↓
2. Empty Droplet Detection (DropletUtils)
    ↓
3. Ambient RNA Correction (SoupX) - optional if <5%
    ↓
4. Doublet Detection (scDblFinder)
    ↓
5. Cell-Level QC Filtering (manual threshold setting)
    ↓
6. Gene-Level QC Filtering
    ↓
7. Normalization & Variable Feature Selection
    ↓
Clean, QC-filtered Sample
(Ready for integration with other samples & clustering)
```

**Apply this workflow to EACH sample separately**, then integrate the clean samples.

### Why Seurat 5 for This Tutorial

Seurat 5 represents a major advancement in single-cell analysis:

- **Layer-based architecture**: Better data organization and memory management
- **Built-in visualization**: Specialized plotting functions for scRNA-seq data
- **Scalability**: Handles datasets from 100 to millions of cells
- **Integration**: Seamless workflow from QC through cell type annotation
- **Active development**: Regular updates with latest methods

> **For Beginners**: Seurat is the most widely used R package for single-cell analysis, with extensive documentation, tutorials, and community support. Learning Seurat 5 gives you a transferable skill applicable to most scRNA-seq projects.

## Dataset Overview: GSE174609 PBMC Study

We’re working with the **GSE174609 dataset**, which examines immune cell responses in periodontitis patients. This study includes **12 samples total**.

**For this tutorial**, we’ll perform QC on **one sample: Healthy\_1** (SRR14575500) to demonstrate the complete workflow.

> **Reference:** See [Part 1 of this tutorial series](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-1-from-fastq-to-count-matrix/) for details on processing FASTQ files through Cell Ranger to generate the count matrices we’ll use here.

### Why PBMCs Are Ideal for Learning

PBMCs are the “gold standard” for single-cell tutorials because they:

- Contain well-characterized cell types
- Have clear marker genes (CD3D for T cells, CD14 for monocytes, etc.)
- Show moderate complexity (not too simple, not too complex)
- Are relatively easy to dissociate (good quality data)
- Have abundant reference datasets for comparison

### Key Quality Metrics to Watch

Throughout this tutorial, we’ll monitor:

- **nCount\_RNA**: Total UMI counts per cell (expect: 1,000-50,000)
- **nFeature\_RNA**: Unique genes detected per cell (expect: 500-5,000)
- **percent.mt**: Mitochondrial gene percentage (expect: <5% for healthy cells)
- **Doublet rate**: Expected ~0.5-2% for our cell count
- **Ambient RNA contamination**: Expected <5% for well-prepared PBMCs

## Setting Up Your R Environment for Seurat 5

### System Requirements

- **R version**: ≥4.0.0 (R ≥4.3.0 recommended)
- **RAM**: 8GB minimum, 16GB+ recommended
- **Storage**: ~2GB for this single sample
- **Operating System**: Windows, macOS, or Linux

### Installing Required Packages

```
#-----------------------------------------------
# Installation: Required packages for single-sample QC
#-----------------------------------------------
 
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))
 
# Install Seurat 5 (the star of this tutorial!)
install.packages("Seurat")
 
# Install SeuratObject (Seurat's data structure package)
install.packages("SeuratObject")
 
# Install visualization packages
install.packages(c(
  "ggplot2",        # For custom plots when needed
  "patchwork",      # For combining Seurat plots
  "dplyr"# Data manipulation
))
 
# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
 
# Install Bioconductor packages for QC
BiocManager::install(c(
  "DropletUtils",           # Empty droplet detection
  "scDblFinder",            # Doublet detection (better than DoubletFinder for Seurat 5)
  "SingleCellExperiment"# Data structure for some QC tools
))
 
# Install SoupX for ambient RNA correction
install.packages("SoupX")
```

### Loading Libraries and Configuration

```
#-----------------------------------------------
# STEP 1: Load libraries and configure environment
#-----------------------------------------------
 
# Core single-cell analysis
library(Seurat)              # Main scRNA-seq toolkit
library(SeuratObject)        # Seurat object infrastructure
 
# QC-specific packages
library(DropletUtils)        # Empty droplet detection
library(scDblFinder)         # Doublet detection
library(SoupX)               # Ambient RNA correction
library(SingleCellExperiment) # SCE objects
 
# Visualization
library(ggplot2)             # Custom plots when needed
library(patchwork)           # Combine plots
library(dplyr)               # Data manipulation
 
# Set working directory
setwd("~/GSE174609_scRNA/QC_analysis")
 
# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("qc_metrics", showWarnings = FALSE)
dir.create("filtered_data", showWarnings = FALSE)
 
# Set random seed for reproducibility
set.seed(42)
 
# Verify Seurat 5 installation
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")
 
# Check for Seurat 5 features
if (packageVersion("Seurat") >= "5.0.0") {
  cat("✓ Seurat 5 detected - layer-based architecture available\n")
} else{
  warning("Please upgrade to Seurat 5 for this tutorial")
}
```

**Output Example:**

```
Seurat version: 5.3.1 
R version: R version 4.5.1 (2025-06-13) 
Working directory: ~/GSE174609_scRNA/QC_analysis 
 
✓ Seurat 5 detected - layer-based architecture available
```

## Loading and Initial Inspection of 10x Data

Cell Ranger (10x Genomics’ processing pipeline) outputs three key files that we’ll load:

- **matrix.mtx.gz**: Sparse UMI count matrix
- **features.tsv.gz**: Gene names and IDs
- **barcodes.tsv.gz**: Cell barcodes

### Understanding Your Cell Ranger Output Structure

```
~/GSE174609_scRNA/cellranger_output/Healthy_1/outs/
├── filtered_feature_bc_matrix/    # Cell Ranger's filtered cells (~3,000-5,000)
│   ├── barcodes.tsv.gz           # Already filtered for empty droplets
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── raw_feature_bc_matrix/         # ALL droplets (~20,000-100,000)
│   ├── barcodes.tsv.gz           # Includes empty droplets + cells
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── web_summary.html               # Cell Ranger QC report
```

**Key difference explained in [Part 1](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-1-from-fastq-to-count-matrix/):**

- **`filtered_feature_bc_matrix/`**: Cell Ranger already removed empty droplets and applied basic ambient RNA correction
- **`raw_feature_bc_matrix/`**: No filtering – all droplets captured during sequencing

**Your choice:**

- Use **raw** → Perform Steps 4-5 (EmptyDrops + SoupX) for maximum control
- Use **filtered** → Skip Steps 4-5, trust Cell Ranger’s filtering

This tutorial uses **raw** for educational completeness.

### Loading the Data into Seurat

```
#-----------------------------------------------
# STEP 2: Load 10x Genomics data
#-----------------------------------------------
 
# Define paths
sample_name <- "Healthy_1"
cellranger_output <- "~/GSE174609_scRNA/cellranger_output/Healthy_1/outs"
 
# OPTION A: Load RAW matrix (what this tutorial demonstrates)
counts <- Read10X(data.dir = file.path(cellranger_output, "raw_feature_bc_matrix"))
 
# OPTION B: Load FILTERED matrix (alternative - skip Steps 4-5 if using this)
# counts <- Read10X(data.dir = file.path(cellranger_output, "filtered_feature_bc_matrix"))
 
# Create Seurat object
# If raw: contains ALL droplets (~20,000-100,000)
# If filtered: contains only cells Cell Ranger called (~3,000-5,000)
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = sample_name,
  min.cells = 0,      # Don't filter genes yet
  min.features = 0    # Don't filter cells yet - we'll do this properly
)
 
# Add sample metadata (consistent with Part 1 tutorial)
seurat_obj$sample_id <- "Healthy_1"
seurat_obj$sra_id <- "SRR14575500"
seurat_obj$condition <- "Healthy"
seurat_obj$patient_id <- "Donor_1"
seurat_obj$time_point <- NA# Not applicable for healthy donors
 
# Quick check
cat("Loaded:", ncol(seurat_obj), "droplets ×", nrow(seurat_obj), "genes\n")
cat("(Most are empty - EmptyDrops will filter in Step 4)\n")
```

**Output Example:**

```
Loaded: 1389510 droplets × 38606 genes
(Most are empty - EmptyDrops will filter in Step 4)
```

> **Using Raw vs Filtered Matrix:**
> 
> **If you loaded RAW (recommended for learning):**
> 
> - Continue with Steps 4-5 (EmptyDrops + SoupX)
> - More control over QC, can detect rare cell types
> - EmptyDrops is statistically more sophisticated than Cell Ranger’s UMI threshold
> 
> **If you loaded FILTERED (faster alternative):**
> 
> - **Skip Steps 4 and 5** – Cell Ranger already did this
> - Jump directly to Step 6 (Doublet Detection)
> - Trust Cell Ranger’s default filtering parameters
> 
> **Important:** Whichever you choose for Healthy\_1, **use the same approach for all 12 samples** to maintain consistency across your dataset.

### Initial Data Exploration

Before filtering, let’s understand the raw data:

```
#-----------------------------------------------
# STEP 3: Initial data exploration
#-----------------------------------------------
 
# Basic statistics
cat("Total droplets:", ncol(seurat_obj), "\n")
cat("Total genes:", nrow(seurat_obj), "\n")
 
# Quick statistics
sparsity <- 1 - (sum(LayerData(seurat_obj, layer = "counts") > 0) / 
                 (nrow(seurat_obj) * ncol(seurat_obj)))
cat("Matrix sparsity:", round(sparsity * 100, 1), "%\n")
 
# UMI distribution - will be bimodal (empty droplets + real cells)
umi_counts <- colSums(LayerData(seurat_obj, layer = "counts"))
cat("Median UMI per droplet:", median(umi_counts), "\n")
cat("Droplets with >500 UMI:", sum(umi_counts > 500), 
    "(likely cells, will be validated by EmptyDrops)\n")
```

**Output Example:**

```
Total droplets: 1389510 
Total genes: 38606 
Matrix sparsity: 99.9 %
Median UMI per droplet: 1 
Droplets with >500 UMI: 9883 (likely cells, will be validated by EmptyDrops)
```

> **Understanding Raw Data**: The raw matrix contains ALL droplets – mostly empty (median UMI ~50-200), plus real cells (UMI >1000). EmptyDrops will statistically separate these based on gene expression profiles, not just UMI counts.

## Empty Droplet Detection with DropletUtils

> **Note:** If you loaded the **filtered matrix** in Step 2, **skip this step** (and Step 5). Cell Ranger has already performed empty droplet removal. Jump directly to Step 6 (Doublet Detection).

Not every droplet captured by the 10x system contains a real cell. Some contain only ambient RNA from lysed cells. While Cell Ranger applies basic UMI-threshold filtering, EmptyDrops uses a more sophisticated statistical approach.

### The Empty Droplet Problem

**Three types of droplets:**

1. **Cell-containing droplets**: High UMI counts, diverse gene expression
2. **Empty droplets**: Low UMI counts, only ambient RNA
3. **Damaged cells**: Intermediate UMI counts (need filtering later)

**The challenge**: Rare cell types might have low RNA content, so simple UMI thresholds can remove real cells.

### Statistical Detection with EmptyDrops

EmptyDrops uses a statistical approach:

1. Estimates the ambient RNA profile from low-count droplets
2. Tests each droplet: “Is this expression profile different from ambient RNA alone?”
3. Assigns FDR (false discovery rate) to each droplet
4. Calls droplets as cells if FDR < 0.01
```
#-----------------------------------------------
# STEP 4: Empty droplet detection
#-----------------------------------------------
 
# Convert Seurat object to SingleCellExperiment for EmptyDrops
# We already loaded the raw matrix in Step 2
sce <- as.SingleCellExperiment(seurat_obj)
 
# Track initial droplet count for QC summary
initial_droplet_count <- ncol(sce)
initial_gene_count <- nrow(sce)
 
# Run EmptyDrops
set.seed(100)
empty_results <- emptyDrops(
  m = counts(sce),
  lower = 100,        # Droplets with <100 UMI assumed empty (for estimating ambient RNA)
  niters = 10000,     # Iterations for p-value calculation
  test.ambient = TRUE# Test if droplets differ from ambient RNA profile
)
 
# Identify cells (FDR < 0.01)
is_cell <- empty_results$FDR < 0.01
is_cell[is.na(is_cell)] <- FALSE# Treat NA as empty (very low UMI droplets)
validated_barcodes <- colnames(sce)[is_cell]
 
# Brief feedback
cat("Droplets tested:", ncol(sce), "\n")
cat("Cells called:", sum(is_cell), 
    "(", round(sum(is_cell)/ncol(sce)*100, 1), "%)\n")
cat("Empty droplets removed:", sum(!is_cell), 
    "(", round(sum(!is_cell)/ncol(sce)*100, 1), "%)\n")
 
# Visual QC - the critical check!
empty_df <- data.frame(
  total_umi = colSums(counts(sce)),
  is_cell = is_cell
)
 
p1 <- ggplot(empty_df, aes(x = log10(total_umi + 1), fill = is_cell)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(
    values = c("TRUE"= "#2E86AB", "FALSE"= "#A23B72"),
    labels = c("Empty Droplet", "Cell"),
    name = "Classification"
  ) +
  labs(
    title = "EmptyDrops: Cell vs Empty Droplet Detection",
    subtitle = paste(sum(is_cell), "cells called from", ncol(sce), "total droplets"),
    x = "log10(UMI + 1)",
    y = "Number of Droplets"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14))
 
ggsave("plots/01_empty_droplets.png", p1, width = 10, height = 6, dpi = 300)
 
# NOW filter the Seurat object to keep only validated cells
# IMPORTANT: Save raw counts BEFORE filtering (needed for SoupX in Step 5)
raw_counts_all_droplets <- LayerData(seurat_obj, layer = "counts")
 
seurat_obj <- subset(seurat_obj, cells = validated_barcodes)
 
cat("After EmptyDrops:", ncol(seurat_obj), "cells retained\n")
```

**Output Example:**

```
Droplets tested: 1389510 
Cells called: 9603 ( 0.7 %)
Empty droplets removed: 1379907 ( 99.3 %)
After EmptyDrops: 9603 cells retained
```

> **Key Concept**: EmptyDrops is more sophisticated than simple UMI thresholds. A droplet with 500 UMI expressing 50 genes in a cell-type-specific pattern is likely a real cell. A droplet with 500 UMI expressing housekeeping genes is likely ambient RNA.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/01_empty_droplets-scaled.png?resize=2048%2C1229&ssl=1)

## Ambient RNA Correction with SoupX

> **Note:** If you loaded the **filtered matrix** in Step 2, **skip this step**. Cell Ranger has already applied basic ambient RNA correction. Jump directly to Step 6 (Doublet Detection).

Even droplets containing cells are contaminated with ambient RNA—the “soup” of RNA floating in the suspension from lysed cells.

### Why Ambient RNA Matters

**The problem:**

- Lysed cells release their RNA into the suspension
- This ambient RNA enters all droplets
- Creates false positive gene expression
- Especially affects highly expressed genes (hemoglobin in blood, myelin in brain)

**The impact:**

- Inflates marker gene expression
- Creates false cell-type assignments
- Reduces specificity of differential expression

### Estimating and Removing Contamination

> **Note:** This step is optional. First check the contamination level – if it’s NULL or <5%, you can skip correction.

```
#-----------------------------------------------
# STEP 5: Ambient RNA correction with SoupX
#-----------------------------------------------
 
# Prepare SoupX channel
# tod (table of droplets) = raw counts from ALL droplets (saved in Step 4)
# toc (table of cells) = filtered counts from validated cells only
tod <- raw_counts_all_droplets
toc <- LayerData(seurat_obj, layer = "counts")
sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)
 
# Quick clustering for contamination estimation
temp_obj <- seurat_obj
temp_obj <- NormalizeData(temp_obj, verbose = FALSE)
temp_obj <- FindVariableFeatures(temp_obj, nfeatures = 2000, verbose = FALSE)
temp_obj <- ScaleData(temp_obj, verbose = FALSE)
temp_obj <- RunPCA(temp_obj, npcs = 30, verbose = FALSE)
temp_obj <- FindNeighbors(temp_obj, dims = 1:30, verbose = FALSE)
temp_obj <- FindClusters(temp_obj, resolution = 0.8, verbose = FALSE)
sc <- setClusters(sc, setNames(as.character(temp_obj$seurat_clusters), colnames(temp_obj)))
 
# Estimate contamination
sc <- tryCatch({
  autoEstCont(sc, verbose = FALSE)
}, error = function(e) {
  cat("Note: autoEstCont failed\n")
  sc$fit$rho <- NULL
  return(sc)
})
 
contamination_fraction <- sc$fit$rho
 
# Check contamination level
cat("Contamination fraction:", 
    ifelse(is.null(contamination_fraction), "NULL (very clean sample)", 
           paste0(round(contamination_fraction * 100, 2), "%")), "\n")
 
# Decision: Skip or correct based on contamination level
if (is.null(contamination_fraction) || contamination_fraction < 0.05) {
  cat("→ Contamination is NULL or <5% - skipping SoupX correction\n")
  cat("   Your sample is clean! Proceeding with original counts.\n")
} else{
  cat("→ Contamination is", round(contamination_fraction * 100, 2), "% (≥5%)\n")
  cat("   Applying SoupX correction...\n")
 
  # Apply correction
  suppressWarnings({
    corrected_counts <- adjustCounts(sc)
  })
  seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = corrected_counts)
  cat("   ✓ SoupX correction applied\n")
}
```

**Output Example:**

```
Contamination fraction: NULL (very clean sample) 
→ Contamination is NULL or &lt;5% - skipping SoupX correction
   Your sample is clean! Proceeding with original counts.
```

> **When to Skip SoupX:** If contamination is NULL (estimation failed) or <5%, your sample is clean and correction is unnecessary. Values >5% indicate ambient RNA contamination worth correcting, especially in tissues prone to cell lysis (brain, tumor).

## Doublet Detection with scDblFinder

Doublets occur when two cells are captured in the same droplet. They appear as artificial “intermediate” cell types and must be removed.

### Understanding Doublets

**How doublets form:**

- 10x captures ~1% doublets per 1,000 cells loaded
- Higher loading = more doublets (10,000 cells ~6 – 9% doublets)
- Completely random—any two cell types can combine

**Why they’re problematic:**

- Create false “hybrid” cell types
- Inflate cell counts
- Confound trajectory analysis
- Mess up differential expression

### Computational Doublet Detection

```
#-----------------------------------------------
# STEP 6: Doublet detection with scDblFinder
#-----------------------------------------------
 
# Join layers and convert to SCE
seurat_obj <- JoinLayers(seurat_obj)
sce <- as.SingleCellExperiment(seurat_obj)
 
# Run scDblFinder (suppress expected warnings about empty layers and xgboost)
set.seed(42)
suppressWarnings({
  sce <- scDblFinder(sce, dbr = NULL, verbose = FALSE)
})
 
# Extract results
seurat_obj$doublet_class <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
seurat_obj$doublet_score <- sce$scDblFinder.score
 
# Brief feedback
n_doublets <- sum(seurat_obj$doublet_class == "Doublet")
doublet_rate <- n_doublets / ncol(seurat_obj) * 100
cat("Doublets:", n_doublets, "(", round(doublet_rate, 1), "%)\n")
 
# Track cells before doublet removal for QC summary
cells_before_doublet_removal <- ncol(seurat_obj)
 
# Quick UMAP for visualization
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
 
# Visualize doublets on UMAP
p3 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "doublet_class",
  cols = c("Singlet"= "#90E0EF", "Doublet"= "#EF233C"),
  pt.size = 1
) +
  labs(title = paste0("Doublet Detection (", round(doublet_rate, 1), "%)")) +
  theme(plot.title = element_text(face = "bold", size = 14))
 
ggsave("plots/02_doublets_umap.png", p3, width = 8, height = 7, dpi = 300)
 
# Remove doublets
seurat_obj <- subset(seurat_obj, subset = doublet_class == "Singlet")
cat("After doublet removal:", ncol(seurat_obj), "cells\n")
```

**Output Example:**

```
Doublets: 1013 ( 10.5 %)
After doublet removal: 8590 cells
```

> **Note on Warnings:** scDblFinder may generate harmless warnings about empty Seurat layers (normal – it only needs raw counts) and xgboost parameter deprecation (package version mismatch). These are suppressed with `suppressWarnings()` for cleaner output but don’t affect doublet detection accuracy.
> 
> **Expected Doublet Rates**: For ~3,000-5,000 cells, expect 0.5-2% doublets. Rates >10% suggest problems with cell loading concentration or the doublet detector being too sensitive. Given that our example sample includes over 9,800 cells, a 10% rate is considered acceptable.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/02_doublets_umap.png?resize=2048%2C1792&ssl=1)

## Cell-Level Quality Control: Manual Threshold Setting

Now we filter individual cells based on their quality metrics. **The key is visual inspection** – look at the distributions and set biologically meaningful thresholds rather than relying on automatic statistical methods.

```
#-----------------------------------------------
# STEP 7: Calculate QC metrics, visualize, and filter cells
#-----------------------------------------------
 
# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
 
# Summary statistics to guide threshold setting
cat("\nQC Metric Distributions:\n")
cat("nCount_RNA (UMI):\n")
cat("  Median:", median(seurat_obj$nCount_RNA), 
    "| Q1-Q3:", quantile(seurat_obj$nCount_RNA, 0.25), "-", 
    quantile(seurat_obj$nCount_RNA, 0.75), "\n")
 
cat("nFeature_RNA (genes):\n")
cat("  Median:", median(seurat_obj$nFeature_RNA),
    "| Q1-Q3:", quantile(seurat_obj$nFeature_RNA, 0.25), "-",
    quantile(seurat_obj$nFeature_RNA, 0.75), "\n")
 
cat("percent.mt:\n")
cat("  Median:", round(median(seurat_obj$percent.mt), 2), "%",
    "| 95th percentile:", round(quantile(seurat_obj$percent.mt, 0.95), 2), "%\n")
 
cat("percent.ribo:\n")
cat("  Median:", round(median(seurat_obj$percent.ribo), 2), "%\n")
 
# Visualize distributions - THE CRITICAL STEP
p4 <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  pt.size = 0.1
) &
  theme(plot.title = element_text(face = "bold"))
 
ggsave("plots/03_qc_violins.png", p4, width = 16, height = 4, dpi = 300)
 
# Scatter plots showing metric relationships
p5 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(title = "UMI vs Genes Detected")
 
p6 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(title = "UMI vs Mitochondrial %")
 
p7 <- FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "percent.ribo") +
  labs(title = "Mitochondrial % vs Ribosomal %")
 
p_scatter <- p5 + p6 + p7
ggsave("plots/04_qc_scatter.png", p_scatter, width = 15, height = 5, dpi = 300)
 
cat("\n→ EXAMINE plots/03_qc_violins.png and plots/04_qc_scatter.png\n")
cat("→ Look for:\n")
cat("   - Bimodal distributions (good vs bad cells)\n")
cat("   - Outliers in scatter plots\n")
cat("   - Relationship between metrics\n\n")
 
# Set thresholds based on visual inspection
# MANUALLY ADJUST THESE based on what you see in the plots!
 
cat("=== Setting Filtering Thresholds ===\n")
cat("Based on visual inspection of the plots above:\n\n")
 
# For our case - adjust these for your specific sample!
nfeature_min <- 500    # Minimum genes per cell
nfeature_max <- 5000   # Maximum genes (above this = likely doublets)
ncount_min <- 800      # Minimum UMI per cell
ncount_max <- 20000    # Maximum UMI (above this = likely doublets)
mt_thresh <- 10        # Maximum mitochondrial % (adjust based on tissue)
 
cat("Set thresholds (adjust based on YOUR data):\n")
cat("  nFeature_RNA: [", nfeature_min, ",", nfeature_max, "]\n")
cat("  nCount_RNA: [", ncount_min, ",", ncount_max, "]\n")
cat("  percent.mt: <", mt_thresh, "%\n\n")
 
cat("Guidelines for threshold setting:\n")
cat("  • nFeature_RNA min: Where does the low-quality tail end? (typically 500-1000)\n")
cat("  • nFeature_RNA max: Where do potential doublets start? (typically 4000-7000)\n")
cat("  • percent.mt: PBMCs <10%, Brain <5%, Tumor <20%\n")
cat("  • If distributions are bimodal, set threshold between modes\n")
cat("  • Aim to remove no more than 25% of cells \n\n")
 
# Visualize thresholds on data
qc_df <- data.frame(
  nCount_RNA = seurat_obj$nCount_RNA,
  nFeature_RNA = seurat_obj$nFeature_RNA,
  percent.mt = seurat_obj$percent.mt
)
 
qc_df$pass_qc <- (
  qc_df$nCount_RNA >= ncount_min &
  qc_df$nCount_RNA <= ncount_max &
  qc_df$nFeature_RNA >= nfeature_min &
  qc_df$nFeature_RNA <= nfeature_max &
  qc_df$percent.mt < mt_thresh
)
 
p8 <- ggplot(qc_df, aes(x = log10(nCount_RNA + 1), y = log10(nFeature_RNA + 1), 
                        color = pass_qc)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = log10(c(ncount_min, ncount_max)), 
             linetype = "dashed", color = "red") +
  geom_hline(yintercept = log10(c(nfeature_min, nfeature_max)), 
             linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE"= "#06D6A0", "FALSE"= "#EF476F")) +
  labs(
    title = "Cell Filtering Thresholds",
    subtitle = paste0(sum(qc_df$pass_qc), " cells pass QC (", 
                      round(sum(qc_df$pass_qc)/nrow(qc_df)*100, 1), "%)"),
    x = "log10(UMI + 1)",
    y = "log10(Genes + 1)",
    color = "Pass QC"
  ) +
  theme_classic()
 
ggsave("plots/05_filtering_thresholds.png", p8, width = 8, height = 7, dpi = 300)
 
# Report filtering impact
removal_pct <- sum(!qc_df$pass_qc) / nrow(qc_df) * 100
cat("Filtering impact:\n")
cat("  Cells before:", nrow(qc_df), "\n")
cat("  Cells passing QC:", sum(qc_df$pass_qc), "\n")
cat("  Cells removed:", sum(!qc_df$pass_qc), 
    "(", round(removal_pct, 1), "%)\n\n")
 
# Breakdown by criterion
cat("Removal breakdown:\n")
cat("  Low genes (<", nfeature_min, "):", 
    sum(qc_df$nFeature_RNA < nfeature_min), "\n")
cat("  High genes (>", nfeature_max, "):", 
    sum(qc_df$nFeature_RNA > nfeature_max), "\n")
cat("  High MT% (>", mt_thresh, "%):", 
    sum(qc_df$percent.mt > mt_thresh), "\n\n")
 
# Sanity check on removal rate
if (removal_pct > 30) {
  cat("⚠ WARNING: Removing", round(removal_pct, 1), "% of cells is high!\n")
  cat("   Consider relaxing thresholds (especially MT%)\n\n")
} elseif (removal_pct < 5) {
  cat("⚠ WARNING: Only removing", round(removal_pct, 1), "% of cells is low.\n")
  cat("   Check if you're retaining low-quality cells.\n\n")
} else{
  cat("✓ Removal rate is reasonable (typical: 10-25%)\n\n")
}
 
# Track cells before cell-level QC filtering for summary
cells_before_cell_qc <- nrow(qc_df)
 
# Apply filters
seurat_obj <- subset(
  seurat_obj,
  subset = nCount_RNA >= ncount_min &
           nCount_RNA <= ncount_max &
           nFeature_RNA >= nfeature_min &
           nFeature_RNA <= nfeature_max &
           percent.mt < mt_thresh
)
 
cat("After filtering:", ncol(seurat_obj), "cells remaining\n")
```

**Output Example:**

```
QC Metric Distributions:
nCount_RNA (UMI):
  Median: 7310.5 | Q1-Q3: 5116.25 - 9578.75 
nFeature_RNA (genes):
  Median: 2604 | Q1-Q3: 2056.5 - 3114.75 
percent.mt:
  Median: 4.17 % | 95th percentile: 18.74 %
percent.ribo:
  Median: 17.53 %
 
→ EXAMINE plots/03_qc_violins.png and plots/04_qc_scatter.png
→ Look for:
   - Bimodal distributions (good vs bad cells)
   - Outliers in scatter plots
   - Relationship between metrics
 
=== Setting Filtering Thresholds ===
Based on visual inspection of the plots above:
 
Set thresholds (adjust based on YOUR data):
  nFeature_RNA: [ 500 , 5000 ]
  nCount_RNA: [ 800 , 20000 ]
  percent.mt: &lt; 10 %
 
Guidelines for threshold setting:
  • nFeature_RNA min: Where does the low-quality tail end? (typically 500-1000)
  • nFeature_RNA max: Where do potential doublets start? (typically 4000-7000)
  • percent.mt: PBMCs &lt;10%, Brain &lt;5%, Tumor &lt;20%
  • If distributions are bimodal, set threshold between modes
  • Aim to remove no more than 25% of cells
 
Filtering impact:
  Cells before: 8590 
  Cells passing QC: 7252 
  Cells removed: 1338 ( 15.6 %)
 
Removal breakdown:
  Low genes (&lt; 500 ): 481 
  High genes (> 5000 ): 110 
  High MT% (> 10 %): 896 
 
✓ Removal rate is reasonable (typical: 10-25%)
 
After filtering: 7252 cells remaining
```

> **Critical:** This manual approach is better than automated statistical thresholds because:
> 
> - You see exactly what you’re filtering
> - You understand WHY certain cells are removed
> - You can adjust for tissue-specific biology
> - You avoid over-filtering from skewed distributions
> 
> **Always examine the plots before setting thresholds!** The violin plots show you where natural breaks occur, and scatter plots reveal quality vs doublet boundaries.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/scRNAseq_qc345-e1765680023431.png?w=1480&ssl=1)

## Gene-Level Quality Control

After filtering cells, we filter genes based on how widely they’re detected across cells.

```
#-----------------------------------------------
# STEP 8: Gene-level QC with detection threshold
#-----------------------------------------------
 
# Calculate gene detection
counts_matrix <- LayerData(seurat_obj, layer = "counts")
gene_detection <- rowSums(counts_matrix > 0)
 
gene_qc <- data.frame(
  gene = rownames(seurat_obj),
  n_cells_detected = gene_detection,
  pct_cells_detected = (gene_detection / ncol(seurat_obj)) * 100,
  is_mt = grepl("^MT-", rownames(seurat_obj)),
  is_ribo = grepl("^RP[SL]", rownames(seurat_obj)),
  is_hb = grepl("^HB[AB]", rownames(seurat_obj))
)
 
cat("Total genes:", nrow(gene_qc), "\n")
cat("MT genes:", sum(gene_qc$is_mt), 
    "| Ribo genes:", sum(gene_qc$is_ribo), 
    "| Hb genes:", sum(gene_qc$is_hb), "\n\n")
 
# Visualize gene detection
p9 <- ggplot(gene_qc, aes(x = pct_cells_detected)) +
  geom_histogram(bins = 50, fill = "#118AB2", alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(
    title = "Gene Detection Across Cells",
    subtitle = "Red line: 0.1% of cells threshold",
    x = "% of Cells Expressing Gene (log scale)",
    y = "Number of Genes"
  ) +
  theme_classic()
 
ggsave("plots/06_gene_detection.png", p9, width = 8, height = 5, dpi = 300)
 
cat("→ EXAMINE plots/06_gene_detection.png\n")
cat("   Choose threshold based on where detection drops off\n\n")
 
# Set filtering threshold - ADJUST BASED ON YOUR DATA
# Common approaches:
#   - 0.1% of cells (very lenient, keeps most genes)
#   - 1% of cells (moderate, standard approach)
#   - 3 cells minimum (absolute count, conservative)
 
min_pct_cells <- 0.1  # Gene must be in ≥0.1% of cells
# Alternatively: min_cells <- 3  # Gene must be in ≥3 cells
 
cat("Setting gene filter threshold:\n")
cat("  Minimum detection: ≥", min_pct_cells, "% of cells\n")
cat("  (Equivalent to ≥", ceiling(ncol(seurat_obj) * min_pct_cells / 100), 
    "cells with current dataset)\n\n")
 
cat("Threshold options:\n")
cat("  • 0.1% of cells: Lenient (keeps rare cell type markers)\n")
cat("  • 1% of cells: Standard (balances detection vs noise)\n")
cat("  • 3-10 cells minimum: Conservative (removes very rare genes)\n\n")
 
# Filter genes
genes_to_keep <- (gene_qc$pct_cells_detected >= min_pct_cells) & !gene_qc$is_hb
 
cat("Genes passing filter:", sum(genes_to_keep), "/", nrow(gene_qc), "\n")
cat("Genes removed:\n")
cat("  Low detection:", sum(gene_qc$pct_cells_detected < min_pct_cells), "\n")
cat("  Hemoglobin:", sum(gene_qc$is_hb), "\n\n")
 
seurat_obj <- seurat_obj[gene_qc$gene[genes_to_keep], ]
 
cat("After gene filtering:", nrow(seurat_obj), "genes remaining\n")
```

**Output Example:**

```
Total genes: 38606 
MT genes: 13 | Ribo genes: 107 | Hb genes: 3 
 
→ EXAMINE plots/06_gene_detection.png
   Choose threshold based on where detection drops off
 
Setting gene filter threshold:
  Minimum detection: ≥ 0.1 % of cells
  (Equivalent to ≥ 8 cells with current dataset)
 
Threshold options:
  • 0.1% of cells: Lenient (keeps rare cell type markers)
  • 1% of cells: Standard (balances detection vs noise)
  • 3-10 cells minimum: Conservative (removes very rare genes)
 
Genes passing filter: 21731 / 38606 
Genes removed:
  Low detection: 16872 
  Hemoglobin: 3 
 
After gene filtering: 21731 genes remaining
```

> **Why Remove Hemoglobin?** In PBMC samples, hemoglobin genes (HBA1, HBA2, HBB) are typically from contaminating red blood cells during sample preparation, not true expression in PBMCs.
> 
> **Detection Threshold:** The percentage-based approach (e.g., 0.1% of cells) scales with dataset size. For 5,000 cells, 0.1% = 5 cells. This is more robust than fixed “3 cells” cutoffs when analyzing datasets of varying sizes.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/06_gene_detection.png?resize=2048%2C1280&ssl=1)

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/scRNAseq_scatter_filtered-e1765914864276.png?w=1480&ssl=1)

## Normalization and Variable Feature Selection

With clean cells and genes, we normalize for downstream analysis.

```
#-----------------------------------------------
# STEP 9: Normalization and variable features
#-----------------------------------------------
 
# Log-normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000, verbose = FALSE)
 
# Find highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
 
top10 <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")
# Top 10 variable genes: IGLC2, IGKC, LYPD2, LINGO2, IGLC3, PTGDS, FCER1A, IGHM, BANK1, S100A9
 
# Visualize variable features
p10 <- VariableFeaturePlot(seurat_obj)
p10 <- LabelPoints(plot = p10, points = top10, repel = TRUE)
ggsave("plots/07_variable_features.png", p10, width = 10, height = 7, dpi = 300)
```

> **Variable Features**: These genes show high variability across cells (not just high expression). They’re likely involved in defining cell types and states, making them most useful for downstream clustering and cell type identification (covered in the next tutorial).

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/07_variable_features-scaled.png?resize=2048%2C1434&ssl=1)

## Saving the Clean Dataset

Now save your QC-filtered data for downstream analysis.

```
#-----------------------------------------------
# STEP 10: Save processed data and create summary
#-----------------------------------------------
 
# Save the clean Seurat object
saveRDS(seurat_obj, file = "filtered_data/Healthy_1_qc_filtered.rds")
 
# Create ONE comprehensive QC summary (instead of many CSV files)
qc_summary <- data.frame(
  sample = sample_name,
 
  # Cell counts through QC pipeline
  initial_droplets = initial_droplet_count,  # From raw matrix before EmptyDrops
  after_emptydrops = length(validated_barcodes),  # Cells called by EmptyDrops
  after_soupx = length(validated_barcodes),       # SoupX doesn't remove cells
  after_doublets = cells_before_doublet_removal - n_doublets,  # After doublet removal
  after_cell_qc = ncol(seurat_obj),  # After cell-level QC filtering
  final_cells = ncol(seurat_obj),    # Final count (same as after_cell_qc)
 
  # Gene counts
  initial_genes = initial_gene_count,  # From raw matrix
  final_genes = nrow(seurat_obj),      # After gene filtering
 
  # Key metrics (from final filtered data)
  median_umi = median(seurat_obj$nCount_RNA),
  median_genes = median(seurat_obj$nFeature_RNA),
  median_mt_pct = median(seurat_obj$percent.mt),
 
  # QC parameters
  contamination_pct = ifelse(is.null(contamination_fraction), 0, 
                             round(contamination_fraction * 100, 2)),
  doublet_rate_pct = round(doublet_rate, 2),
 
  # Filtering thresholds used
  ncount_min = ncount_min,
  ncount_max = ncount_max,
  nfeature_min = nfeature_min,
  nfeature_max = nfeature_max,
  mt_threshold = mt_thresh
)
 
# Save comprehensive summary
write.csv(qc_summary, "qc_metrics/QC_summary.csv", row.names = FALSE)
write.csv(seurat_obj@meta.data, "filtered_data/cell_metadata.csv")
 
cat("\n=== QC Pipeline Complete ===\n")
cat("Final cells:", ncol(seurat_obj), "\n")
cat("Final genes:", nrow(seurat_obj), "\n")
cat("Plots saved:", length(list.files("plots")), "files\n\n")
 
cat("Files created:\n")
cat("  • filtered_data/Healthy_1_qc_filtered.rds - Clean Seurat object\n")
cat("  • qc_metrics/QC_summary.csv - Comprehensive QC summary\n")
cat("  • filtered_data/cell_metadata.csv - Cell-level metadata\n")
```

**Output Example:**

```
=== QC Pipeline Complete ===
Final cells: 7252 
Final genes: 21731 
Plots saved: 7 files
 
Files created:
  • filtered_data/Healthy_1_qc_filtered.rds - Clean Seurat object
  • qc_metrics/QC_summary.csv - Comprehensive QC summary
  • filtered_data/cell_metadata.csv - Cell-level metadata
```

> **Important:** This QC workflow should be **applied to EACH sample in your dataset independently** before integration. Different samples may have different quality profiles and require slightly different thresholds. Always inspect the QC plots for each sample individually.

## Understanding Your QC Results

### What Good Data Looks Like

After successful QC, check your `QC_summary.csv` for:

**Cell Metrics:**

- Final cells: 3,000-5,000 (for this PBMC sample)
- Median UMI: 2,000-10,000
- Median genes: 1,000-3,000
- Median MT%: <3-5%
- Doublet rate: <2%
- Contamination: <5%

**QC Removal Rates (what to expect):**

- **10-25% removal**: Normal, healthy QC
- **<5% removal**: Too lenient – may retain low-quality cells
- **30-40% removal**: Too strict – check thresholds (especially MT%)
- **\>50% removal**: Major issue – sample quality or threshold problems

### Key Files Generated

After running this tutorial, you’ll have:

```
your_project/
├── plots/                          # Visual QC checks
│   ├── 01_empty_droplets.png      # EmptyDrops detection
│   ├── 02_doublets_umap.png       # Doublet detection
│   ├── 03_qc_violins.png          # Cell QC metrics
│   ├── 04_qc_scatter.png          # Cell QC relationships
│   ├── 05_filtering_thresholds.png # Filter visualization
│   ├── 06_gene_detection.png      # Gene QC
│   └── 07_variable_features.png   # Highly variable genes
│
├── qc_metrics/
│   └── QC_summary.csv              # ONE comprehensive summary
│
└── filtered_data/
    ├── Healthy_1_qc_filtered.rds   # Clean Seurat object (ready for the next step!)
    └── cell_metadata.csv           # Final cell information
```

### Red Flags to Watch For

**Problem signs:**

- Very high doublet rate (>10%)
- Most cells have high MT% (>10%)
- Bimodal UMI distribution without clear biological reason
- Removing >30% of cells (check MT% threshold!)
- \>30% ambient contamination

## Important: Applying QC to All 12 Samples in GSE174609

**This QC workflow must be applied to EACH of the 12 samples independently.**

### Why Per-Sample QC Matters

- Different samples may have different quality profiles
- Batch effects can create sample-specific distributions
- Thresholds should be adjusted for each sample based on visual inspection
- EmptyDrops and SoupX parameters are sample-specific
- **Matrix choice (raw vs filtered) should be consistent across all samples**

### Key Considerations

**Consistency is crucial:**

- Use the **same matrix type** (raw or filtered) for all samples
- If using raw → Run EmptyDrops + SoupX on all samples
- If using filtered → Skip EmptyDrops + SoupX on all samples

**Sample-specific adjustments:**

- MT% thresholds may vary (e.g., Patient samples might need more lenient thresholds)
- Gene/UMI thresholds should be guided by each sample’s distribution
- Doublet rates depend on cell loading (should be similar across samples)

**What NOT to do:**

- Don’t use the exact same numerical thresholds for all samples
- Don’t skip visual inspection of QC plots
- Don’t mix raw and filtered matrices across samples
- Don’t proceed with integration if any sample looks problematic

## Troubleshooting Common Issues

### Issue 1: Too Few Cells After Filtering

**Symptoms**: <500 cells remaining

**Causes:**

- Over-stringent filtering thresholds
- Poor sample quality
- Wrong tissue type assumptions

**Solutions:**

- Check Cell Ranger web\_summary.html for sequencing issues
- Review MT% threshold for your tissue type (may need to be more lenient)
- Check if your percentile-based thresholds have appropriate sanity checks

### Issue 1b: Removing Too Many Cells (>30%)

**Symptoms**: QC removes 30-40%+ of cells

**Most Common Cause**: **MT% threshold too strict**

**Diagnosis:**

```
# Check what's causing the removals
summary(seurat_obj$percent.mt)
cat("Cells with >5% MT:", sum(seurat_obj$percent.mt > 5), "\n")
cat("Cells with >10% MT:", sum(seurat_obj$percent.mt > 10), "\n")
 
# Break down removals
cat("Low genes:", sum(seurat_obj$nFeature_RNA < nfeature_min), "\n")
cat("High MT%:", sum(seurat_obj$percent.mt > mt_thresh), "\n")
```

**Solutions:**

- If most removals are MT%: **Use 10% threshold instead of 5%** for PBMCs
- If median MT% is 3-6%: Your sample is healthy, use percentile-based + sanity checks
- Check the breakdown – if >90% of removals are from one criterion, that threshold is wrong
- The updated Step 8 includes automatic sanity checks to prevent this

### Issue 2: High Doublet Rate

**Symptoms**: >10% doublets detected

**Causes:**

- High cell loading concentration
- scDblFinder being too sensitive
- Actual high doublet rate

**Solutions:**

```
# Use stricter doublet score threshold
seurat_obj$doublet_class <- ifelse(
  sce$scDblFinder.score > 0.5,  # Stricter than default
  "Doublet",
  "Singlet"
)
```

### Issue 3: Cluster Driven by QC Metrics

**Symptoms**: One cluster shows uniformly high MT% or low UMI

**Causes:**

- Incomplete filtering
- True biological stressed population

**Solutions:**

- Re-examine filtering thresholds
- Check if cluster expresses stress markers (HSP genes)
- If technical, increase filtering stringency
- If biological (e.g., apoptotic cells in tumor), annotate appropriately

### Issue 4: SoupX Fails (Often Not a Problem!)

**Symptoms**: `contamination_fraction` is NULL or error during estimation

**Causes:**

- Too few cells for reliable clustering (<100 cells)
- Very homogeneous population (all same cell type)
- Insufficient empty droplets in raw matrix
- **Sample is very clean (actually good!)**

**What Happens:**

- Tutorial automatically handles this with error handling
- Uses original counts if SoupX fails
- Continues to next step without issue

**Solutions:**

- If contamination is NULL: **Just continue** – the tutorial skips correction
- This often means your sample prep was excellent
- You can manually set contamination if you suspect it:
```
# If you know contamination exists but SoupX failed
sc$fit$rho &lt;- 0.03  # Set 3% contamination manually
corrected_counts &lt;- adjustCounts(sc)
```
- Or skip SoupX entirely – it’s optional for clean samples

## Best Practices Summary

### Essential QC Principles

**1\. Visual Inspection is Key**

- Plots are your primary QC tool
- Don’t just print metrics – visualize them
- One good plot > multiple summary statistics

**2\. Streamlined Reporting**

- Brief console feedback during analysis
- Visual QC checks via plots
- ONE comprehensive summary at the end

**3\. Data-Driven Thresholds with Sanity Checks**

- Use percentile-based thresholds
- **Apply biological sanity checks** – prevent unrealistic thresholds
- PBMCs: ≥500 genes, ≥800 UMI, <10% MT
- Warn if removal rate is >30% (too strict) or <5% (too lenient)
- Document all decisions and adjustments in your summary

**4\. Iterative QC**

- QC is not one-and-done
- Return to QC plots after clustering
- Re-evaluate if results seem off

**5\. Biological Context**

- Know your tissue’s expected metrics
- Different cell types have different profiles
- Stressed cells aren’t always technical artifacts

## Conclusion: Building a Foundation

You’ve now completed comprehensive quality control on a single PBMC sample, learning:

**Every QC step in depth** – from empty droplets to normalization  
**Seurat 5’s architecture** – layers, functions, and workflows  
**Built-in visualization tools** – when and how to use them  
**Data-driven filtering** – thresholds and validation  
**Problem diagnosis** – recognizing and fixing QC issues

### Resources for Continued Learning

**Official Documentation:**

- Seurat 5 vignettes: [https://satijalab.org/seurat/](https://satijalab.org/seurat/)
- Seurat GitHub: [https://github.com/satijalab/seurat](https://github.com/satijalab/seurat)
- DropletUtils GitHub: [https://github.com/MarioniLab/DropletUtils?tab=readme-ov-file](https://github.com/MarioniLab/DropletUtils?tab=readme-ov-file)
- SoupX GitHub: [https://github.com/constantAmateur/SoupX](https://github.com/constantAmateur/SoupX)
- scDblFinder GitHub: [https://github.com/plger/scDblFinder](https://github.com/plger/scDblFinder)

## References

1. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol. 2024;42(2):293-304. \[Seurat 5 paper\]
2. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol. 2019;15(6):e8746.
3. Lun ATL, Riesenfeld S, Andrews T, et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biol. 2019;20(1):63.
4. Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience. 2020;9(12):giaa151.
5. Germain PL, Lun A, Garcia Meixide C, Macnair W, Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. F1000Res. 2021;10:979.
6. Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control. Bioinformatics. 2021;37(7):963-967.
7. McCarthy DJ, Campbell KR, Lun ATL, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics. 2017;33(8):1179-1186.
8. Amezquita RA, Lun ATL, Becht E, et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods. 2020;17(2):137-145.
9. Zheng GXY, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun. 2017;8:14049.
10. Lee H, Joo JY, Sohn DH, et al. Single-cell RNA sequencing reveals rebalancing of immunological response in patients with periodontitis after non-surgical periodontal therapy. J Transl Med. 2022;20(1):504. \[Dataset source\]

---

*This tutorial is part of the NGS101.com series on single cell sequencing analysis. If this tutorial helped advance your research, **please comment and share your experience** to help other researchers! **Subscribe to stay updated** with our latest bioinformatics tutorials and resources.*