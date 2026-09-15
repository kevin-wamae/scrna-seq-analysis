---
title: "How to Analyze Single-Cell RNA-seq Data – Complete Beginner’s Guide Part 3: Integration and Clustering"
source: "https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-3-integration-and-clustering/"
author:
  - "[[Lei]]"
published: 2025-12-25
created: 2026-07-12
description: "Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis In Part 1 and Part 2 of this tutorial series, we processed PBMC samples from the GSE174609 dataset through the complete pipeline: from raw FASTQ files to quality-controlled count matrices. Now we face a critical question: How do we analyze multiple samples together to identify cell types […]"
tags:
  - "clippings"
---
## Table of Contents
- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis|1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.1 The Batch Effect Problem in Single-Cell Data|1.1 The Batch Effect Problem in Single-Cell Data]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.2 What is Data Integration?|1.2 What is Data Integration?]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?|1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.1 The Bulk RNA-seq vs Single-Cell RNA-seq Visualisation Paradigm Shift|1.3.1 The Bulk RNA-seq vs Single-Cell RNA-seq Visualisation Paradigm Shift]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.2 What You’re Actually Seeing When Samples “Mix” on the UMAP|1.3.2 What You’re Actually Seeing When Samples “Mix” on the UMAP]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.3 Where Treatment Effects Are Actually Preserved (And Why You Can’t See Them on the UMAP)|1.3.3 Where Treatment Effects Are Actually Preserved (And Why You Can’t See Them on the UMAP)]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.4 How We Actually Detect and Measure Treatment Effects|1.3.4 How We Actually Detect and Measure Treatment Effects]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.5 Quantitative Proof That Biology Isn’t Lost|1.3.5 Quantitative Proof That Biology Isn’t Lost]]
		- [[#1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?#1.3.6 What Actual Overcorrection Looks Like (And How to Recognise It)|1.3.6 What Actual Overcorrection Looks Like (And How to Recognise It)]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.4 How Integration Works: The Conceptual Framework|1.4 How Integration Works: The Conceptual Framework]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.5 Evaluating Integration Quality: How to Know It Worked|1.5 Evaluating Integration Quality: How to Know It Worked]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.6 Understanding Clustering: Finding Cell Populations|1.6 Understanding Clustering: Finding Cell Populations]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.7 The Complete Multi-Sample Analysis Workflow|1.7 The Complete Multi-Sample Analysis Workflow]]
	- [[#1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis#1.8 What You’ll Learn in This Tutorial|1.8 What You’ll Learn in This Tutorial]]
- [[#2. Setting Up Your R Environment|2. Setting Up Your R Environment]]
	- [[#2. Setting Up Your R Environment#2.1 Required Packages for Integration and Clustering|2.1 Required Packages for Integration and Clustering]]
	- [[#2. Setting Up Your R Environment#2.2 Loading Libraries and Configuration|2.2 Loading Libraries and Configuration]]
- [[#3. Loading Quality-Controlled Data from Part 2|3. Loading Quality-Controlled Data from Part 2]]
	- [[#3. Loading Quality-Controlled Data from Part 2#3.1 Understanding “batches” in this dataset|3.1 Understanding “batches” in this dataset]]
	- [[#3. Loading Quality-Controlled Data from Part 2#3.2 Loading All QC’d Samples|3.2 Loading All QC’d Samples]]
	- [[#3. Loading Quality-Controlled Data from Part 2#3.3 Quick Quality Check|3.3 Quick Quality Check]]
- [[#4. The Naive Approach: Direct Merge Without Integration|4. The Naive Approach: Direct Merge Without Integration]]
	- [[#4. The Naive Approach: Direct Merge Without Integration#4.1 Merging Samples Directly|4.1 Merging Samples Directly]]
	- [[#4. The Naive Approach: Direct Merge Without Integration#4.2 Visualizing Batch Effects|4.2 Visualizing Batch Effects]]
	- [[#4. The Naive Approach: Direct Merge Without Integration#4.3 Quantifying Batch Effects|4.3 Quantifying Batch Effects]]
- [[#5. Data Integration: Five Methods Compared|5. Data Integration: Five Methods Compared]]
	- [[#5. Data Integration: Five Methods Compared#5.1 Overview of Integration Methods|5.1 Overview of Integration Methods]]
	- [[#5. Data Integration: Five Methods Compared#5.2 Why This Tutorial Focuses on Four R-Based Methods|5.2 Why This Tutorial Focuses on Four R-Based Methods]]
	- [[#5. Data Integration: Five Methods Compared#5.3 Preparing Data for Integration|5.3 Preparing Data for Integration]]
	- [[#5. Data Integration: Five Methods Compared#5.4 Method 1: CCA Integration|5.4 Method 1: CCA Integration]]
	- [[#5. Data Integration: Five Methods Compared#5.5 Method 2: RPCA Integration|5.5 Method 2: RPCA Integration]]
	- [[#5. Data Integration: Five Methods Compared#5.6 Method 3: Harmony Integration|5.6 Method 3: Harmony Integration]]
	- [[#5. Data Integration: Five Methods Compared#5.7 Method 4: FastMNN Integration|5.7 Method 4: FastMNN Integration]]
	- [[#5. Data Integration: Five Methods Compared#5.8 Method 5: scVI Integration (Advanced – Not Covered)|5.8 Method 5: scVI Integration (Advanced – Not Covered)]]
- [[#6. Comparing Integration Methods|6. Comparing Integration Methods]]
	- [[#6. Comparing Integration Methods#6.1 Visual Comparison: UMAPs Side by Side|6.1 Visual Comparison: UMAPs Side by Side]]
	- [[#6. Comparing Integration Methods#6.2 Quantitative Comparison: Mixing and Separation Metrics|6.2 Quantitative Comparison: Mixing and Separation Metrics]]
	- [[#6. Comparing Integration Methods#6.3 Assessing Biological Signal Preservation|6.3 Assessing Biological Signal Preservation]]
	- [[#6. Comparing Integration Methods#6.4 Selecting the Best Integration Method|6.4 Selecting the Best Integration Method]]
- [[#7. Understanding Clustering: From Integrated Data to Cell Populations|7. Understanding Clustering: From Integrated Data to Cell Populations]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.1 The Clustering Algorithm: Graph-Based Community Detection|7.1 The Clustering Algorithm: Graph-Based Community Detection]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.2 The Resolution Parameter: Controlling Cluster Granularity|7.2 The Resolution Parameter: Controlling Cluster Granularity]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.3 Common Misconception: Cluster Number Should Match Expected Cell Types|7.3 Common Misconception: Cluster Number Should Match Expected Cell Types]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.4 Clustering at Multiple Resolutions|7.4 Clustering at Multiple Resolutions]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.5 Visualizing Multi-Resolution Clustering|7.5 Visualizing Multi-Resolution Clustering]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.6 Determining Optimal Resolution|7.6 Determining Optimal Resolution]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.7 Cluster Quality Assessment|7.7 Cluster Quality Assessment]]
	- [[#7. Understanding Clustering: From Integrated Data to Cell Populations#7.8 Final Integrated and Clustered Visualization|7.8 Final Integrated and Clustered Visualization]]
- [[#8. Saving Integrated and Clustered Data|8. Saving Integrated and Clustered Data]]
- [[#9. Understanding Your Integration Results: Interpreting the Plots|9. Understanding Your Integration Results: Interpreting the Plots]]
	- [[#9. Understanding Your Integration Results: Interpreting the Plots#9.1 What Good Integration Looks Like|9.1 What Good Integration Looks Like]]
	- [[#9. Understanding Your Integration Results: Interpreting the Plots#9.2 Interpreting Multi-Resolution Clustering|9.2 Interpreting Multi-Resolution Clustering]]
	- [[#9. Understanding Your Integration Results: Interpreting the Plots#9.3 Common Integration Issues and Diagnoses|9.3 Common Integration Issues and Diagnoses]]
- [[#10. Best Practices for Integration and Clustering|10. Best Practices for Integration and Clustering]]
	- [[#10. Best Practices for Integration and Clustering#10.1 General Principles|10.1 General Principles]]
	- [[#10. Best Practices for Integration and Clustering#10.2 Method Selection Guidelines|10.2 Method Selection Guidelines]]
	- [[#10. Best Practices for Integration and Clustering#10.3 Resolution Selection Strategy|10.3 Resolution Selection Strategy]]
	- [[#10. Best Practices for Integration and Clustering#10.4 Common Pitfalls to Avoid|10.4 Common Pitfalls to Avoid]]
	- [[#10. Best Practices for Integration and Clustering#10.5 When Integration May Not Be Necessary|10.5 When Integration May Not Be Necessary]]
	- [[#10. Best Practices for Integration and Clustering#10.6 Issue-Specific Solutions|10.6 Issue-Specific Solutions]]
- [[#11. Conclusion: From Raw Data to Biological Insights|11. Conclusion: From Raw Data to Biological Insights]]
	- [[#11. Conclusion: From Raw Data to Biological Insights#11.1 Critical Insights|11.1 Critical Insights]]
	- [[#11. Conclusion: From Raw Data to Biological Insights#11.2 What Makes Good Integration|11.2 What Makes Good Integration]]
	- [[#11. Conclusion: From Raw Data to Biological Insights#11.3 Resources for Continued Learning|11.3 Resources for Continued Learning]]
- [[#12. References|12. References]]

## 1. Introduction: Why Integration Matters in Multi-Sample scRNA-seq Analysis

In [Part 1](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-1-from-fastq-to-count-matrix/) and [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/) of this tutorial series, we processed PBMC samples from the [GSE174609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174609) dataset through the complete pipeline: from raw FASTQ files to quality-controlled count matrices. Now we face a critical question: **How do we analyse multiple samples together to identify cell types and compare biological conditions?**

For this tutorial, we’ll work with **8 samples** (*<u>4 Healthy + 4 Post-Treatment</u>*) from the full dataset. This subset keeps the analysis manageable for learning while still demonstrating all integration concepts. The same workflow scales to any number of samples.

The intuitive approach might be to simply merge all samples into one large dataset and proceed with clustering. However, this “naive merge” often produces misleading results where cells cluster by **sample identity** or **processing batch** rather than by **biological cell type**. This technical variation obscures the true biological signal we’re trying to detect.

### 1.1 The Batch Effect Problem in Single-Cell Data

**Batch effects** are systematic technical differences between samples that arise from:
1. **Different processing days**: Samples prepared on different days experience varying environmental conditions (temperature, humidity, reagent age)
2. **Different operators**: Even with standardised protocols, human technique introduces variability
3. **10x chip variation**: Each Chromium chip has slight manufacturing differences
4. **Sequencing runs**: Different flow cells or sequencing dates create systematic biases
5. **Sample preparation timing**: Hours between sample collection and processing affect cell stress responses

These technical variations can be **larger than biological differences** between cell types or treatment conditions, making it impossible to distinguish true biology from technical artefacts.

### 1.2 What is Data Integration?

==**Data integration** (also called ***batch correction*** or ***data harmonisation***) is the computational process of removing technical variation between samples while preserving biological variation==. Think of it as aligning samples so that cells of the same type cluster together regardless of which sample they came from.

**Key principle**: Integration corrects for **who prepared the sample**, not **what the sample represents biologically**.

### 1.3 The Critical Misconception: Will Integration Erase My Treatment Effects?

This is the most common concern when analysing multi-sample experiments with different conditions, and it manifests in a question we see repeatedly from beginners:

> **“After integration, I don’t see my healthy and diseased samples separating on the UMAP anymore. They’re all mixed together. Doesn’t this mean integration overcorrected and removed all the biological differences between my conditions?”**

**Short answer: No. The mixing you see is exactly what successful integration should look like.**

This confusion stems from a fundamental misunderstanding about what UMAP visualisation represents in single-cell analysis versus what PCA plots represent in bulk RNA-seq. Let’s break down why this intuition from bulk RNA-seq doesn’t translate to single-cell data.

#### 1.3.1 The Bulk RNA-seq vs Single-Cell RNA-seq Visualisation Paradigm Shift

**In Bulk RNA-seq PCA plots (what most people learned first):**

When you perform bulk RNA-seq and create a PCA plot, you **expect and want** to see your treatment groups separate into distinct clusters. If your healthy samples cluster on the left and your diseased samples cluster on the right, you celebrate; this visual separation confirms that your treatment has a strong transcriptional effect. The distance between sample clusters on the PCA plot directly reflects the magnitude of treatment-induced changes.

This makes intuitive sense: each sample is one data point representing millions of pooled cells, and the PCA shows you whether samples from different conditions have different overall gene expression profiles.

**In scRNA-seq UMAP plots after integration (the paradigm shift):**

When you perform scRNA-seq integration and create a UMAP plot, you **expect and want** to see your samples **NOT** separating by condition. If cells from healthy donors intermingle completely with cells from diseased patients, you celebrate; this visual mixing confirms that integration successfully removed technical batch effects while preserving biological cell type structure.

This seems counterintuitive at first, but here’s why it makes sense: each data point is now an individual cell, not a sample. You have thousands of cells from each sample, representing multiple different cell types (T cells, B cells, monocytes, etc.). The UMAP is designed to show you **cell type relationships**, not **condition differences**.

**The key insight that resolves the confusion:**

- **Bulk RNA-seq PCA**: Shows sample-level differences (treatment effects)
- **scRNA-seq UMAP after integration**: Shows cell-level similarities (cell type structure)

These are fundamentally different questions being answered by fundamentally different visualisations.

#### 1.3.2 What You’re Actually Seeing When Samples “Mix” on the UMAP

Let’s use a concrete example from your periodontitis study to understand what’s happening:

**Before Integration (the problem):**

Your UMAP shows clear separation, but for the **wrong reasons**:

- All cells from Healthy Donor 1 cluster together (T cells + B cells + monocytes from this donor)
- All cells from Healthy Donor 2 cluster in a different region (T cells + B cells + monocytes from this donor)
- All cells from Post-Treatment Patient 1 cluster in yet another region
- All cells from Post-Treatment Patient 2 cluster separately

This separation is a **technical artefact**, not biology. Why? Because T cells from Healthy Donor 1 should be more similar to T cells from Post-Treatment Patient 1 (same cell type) than to B cells from Healthy Donor 1 (same sample but different cell type). But they’re not clustering that way—they’re clustering by which sample they came from, which indicates batch effects.

**After Integration (the solution):**

Your UMAP now shows mixing:

- All T cells cluster together (from healthy donors + post-treatment patients intermixed)
- All B cells cluster together (from healthy donors + post-treatment patients intermixed)
- All monocytes cluster together (from healthy donors + post-treatment patients intermixed)

This mixing is **biological truth**: cells of the same type ARE more similar to each other than cells of different types, regardless of which donor they came from or which condition that donor represents.

**The critical realisation:** The treatment differences haven’t disappeared—they’re just not visible in the UMAP because UMAP shows cell type structure, not treatment effects.

#### 1.3.3 Where Treatment Effects Are Actually Preserved (And Why You Can’t See Them on the UMAP)

Here’s the resolution to the paradox: **Treatment effects exist in gene expression space, not in UMAP visualisation space.**

**What Integration Actually Does:**

Integration algorithms are sophisticated enough to distinguish between two types of variation in your data:

**Type 1: Sample-level technical variation (what integration REMOVES)**

- Healthy Donor 1 was processed on Monday morning by Technician A
- Healthy Donor 2 was processed on Tuesday afternoon by Technician B
- The technical differences in sample processing create systematic biases
- These biases affect all cell types equally within each sample
- Integration identifies and removes these sample-specific technical effects

**Type 2: Condition-level biological variation (what integration PRESERVES)**

- Post-treatment monocytes express different genes than healthy monocytes
- Post-treatment samples have different proportions of cell types than healthy samples
- Post-treatment CD4+ T cells show different activation states than healthy CD4+ T cells
- These biological differences are consistent across all post-treatment patients
- Integration recognises these as real biology and preserves them

**Why can’t you see treatment effects on the UMAP?**

The UMAP is a 2-dimensional projection of a 20,000-dimensional gene expression space. It’s specifically optimised to show the **largest sources of variation** in your data. In single-cell data, the largest source of variation is **cell type identity** —the differences between T cells, B cells, and monocytes dwarf the differences between healthy and diseased states within any given cell type.

Think of it this way: the difference between a car and a bicycle (cell types) is much larger than the difference between a red car and a blue car (treatment conditions within a cell type). The UMAP emphasises the car vs bicycle distinction, but the red vs blue information is still preserved in the full dataset—you just need different analyses to see it.

#### 1.3.4 How We Actually Detect and Measure Treatment Effects

This is where Part 5 of the tutorial series becomes critical. Treatment effects are detected through **three complementary approaches**, none of which rely on visual UMAP separation:

**Approach 1: Within-Cell-Type Differential Expression**

After identifying cell types, you compare gene expression between conditions within each specific cell type. For example:

- Compare monocytes from healthy donors vs monocytes from post-treatment patients
- This reveals: Post-treatment monocytes express 2-fold higher levels of anti-inflammatory cytokines
- The treatment effect is clear and quantified, even though healthy and diseased monocytes mixed on the UMAP

**Approach 2: Cell Type Composition Changes**

Calculate the proportion of each cell type in healthy vs diseased samples:

- Healthy samples: 60% T cells, 15% B cells, 20% monocytes, 5% other
- Post-treatment samples: 50% T cells, 15% B cells, 30% monocytes, 5% other
- Treatment effect: 50% increase in monocyte proportion
- Again, this analysis doesn’t require UMAP separation—it’s a simple counting exercise

**Approach 3: Cell State Distribution Shifts**

Within a cell type, examine how cells distribute across different activation states:

- Healthy CD4+ T cells: 70% resting, 20% activated, 10% memory
- Post-treatment CD4+ T cells: 50% resting, 30% activated, 20% memory
- Treatment effect: Shift toward activated and memory states
- These shifts are detected through gene expression signatures, not UMAP positions

**The key insight:** All three approaches reveal treatment effects that are completely invisible on the UMAP but perfectly preserved in the underlying data.

#### 1.3.5 Quantitative Proof That Biology Isn’t Lost

Remember when we calculated integration quality metrics (*come later in Steps [[#6.2 Quantitative Comparison Mixing and Separation Metrics|12]] and [[#6.3 Assessing Biological Signal Preservation|13]]*)? We measured **two independent scores**:

**Score 1: Sample Mixing (Technical Correction)**

- Measures how well cells from different samples intermingle within neighbourhoods
- High score = batch effects successfully removed
- This is what makes the UMAP look “mixed”

**Score 2: Biological Preservation (Signal Retention)**

- Measures whether conditions remain distinguishable using appropriate metrics
- High score = biological differences retained despite visual mixing
- This proves the treatment signal wasn’t destroyed

**Your integration method was chosen because it had BOTH high mixing AND high biological preservation.** If integration had overcorrected, you would see:

- High mixing score (samples mix) ✓
- **Low biological preservation score** (conditions indistinguishable) ✗
- But you don’t see this—your biological preservation score is also high

This quantitative evidence confirms that the visual mixing on the UMAP does not indicate lost biology.

#### 1.3.6 What Actual Overcorrection Looks Like (And How to Recognise It)

True overcorrection has specific warning signs that go beyond just visual UMAP mixing:

**Warning Sign 1: Known Marker Genes Lose Specificity**

- CD3D (T cell marker) becomes equally expressed in B cells and monocytes
- CD14 (monocyte marker) shows up strongly in T cells
- Cell-type-defining genes become ambiguous across all clusters

**Warning Sign 2: Biological Preservation Score Drops**

- Quantitative metrics show conditions are now indistinguishable
- Treatment groups have become statistically inseparable
- Known disease signatures disappear from the data

**Warning Sign 3: Cell Type Proportions Become Identical**

- Healthy samples: 60% T cells, 15% B cells, 25% monocytes
- Diseased samples: 60% T cells, 15% B cells, 25% monocytes
- All biological differences in cell composition artificially equalised

**Warning Sign 4: Within-Cell-Type Differential Expression Fails**

- Testing for treatment effects within monocytes returns zero differentially expressed genes
- Known disease markers (validated by qPCR or flow cytometry) no longer show differences
- Literature-validated biology cannot be recapitulated

**Warning Sign 5: Results Contradict Orthogonal Validation**

- Flow cytometry shows 50% more activated T cells in disease
- Your integrated data shows no activation difference
- Integration has removed real biology

**The critical distinction:** Visual UMAP mixing alone is NOT evidence of overcorrection. You need to see multiple warning signs, particularly loss of known biology and failed differential expression, before concluding that integration was too aggressive.

### 1.4 How Integration Works: The Conceptual Framework

Modern integration methods use **anchoring strategies** to align datasets:

1. **Identify anchor cells**: Find cells from different samples that are clearly the same cell type or cluster (e.g., T cells in Sample 1 match T cells in Sample 2)
2. **Learn technical differences**: Calculate how gene expression differs between anchored cells from different samples
3. **Apply corrections**: Adjust all cells to remove sample-specific technical variation
4. **Preserve biology**: Retain variation that’s consistent across samples (true biological differences)

Different integration methods use variations of this strategy with different assumptions and strengths.

### 1.5 Evaluating Integration Quality: How to Know It Worked

Good integration should produce specific visual and quantitative outcomes:

**Visual indicators:**

1. **UMAP mixing**: Cells from different samples intermingle within each cluster
2. **Preserved cell types**: Distinct clusters corresponding to known biological cell types remain separated
3. **Condition-specific populations**: Rare or condition-specific cell types are still detectable
4. **No technical clusters**: No clusters driven purely by sample identity

**Quantitative metrics:**

1. **Mixing metric**: Measures how well samples mix within neighbourhoods
2. **Silhouette score**: Assesses cluster separation
3. **kBET (k-nearest neighbour Batch Effect Test)**: Tests for batch effects in local neighbourhoods
4. **Biological conservation**: Marker genes remain differentially expressed

**Red flags indicating over-integration:**

- All samples perfectly mixed even in rare cell populations
- Loss of known marker gene expression
- Disappearance of expected biological variation
- Creation of artificial intermediate cell types

**Red flags indicating under-integration:**

- Samples form distinct clusters instead of mixing
- Same cell-type splits across multiple clusters by sample
- Technical variation dominates biological signal

### 1.6 Understanding Clustering: Finding Cell Populations

After integration removes technical variation, **clustering** groups cells based on transcriptional similarity. This allows us to:

1. **Identify cell types**: T cells, B cells, monocytes, etc.
2. **Discover cell states**: Activated vs resting, proliferating vs quiescent
3. **Find rare populations**: Dendritic cells, stem cells, transitional states
4. **Enable downstream analysis**: Differential expression requires knowing which cells to compare

**Clustering workflow:**

1. **Dimensionality reduction**: Compress 20,000 genes into ~30 principal components capturing major variation
2. **Nearest neighbour graph**: Connect each cell to its most similar neighbours
3. **Community detection**: Find groups of densely connected cells (clusters)
4. **Visualisation**: Project onto 2D UMAP for interpretation

The **resolution parameter** controls clustering granularity:

- **Low resolution (0.2-0.5)**: Broad cell types (T cells, B cells, monocytes)
- **Medium resolution (0.6-1.0)**: Cell subtypes (CD4+ T cells, CD8+ T cells, naive vs memory)
- **High resolution (1.2-2.0)**: Fine-grained states (Th1 vs Th2 vs Th17 CD4+ T cells)

### 1.7 The Complete Multi-Sample Analysis Workflow

```
8 QC'd Samples (from Part 2: 4 Healthy + 4 Post-Treatment)
    ↓
Attempt 1: Naive Merge (no integration)
    ↓ → Batch effects dominate
    ↓ → Cells cluster by sample, not biology
    ↓
Attempt 2: Four Integration Methods (R-based)
    ↓
    ├─→ CCA Integration (canonical correlation)
    ├─→ RPCA Integration (reciprocal PCA)
    ├─→ Harmony (probabilistic, fast)
    └─→ FastMNN (mutual nearest neighbors)
        ↓
Compare Integration Results
    ↓ (mixing scores, biological preservation)
    ↓
Select Best Method
    ↓
Clustering at Multiple Resolutions
    ↓
Evaluate Cluster Quality
    ↓
Cell Type Annotation
    ↓
Differential Expression
```

### 1.8 What You’ll Learn in This Tutorial

By the end of Part 3, you’ll be able to:

✓ **Understand** why naive merging fails and integration is necessary
✓ **Apply** four different integration methods to your data
✓ **Compare** integration results visually and quantitatively
✓ **Select** the optimal method for your specific dataset
✓ **Perform** clustering at appropriate resolutions
✓ **Evaluate** whether clusters represent true biology or artefacts
✓ **Prepare** data for cell type annotation and downstream analysis

Let’s begin by setting up our environment and loading the quality-controlled data from [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/).

## 2. Setting Up Your R Environment

We’ll build on the environment established in [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/), adding packages specifically for integration methods.

### 2.1 Required Packages for Integration and Clustering

```R
#-----------------------------------------------
# Installation: Integration-specific packages
#-----------------------------------------------

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core packages (already installed in Part 2)
# Seurat, ggplot2, dplyr, patchwork

# Install Bioconductor packages for advanced integration
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "batchelor",  # FastMNN integration (used by Seurat's IntegrateLayers)
  "scran"       # Additional normalisation methods
))

# Install Harmony for integration
install.packages("harmony")

# Install additional visualisation and utility packages
install.packages(c(
  "remotes",         # Install R packages stored in GitHub
  "ggrepel",         # Non-overlapping text labels
  "RColorBrewer",    # Color palettes
  "viridis"# Perceptually uniform color scales
))

# Wrapper functions for integration methods
remotes::install_github('satijalab/seurat-wrappers')

# Install future for parallel processing (used by RPCA integration)
install.packages("future")

# Install FNN for k-nearest neighbour calculations (used in quality metrics)
install.packages("FNN")

# Install cluster for clustering quality metrics (silhouette scores)
install.packages("cluster")

# Install reshape2 for data reshaping (used in cluster quality assessment)
install.packages("reshape2")
```

### 2.2 Loading Libraries and Configuration

```R
#-----------------------------------------------
# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Integration methods
library(harmony)
library(batchelor)
library(SeuratWrappers)

# Visualization and data manipulation
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Quality metrics and utilities
library(FNN)              # K-nearest neighbor calculations for mixing metrics
library(cluster)          # Silhouette scores for clustering quality
library(reshape2)         # Data reshaping for visualization

# Parallel processing
library(future)
options(future.globals.maxSize = 20 * 1024^3)  # Increase to 20GB for large datasets

# Set working directory (adjust to your path)
setwd("~/GSE174609_scRNA/integration_analysis")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/integration_comparison", showWarnings = FALSE)
dir.create("plots/clustering", showWarnings = FALSE)
dir.create("integrated_data", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))

cat("Environment setup complete\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
# Seurat version: 5.3.1
```

## 3. Loading Quality-Controlled Data from Part 2

In [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/), we performed comprehensive quality control on samples from [GSE174609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174609). Each sample was independently filtered for empty droplets, doublets, low-quality cells, and uninformative genes. Now we’ll load these clean datasets for integration.

For this tutorial, we use **8 samples** (4 Healthy + 4 Post-Treatment) to keep the analysis manageable while demonstrating all key concepts.

### 3.1 Understanding “batches” in this dataset

In scRNA-seq, each individual sample typically represents its own technical batch due to separate processing. Even though samples may be processed on the same day, slight variations in handling time, reagent temperature, cell dissociation efficiency, and other technical factors create sample-specific effects. Therefore:

- **Technical batch = individual sample** (8 batches in our analysis)
- **Biological replicate = patient/donor** (different individuals within each condition)

Integration corrects for the technical differences between these 8 samples while preserving the biological differences between the 2 conditions (Healthy vs Post-treatment).

**Expected outcomes after integration:**

- Samples should mix within each cell type
- Treatment effects should remain visible in:
- Cell type proportions
- Gene expression within cell types
- Cell state distributions

### 3.2 Loading All QC’d Samples

```R
#-----------------------------------------------
# STEP 2: Load QC-filtered samples from Part 2
#-----------------------------------------------

# Define sample metadata
sample_metadata <- data.frame(
  sample_id = c(
    "Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4",
    "Post_Patient_1", "Post_Patient_2", "Post_Patient_3", "Post_Patient_4"
  ),
  condition = c(
    rep("Healthy", 4),
    rep("Post_Treatment", 4)
  ),
  patient_id = c(
    "Donor_1", "Donor_2", "Donor_3", "Donor_4",
    "Patient_1", "Patient_2", "Patient_3", "Patient_4"
  ),
  stringsAsFactors = FALSE
)

# Load QC-filtered Seurat objects
qc_data_path <- "~/GSE174609_scRNA/QC_analysis/filtered_data"

seurat_list <- lapply(sample_metadata$sample_id, function(sample_id) {
  obj <- readRDS(file.path(qc_data_path, paste0(sample_id, "_qc_filtered.rds")))

  # Ensure metadata is consistent
  obj$sample_id <- sample_id
  obj$condition <- sample_metadata$condition[sample_metadata$sample_id == sample_id]
  obj$patient_id <- sample_metadata$patient_id[sample_metadata$sample_id == sample_id]

  return(obj)
})

names(seurat_list) <- sample_metadata$sample_id

# Ensure all samples have the same genes (critical for FastMNN integration)
all_genes <- lapply(seurat_list, rownames)
common_genes <- Reduce(intersect, all_genes)

# Subset all samples to common genes
seurat_list <- lapply(seurat_list, function(obj) {
  obj[common_genes, ]
})

# Report dimensions
cat("\nLoaded", length(seurat_list), "QC-filtered samples:\n")
for (sample_id in names(seurat_list)) {
  cat(sprintf("  %s: %d cells × %d genes\n",
              sample_id,
              ncol(seurat_list[[sample_id]]),
              nrow(seurat_list[[sample_id]])))
}

cat("\nTotal cells across all samples:",
    sum(sapply(seurat_list, ncol)), "\n")

#   Healthy_1: 8018 cells × 18861 genes
#   Healthy_2: 8249 cells × 18861 genes
#   Healthy_3: 9057 cells × 18861 genes
#   Healthy_4: 10175 cells × 18861 genes
#   Post_Patient_1: 8132 cells × 18861 genes
#   Post_Patient_2: 8465 cells × 18861 genes
#   Post_Patient_3: 9339 cells × 18861 genes
#   Post_Patient_4: 11214 cells × 18861 genes
```

> **Important**: These Seurat objects were already normalized in [Part 2](https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/) using `NormalizeData()` and have variable features identified with `FindVariableFeatures()`. We’ll use this normalized data for integration.

### 3.3 Quick Quality Check

Before proceeding, let’s verify our data quality across all samples. This is a good checkpoint before starting integration.

## 4. The Naive Approach: Direct Merge Without Integration

Before exploring integration methods, let’s see what happens when we simply merge samples without correction. This demonstrates **why integration is necessary**.

### 4.1 Merging Samples Directly

```R
#-----------------------------------------------
# STEP 3: Naive merge without integration
#-----------------------------------------------

# Merge all samples
merged_naive <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE174609_Naive_Merge"
)

cat("Merged dataset:", ncol(merged_naive), "cells ×", nrow(merged_naive), "genes\n")
# Merged dataset: 72649 cells × 18861 genes

# Standard Seurat workflow
merged_naive <- NormalizeData(merged_naive, verbose = FALSE)
merged_naive <- FindVariableFeatures(merged_naive, nfeatures = 2000, verbose = FALSE)
merged_naive <- ScaleData(merged_naive, verbose = FALSE)
merged_naive <- RunPCA(merged_naive, npcs = 50, verbose = FALSE)
merged_naive <- RunUMAP(merged_naive, dims = 1:30, reduction = "pca", verbose = FALSE)
merged_naive <- FindNeighbors(merged_naive, dims = 1:30, verbose = FALSE)
merged_naive <- FindClusters(merged_naive, resolution = 0.6, verbose = FALSE)

cat("Clustering complete:", length(unique(merged_naive$seurat_clusters)), "clusters identified\n")
# Clustering complete: 19 clusters identified
```

### 4.2 Visualizing Batch Effects

```R
#-----------------------------------------------
# STEP 4: Visualize batch effects in naive merge
#-----------------------------------------------

# Define color palettes for visualization
# 8 samples need 8 distinct, visible colors
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",  # Healthy 1-4: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#999999"# Post_Patient 1-4: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# 2 conditions need 2 distinct colors
condition_colors <- c(
  "Healthy"= "#2E86AB",       # Blue
  "Post_Treatment"= "#F18F01"# Orange
)

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

# Combine plots
combined_naive <- (p1_naive | p2_naive) / (p3_naive | p4_naive)
ggsave("plots/01_naive_merge_batch_effects.png", combined_naive,
       width = 16, height = 12, dpi = 300)
```

> **Understanding Batch Effects in Large Datasets**: With 72,649 cells in our example, the UMAP plot becomes very dense with overlapping points, making it difficult to visually assess batch effects by color alone. While smaller datasets (<10,000 cells) show clear visual separation when samples cluster by batch rather than biology, large datasets require **quantitative metrics** (mixing scores, silhouette scores) to properly evaluate batch effects. The split view (bottom right panel) and the quantification section below provide clearer evidence of batch-driven clustering.
>
> **What Batch Effects Look Like**: In a naive merge, cells often cluster by **sample identity** rather than **cell type**. This means Healthy\_1 and Healthy\_2 may form distinct clusters even though they contain the same cell types. Integration corrects this by aligning samples based on shared biological variation while removing technical variation.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/01_naive_merge_batch_effects-scaled.png?resize=2048%2C1536&ssl=1)

**Here is another dataset with far fewer cells. The UMAPs from the naïve merge clearly reveal the batch effects discussed above.**

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/01_naive_merge_batch_effects-1-scaled.png?resize=2048%2C1536&ssl=1)

### 4.3 Quantifying Batch Effects

```R
#-----------------------------------------------
# STEP 5: Quantify batch effects with mixing metrics
#-----------------------------------------------

# Calculate per-cluster sample composition
cluster_composition <- table(merged_naive$seurat_clusters, merged_naive$sample_id)
cluster_composition_pct <- prop.table(cluster_composition, margin = 1) * 100

# Find clusters dominated by single samples (>50% from one sample = batch-driven)
dominant_sample_clusters <- apply(cluster_composition_pct, 1, max) > 50
n_batch_clusters <- sum(dominant_sample_clusters)

cat("\nBatch effect assessment:\n")
cat("  Clusters dominated by single sample (>50%):", n_batch_clusters,
    "/", nrow(cluster_composition), "\n")
# Clusters dominated by single sample (>50%): 1 / 19

if (n_batch_clusters > 0) {
  cat("  ⚠  Naive merge shows significant batch effects\n")
  cat("  → Integration is necessary\n")
}
# → Integration is necessary
```

## 5. Data Integration: Five Methods Compared

Now we’ll apply four different integration methods to remove batch effects while preserving biological variation. Each method has different strengths and is suited to different scenarios.

### 5.1 Overview of Integration Methods

**1\. CCA Integration (Canonical Correlation Analysis)**

- **Algorithm**: Identifies shared correlation structures across datasets
- **Strengths**:
- Excellent for diverse cell types
- Handles datasets with different cell type compositions
- Robust to varying cell numbers per sample
- **Best for**: Experiments where samples may have different cell type proportions
- **Computational cost**: Medium
- **When to use**: Default choice for most experiments

**2\. RPCA Integration (Reciprocal PCA)**

- **Algorithm**: Uses reciprocal PCA to identify shared variation
- **Strengths**:
	- *Faster than CCA*
	- *Works well with large datasets (>100,000 cells)*
	- *Less memory-intensive*
- **Best for**: Large-scale studies, similar cell type compositions across samples
- **Computational cost**: Low-Medium
- **When to use**: When computational efficiency is critical, or datasets are very large

**3\. Harmony**

- **Algorithm**: Iteratively adjusts cell embeddings to remove batch effects
- **Strengths**:
	- *Very fast*
	- *Preserves global structure well*
	- *Simple single-step process*
	- *Scales to very large datasets*
- **Best for**: Large datasets, soft batch correction, preserving broad structure
- **Computational cost**: Low
- **When to use**: First-pass integration, large datasets, when speed matters

**4\. FastMNN (Fast Mutual Nearest Neighbours)**

- **Algorithm**: Finds mutual nearest neighbours across batches and corrects in low-dimensional space
- **Strengths**:
	- *Preserves biological variation well*
	- *Handles complex batch structures*
	- *Less aggressive correction (conservative)*
- **Best for**: Datasets with subtle biological differences, hierarchical batch structures
- **Computational cost**: Medium
- **When to use**: When biological signal is subtle, and you want conservative correction

**5\. scVI (Single-Cell Variational Inference)**

- **Algorithm**: Deep learning-based method using variational autoencoders
- **Strengths**:
	- *State-of-the-art performance on complex datasets*
	- *Can integrate millions of cells*
	- *Handles very large batch effects*
	- *Learns data-specific corrections*
- **Best for**: Very large datasets, complex batch structures, when computational resources allow
- **Computational cost**: High (GPU recommended)
- **When to use**: When R-based methods show insufficient correction and computational resources are available

### 5.2 Why This Tutorial Focuses on Four R-Based Methods

The four R-based methods (CCA, RPCA, Harmony, FastMNN) we cover are:
	*- Easy to install through CRAN/Bioconductor*
	*- Work reliably across different systems*
	*- Often perform comparably to scVI for most datasets*
	*- Sufficient for high-quality integration*
	*- No GPU required*

These methods will meet the needs of >95% of scRNA-seq integration tasks. We provide detailed scVI information later for advanced users who need it.

### 5.3 Preparing Data for Integration

All integration methods in Seurat 5 use a unified workflow through `IntegrateLayers()`. First, we need to prepare our merged object:

```R
#-----------------------------------------------
# STEP 6: Prepare data for integration
#-----------------------------------------------

# Merge all samples with split layers (required for Seurat 5 integration)
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE174609_Integration"
)

# Add all metadata
# Extract sample_id from cell names (added via add.cell.ids)
# Cell names format: "SampleID_BARCODE-1"
cell_names <- colnames(merged_seurat)
merged_seurat$sample_id <- gsub("_[ACGT].*$", "", cell_names)  # Remove barcode, keep sample_id

# Add other metadata by matching sample_id
merged_seurat$condition <- sample_metadata$condition[match(merged_seurat$sample_id, sample_metadata$sample_id)]
merged_seurat$patient_id <- sample_metadata$patient_id[match(merged_seurat$sample_id, sample_metadata$sample_id)]

# Check if layers are already split
current_layers <- Layers(merged_seurat[["RNA"]])
if (length(current_layers) > 1) {
  cat("Layers already split by merge operation (", length(current_layers), " layers)\n", sep = "")
} else{
  # Split layers by sample if not already split
  cat("Splitting layers by sample\n")
  merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$sample_id)
}

# Normalize and find variable features on split layers
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 2000, verbose = FALSE)

# Scale data and run PCA on split layers
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE)

cat("Data prepared for integration\n")
cat("Layers:", length(Layers(merged_seurat, search = "data")), "samples\n")
# Layers: 8 samples
```

### 5.4 Method 1: CCA Integration

```R
#-----------------------------------------------
# STEP 7: CCA Integration
#-----------------------------------------------

# Integrate using CCA
integrated_cca <- IntegrateLayers(
  object = merged_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

# Standard downstream workflow
integrated_cca <- FindNeighbors(integrated_cca, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
integrated_cca <- FindClusters(integrated_cca, resolution = 0.6, verbose = FALSE)
integrated_cca <- RunUMAP(integrated_cca, reduction = "integrated.cca", dims = 1:30,
                          reduction.name = "umap.cca", verbose = FALSE)

cat("CCA integration complete:", length(unique(integrated_cca$seurat_clusters)), "clusters\n")
# CCA integration complete: 23 clusters
```

### 5.5 Method 2: RPCA Integration

RPCA is computationally efficient and works well with large datasets.

> **Note on Parallel Processing**: RPCA integration uses the `future` package for parallel processing. We configured `future.globals.maxSize = 20GB` in Step 1 to handle large object transfers. If you still encounter memory errors, increase this value or reduce the dataset size.

```R
#-----------------------------------------------
# STEP 8: RPCA Integration
#-----------------------------------------------

# Integrate using RPCA
integrated_rpca <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

# Downstream workflow
integrated_rpca <- FindNeighbors(integrated_rpca, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
integrated_rpca <- FindClusters(integrated_rpca, resolution = 0.6, verbose = FALSE)
integrated_rpca <- RunUMAP(integrated_rpca, reduction = "integrated.rpca", dims = 1:30,
                           reduction.name = "umap.rpca", verbose = FALSE)

cat("RPCA integration complete:", length(unique(integrated_rpca$seurat_clusters)), "clusters\n")
# RPCA integration complete: 22 clusters
```

### 5.6 Method 3: Harmony Integration

Harmony is fast and effective, working directly on the PCA embedding rather than requiring layer-based operations.

```R
#-----------------------------------------------
# STEP 9: Harmony Integration
#-----------------------------------------------

# Harmony works directly on PCA embedding
# First, join layers and run standard workflow
merged_harmony <- JoinLayers(merged_seurat)
merged_harmony <- RunPCA(merged_harmony, npcs = 50, verbose = FALSE)

# Run Harmony (corrects batch effects on PCA embedding)
integrated_harmony <- RunHarmony(merged_harmony, "sample_id")

# Downstream workflow
integrated_harmony <- FindNeighbors(integrated_harmony, reduction = "harmony", dims = 1:30, verbose = FALSE)
integrated_harmony <- FindClusters(integrated_harmony, resolution = 0.6, verbose = FALSE)
integrated_harmony <- RunUMAP(integrated_harmony, reduction = "harmony", dims = 1:30,
                                reduction.name = "umap.harmony", verbose = FALSE)

cat("Harmony integration complete:", length(unique(integrated_harmony$seurat_clusters)), "clusters\n")
# Harmony integration complete: 21 clusters
```

### 5.7 Method 4: FastMNN Integration

FastMNN uses mutual nearest neighbors to align datasets. It works directly on the split layers in the merged Seurat object using `IntegrateLayers()`.

```R
#-----------------------------------------------
# STEP 10: FastMNN Integration
#-----------------------------------------------

# Integrate using FastMNN (works on split layers)
integrated_fastmnn <- IntegrateLayers(
  object = merged_seurat,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

cat("FastMNN integration complete\n")

# Downstream workflow
integrated_fastmnn <- FindNeighbors(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
integrated_fastmnn <- FindClusters(integrated_fastmnn, resolution = 0.6, verbose = FALSE)
integrated_fastmnn <- RunUMAP(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30,
                              reduction.name = "umap.mnn", verbose = FALSE)

cat("Clustering complete:", length(unique(integrated_fastmnn$seurat_clusters)), "clusters\n")
# FastMNN integration complete: 20 clusters
```

### 5.8 Method 5: scVI Integration (Advanced – Not Covered)

scVI (single-cell Variational Inference) is a powerful deep learning-based integration method that uses variational autoencoders to learn a low-dimensional representation of scRNA-seq data. While it often achieves excellent integration quality, **we skip this method in this tutorial due to significant setup challenges**.

**Why scVI is Powerful:**

- Uses deep learning to model complex batch effects
- Can integrate very large datasets (millions of cells)
- Often produces superior integration compared to linear methods
- Can model both technical and biological variation simultaneously

**Why We’re Skipping It:**

- **Python dependency**: Requires Python 3.8+ and PyTorch installation
- **Complex environment**: Setting up scvi-tools and its dependencies is error-prone
- **Version conflicts**: Frequent incompatibilities between PyTorch, CUDA, and scvi-tools versions
- **Resource intensive**: Requires significant computational resources (GPU recommended)
- **Debugging difficulty**: Python-R interface issues through reticulate are hard to troubleshoot
- **Time investment**: Can take hours to set up correctly, even for experienced users

**Our Recommendation:**
The four R-based methods we covered (CCA, RPCA, Harmony, FastMNN) are:

- **Easier to install** – pure R packages through CRAN/Bioconductor
- **Faster to run** – no deep learning training required
- **Equally effective** – often perform comparably to scVI for most datasets
- **More reproducible** – fewer version dependency issues

**If You Want to Try scVI Anyway:**
For those comfortable with Python environments, see:

- scVI documentation: https://docs.scvi-tools.org/
- Seurat-scVI integration: https://satijalab.org/seurat/articles/seurat5\_integration\_bridge

For this tutorial, **we proceed with our 4 reliable R-based methods**, which are more than sufficient for high-quality integration.

## 6. Comparing Integration Methods

Now that we’ve applied four different integration methods, let’s systematically compare their results to determine which works best for our dataset.

**Methods we’re comparing:**

1. **Naive Merge** (no integration – baseline)
2. **CCA Integration** (Canonical Correlation Analysis)
3. **RPCA Integration** (Reciprocal PCA)
4. **Harmony Integration** (Fast, probabilistic)
5. **FastMNN Integration** (Mutual nearest neighbors)

### 6.1 Visual Comparison: UMAPs Side by Side

```R
#-----------------------------------------------
# STEP 11: Generate comparison UMAPs
#-----------------------------------------------

cat("\n=== Generating Comparison Plots ===\n")

# Create a function to make standardized UMAP plots
make_comparison_plot <- function(seurat_obj, reduction_name, title, group_by = "sample_id", colors = NULL) {
  p <- DimPlot(seurat_obj, reduction = reduction_name, group.by = group_by, pt.size = 0.05) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"))

  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  return(p)
}

# UMAPs colored by sample (assesses mixing)
p_sample_naive <- make_comparison_plot(merged_naive, "umap", "Naive Merge",
                                       colors = sample_colors)
p_sample_cca <- make_comparison_plot(integrated_cca, "umap.cca", "CCA",
                                     colors = sample_colors)
p_sample_rpca <- make_comparison_plot(integrated_rpca, "umap.rpca", "RPCA",
                                      colors = sample_colors)
p_sample_harmony <- make_comparison_plot(integrated_harmony, "umap.harmony", "Harmony",
                                         colors = sample_colors)
p_sample_mnn <- make_comparison_plot(integrated_fastmnn, "umap.mnn", "FastMNN",
                                     colors = sample_colors)

# Combine sample-colored UMAPs
combined_samples <- (p_sample_naive | p_sample_cca | p_sample_rpca) /
                    (p_sample_harmony | p_sample_mnn | plot_spacer())

ggsave("plots/integration_comparison/02_integration_by_sample.png",
       combined_samples, width = 16, height = 12, dpi = 300)

# UMAPs colored by condition (assesses biological preservation)
p_cond_naive <- make_comparison_plot(merged_naive, "umap", "Naive Merge", "condition",
                                     colors = condition_colors)
p_cond_cca <- make_comparison_plot(integrated_cca, "umap.cca", "CCA", "condition",
                                   colors = condition_colors)
p_cond_rpca <- make_comparison_plot(integrated_rpca, "umap.rpca", "RPCA", "condition",
                                    colors = condition_colors)
p_cond_harmony <- make_comparison_plot(integrated_harmony, "umap.harmony", "Harmony", "condition",
                                       colors = condition_colors)
p_cond_mnn <- make_comparison_plot(integrated_fastmnn, "umap.mnn", "FastMNN", "condition",
                                   colors = condition_colors)

# Combine condition-colored UMAPs
combined_conditions <- (p_cond_naive | p_cond_cca | p_cond_rpca) /
                       (p_cond_harmony | p_cond_mnn | plot_spacer())

ggsave("plots/integration_comparison/03_integration_by_condition.png",
       combined_conditions, width = 16, height = 12, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/scRNAseq_integration_comparison-e1766290127773.png?w=1480&ssl=1)

Here is another dataset with far fewer cells, which better demonstrates the effects of integration visually. The UMAP projections compare different integration methods. In contrast to the naïve merge, the integrated datasets show that cells from different samples are well mixed across clusters, with no cluster dominated by cells from a single sample. This indicates successful removal of batch effects associated with sample identity.

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/scRNAseq_integration_comparison_better_viz-e1766501285787.png?w=1480&ssl=1)

### 6.2 Quantitative Comparison: Mixing and Separation Metrics

To objectively compare integration methods, we calculate **mixing scores** that quantify how well cells from different samples mix within biological clusters. Higher scores indicate better integration.

> **How the mixing metric works**: For each cell, we find its k=50 nearest neighbors in the UMAP space and calculate the diversity of samples among those neighbors using the inverse Simpson’s Index. If all neighbours are from different samples (perfect mixing), the score is high. If all neighbours are from the same sample (poor integration), the score is low.

```R
#-----------------------------------------------
# STEP 12: Calculate integration quality metrics
#-----------------------------------------------

# Function to calculate mixing metric (local inverse Simpson's Index)
# This measures sample diversity in each cell's k-nearest neighborhood
calculate_mixing_metric <- function(seurat_obj, group_by = "sample_id", reduction = "umap", k = 50) {
  # Get embedding
  embedding <- Embeddings(seurat_obj, reduction = reduction)

  # Calculate k-nearest neighbors using FNN package
  nn_result <- FNN::get.knn(embedding[, 1:2], k = k)

  # For each cell, calculate diversity of samples in its neighborhood
  group_labels <- seurat_obj@meta.data[[group_by]]
  mixing_scores <- sapply(1:nrow(nn_result$nn.index), function(i) {
    neighbors <- nn_result$nn.index[i, ]
    neighbor_groups <- group_labels[neighbors]
    props <- table(neighbor_groups) / length(neighbor_groups)
    # Inverse Simpson's Index (higher = more mixed)
    1 / sum(props^2)
  })

  return(mean(mixing_scores))
}

# Calculate for all methods
methods_list <- list(
  "Naive"= list(obj = merged_naive, reduction = "umap"),
  "CCA"= list(obj = integrated_cca, reduction = "umap.cca"),
  "RPCA"= list(obj = integrated_rpca, reduction = "umap.rpca"),
  "Harmony"= list(obj = integrated_harmony, reduction = "umap.harmony"),
  "FastMNN"= list(obj = integrated_fastmnn, reduction = "umap.mnn")
)

# Calculate mixing metrics
mixing_results <- data.frame(
  method = names(methods_list),
  mixing_score = sapply(methods_list, function(x) {
    calculate_mixing_metric(x$obj, group_by = "sample_id", reduction = x$reduction)
  }),
  n_clusters = sapply(methods_list, function(x) {
    length(unique(x$obj$seurat_clusters))
  })
)

# Display results
cat("\nIntegration Quality Metrics:\n")
cat("(Higher mixing score = better sample mixing)\n\n")
print(mixing_results)

# Save metrics
write.csv(mixing_results, "metadata/integration_comparison_metrics.csv", row.names = FALSE)

# Visualize mixing scores
p_mixing <- ggplot(mixing_results, aes(x = reorder(method, -mixing_score), y = mixing_score, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mixing_score, 2)), vjust = -0.5) +
  labs(title = "Integration Method Comparison: Sample Mixing",
       subtitle = "Higher score = better integration (samples mix within cell types)",
       x = "Integration Method",
       y = "Mixing Score (Local Inverse Simpson's Index)") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

ggsave("plots/integration_comparison/04_mixing_scores.png", p_mixing,
       width = 10, height = 6, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/04_mixing_scores-scaled.png?resize=2048%2C1229&ssl=1)

### 6.3 Assessing Biological Signal Preservation

A good integration method should **remove technical variation** (samples mix) while **preserving biological variation** (conditions remain distinguishable).

**How This Metric Works (For Beginners):**

Integration is a balancing act. We want to remove batch effects (technical variation) but keep biological differences (biological variation). This metric helps us check if we’ve struck the right balance.

**The Measurement Process:**

1. **Find condition centers**: Calculate the average UMAP position for all Healthy cells and all Post\_Treatment cells
2. **Measure separation**: Calculate the distance between these two centers
3. **Interpret the distance**:
- **High distance** = Conditions are well-separated (biological signal preserved) ✓
- **Low distance** = Conditions overlap (might have over-integrated, losing biology)

**What We’re Looking For:**

- **High mixing score** (from previous section) → Samples within each condition mix well
- **Moderate-to-high separation** (this section) → Conditions remain distinguishable

**Examples:**

- **Good integration**: Healthy cells mix together (mixing score high), BUT Healthy cells and Post\_Treatment cells remain in distinct regions (separation high)
- **Over-integration**: Everything mixes together, losing biological differences (mixing high, separation low)
- **Under-integration**: Samples stay separated, batch effects remain (mixing low, separation high)
```R
#-----------------------------------------------
# STEP 13: Assess biological signal preservation
#-----------------------------------------------

# Function to test if conditions are still distinguishable after integration
# Uses UMAP distances between conditions
assess_condition_separation <- function(seurat_obj, reduction) {
  umap_coords <- seurat_obj[[reduction]]@cell.embeddings[, 1:2]
  conditions <- seurat_obj$condition

  # Calculate mean UMAP position for each condition
  condition_centers <- aggregate(umap_coords, by = list(conditions), FUN = mean)

  # Calculate pairwise distances between condition centers
  dist_matrix <- dist(condition_centers[, -1])
  mean_dist <- mean(dist_matrix)

  return(mean_dist)
}

# Calculate for all methods
separation_results <- data.frame(
  method = names(methods_list),
  condition_separation = sapply(methods_list, function(x) {
    assess_condition_separation(x$obj, x$reduction)
  })
)

mixing_results$condition_separation <- separation_results$condition_separation

cat("\nBiological Signal Preservation:\n")
cat("(Higher separation = conditions remain distinguishable)\n\n")
print(mixing_results[, c("method", "mixing_score", "condition_separation", "n_clusters")])

# Ideal integration: High mixing + Moderate separation
# Over-integration: High mixing + Low separation (conditions lost)
# Under-integration: Low mixing + High separation (batches not corrected)

# Visualize the trade-off
p_tradeoff <- ggplot(mixing_results, aes(x = mixing_score, y = condition_separation, label = method)) +
  geom_point(size = 4, aes(color = method)) +
  geom_text_repel(size = 4, box.padding = 0.5) +
  labs(title = "Integration Quality: Mixing vs Biological Preservation",
       subtitle = "Ideal: Upper right (high mixing + preserved biology)",
       x = "Sample Mixing Score (higher = better integration)",
       y = "Condition Separation (higher = preserved biology)") +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

ggsave("plots/integration_comparison/05_mixing_vs_preservation.png", p_tradeoff,
       width = 10, height = 7, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/05_mixing_vs_preservation-scaled.png?resize=2048%2C1434&ssl=1)

### 6.4 Selecting the Best Integration Method

```
#-----------------------------------------------
# STEP 14: Select optimal integration method
#-----------------------------------------------

# Rank methods by mixing score (primary criterion)
mixing_results <- mixing_results[order(-mixing_results$mixing_score), ]

cat("\nMethod ranking (by sample mixing):\n")
for (i in1:nrow(mixing_results)) {
  cat(sprintf("  %d. %s (mixing: %.2f, separation: %.2f)\n",
              i,
              mixing_results$method[i],
              mixing_results$mixing_score[i],
              mixing_results$condition_separation[i]))
}

# Select top method
best_method <- mixing_results$method[1]
cat("\nRecommended method:", best_method, "\n")

# For subsequent analysis, we'll use the best method
# In most cases, this will be Harmony or CCA
if (best_method == "Harmony") {
  integrated_final <- integrated_harmony
  reduction_final <- "umap.harmony"
} elseif (best_method == "CCA") {
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
} elseif (best_method == "RPCA") {
  integrated_final <- integrated_rpca
  reduction_final <- "umap.rpca"
} elseif (best_method == "FastMNN") {
  integrated_final <- integrated_fastmnn
  reduction_final <- "umap.mnn"
} else{
  # Default to CCA if method not found
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
  cat("Warning: Method not recognized, defaulting to CCA\n")
}

cat("\nUsing", best_method, "for downstream clustering analysis\n")
# Using CCA for downstream clustering analysis
```

> **Interpretation Guide**:
>
> - **High mixing score (>3.0)**: Good integration, samples mix well within cell types
> - **Moderate separation**: Biological conditions remain distinguishable
> - **Balanced scores**: Ideal integration removes technical variation while preserving biology
>
> **Warning signs**:
>
> - Very high mixing + very low separation: Over-integration (lost biological signal)
> - Very low mixing + very high separation: Under-integration (batches not corrected)

## 7. Understanding Clustering: From Integrated Data to Cell Populations

With successfully integrated data, we can now cluster cells to identify distinct populations. Clustering groups cells with similar gene expression profiles, allowing us to discover cell types, cell states, and rare populations.

### 7.1 The Clustering Algorithm: Graph-Based Community Detection

Seurat uses the **Louvain algorithm** (or its faster variant, Leiden) for clustering:

**Step 1: Build k-nearest neighbour (KNN) graph**
- Each cell connects to its k most similar neighbours
- Similarity measured in integrated PCA space (not raw gene expression)
- Default: k=20 neighbours

**Step 2: Calculate shared nearest neighbours (SNN)**
- Weight edges based on shared neighbours between cells
- Cells with many common neighbours = stronger connection
- Creates more robust clusters

**Step 3: Community detection**
- Apply Louvain/Leiden algorithm to find communities
- Communities = groups of densely connected cells
- Resolution parameter controls granularity

### 7.2 The Resolution Parameter: Controlling Cluster Granularity

The **resolution parameter** controls how many clusters are identified:

**Low resolution (0.2-0.5):**
- Fewer, larger clusters
- Broad cell types (T cells, B cells, monocytes)
- Good for initial exploration
- **Use when**: You want major cell type categories

**Medium resolution (0.6-1.0):**
- Moderate number of clusters
- Cell subtypes (CD4+ T cells, CD8+ T cells, classical monocytes)
- Standard choice for most analyses
- **Use when**: Balanced between detail and interpretability

**High resolution (1.2-2.0):**
- Many small clusters
- Fine-grained cell states (Th1, Th2, Th17 CD4+ T cells)
- May split biological populations unnecessarily
- **Use when**: You need maximum detail or are focusing on specific cell types

==**Rules of thumb:**==
- ==**~3,000 cells**: resolution 0.4-0.6==
- ==**~10,000 cells**: resolution 0.6-0.8==
- ==**~50,000 cells**: resolution 0.8-1.2==
- ==**\>100,000 cells**: resolution 1.0-1.5==

**NB:** *For our dataset (~72,000 cells), we’ll test resolutions from 0.4 to 1.2.*

### 7.3 Common Misconception: Cluster Number Should Match Expected Cell Types

A frequent misunderstanding among beginners is that **the number of clusters should exactly equal the number of cell types you expect to find**. For example, if you expect to find 8 major PBMC cell types (T cells, B cells, monocytes, NK cells, dendritic cells, etc.), you might think you should adjust the resolution until you get exactly 8 clusters.

**This is incorrect.** Here’s why:

**1\. Cell types exist at multiple hierarchical levels**
- At low resolution: “T cells” = 1 cluster
- At medium resolution: “CD4+ T cells” + “CD8+ T cells” = 2 clusters
- At high resolution: “Naive CD4+ T” + “Memory CD4+ T” + “Naive CD8+ T” + “Effector CD8+ T” = 4+ clusters

**2\. Cell states add complexity beyond cell types**
- Activated vs resting cells within the same type
- Proliferating vs quiescent cells
- Stressed or transitional states
- Cell cycle phases

**3\. The “right” number of clusters depends on your question**
- Broad comparison? Use fewer clusters (major types)
- Detailed analysis? Use more clusters (subtypes and states)
- There is no single “correct” answer

**4\. Biological reality is continuous, not discrete**
- Cells exist on spectra and gradients
- Clustering imposes artificial boundaries
- Two similar cell types may legitimately appear as one cluster
- One heterogeneous cell type may split into multiple clusters

==**The proper approach:**==
- ==Test multiple resolutions systematically==
- ==Evaluate biological interpretability at each resolution==
- ==Use marker genes to validate cluster identities==
- ==Choose resolution based on your analysis goals, not preconceived cluster counts==
- ==Remember: 12 clusters representing 8 cell types (with subtypes) is perfectly valid==

**NB:** *For PBMCs, you might expect 8-10 major cell types, but **12-15 clusters is common and appropriate** when accounting for T cell subtypes, monocyte subsets, and cell states. Don’t force your data into an arbitrary number of clusters.*

### 7.4 Clustering at Multiple Resolutions

```
#-----------------------------------------------
# STEP 15: Perform clustering at multiple resolutions
#-----------------------------------------------

# Test multiple resolutions
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)

for (res in resolutions) {
  integrated_final <- FindClusters(
    integrated_final,
    resolution = res,
    verbose = FALSE
  )

  n_clusters <- length(unique(integrated_final@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("Resolution %.1f: %d clusters\n", res, n_clusters))
}

# Rename cluster columns for clarity
colnames(integrated_final@meta.data) <- gsub("RNA_snn_res\\.", "clusters_res_",
                                              colnames(integrated_final@meta.data))
```

### 7.5 Visualizing Multi-Resolution Clustering

```
#-----------------------------------------------
# STEP 16: Visualize clustering results across resolutions
#-----------------------------------------------

# Create UMAP plots for each resolution
plot_list <- lapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)
  n_clusters <- length(unique(integrated_final@meta.data[[cluster_col]]))

  DimPlot(integrated_final,
          reduction = reduction_final,
          group.by = cluster_col,
          label = TRUE,
          label.size = 4,
          pt.size = 0.3) +
    ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
    NoLegend()
})

# Combine plots
combined_resolutions <- wrap_plots(plot_list, ncol = 2)
ggsave("plots/clustering/06_multi_resolution_clustering.png",
       combined_resolutions, width = 14, height = 14, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/06_multi_resolution_clustering-scaled.png?resize=2048%2C2048&ssl=1)

### 7.6 Determining Optimal Resolution

**The Silhouette Score Metric:**

We use **silhouette scores** to objectively assess cluster quality at each resolution:

1. **What it measures**: How well each cell fits into its assigned cluster
- For each cell: Calculate average distance to cells in its own cluster vs other clusters
- Score ranges from -1 to +1
1. **Interpretation**:
- **High score (~0.5-0.7)**: Clusters are well-separated, cells clearly belong to their clusters ✓
- **Moderate score (~0.3-0.5)**: Reasonable clusters, some overlap between groups
- **Low score (<0.3)**: Clusters poorly defined, may be over-clustering or under-clustering
1. **The trade-off**:
- Too few clusters (low resolution) → Biology is oversimplified, different cell types lumped together
- Too many clusters (high resolution) → Over-splitting of natural cell populations

**Our Strategy:**

1. Test multiple resolutions (0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
2. Calculate silhouette score for each (using 5,000-cell subsample)
3. Look for resolution with high silhouette score + interpretable cluster number
4. Visually validate the chosen resolution makes biological sense
```
#-----------------------------------------------
# STEP 17: Evaluate optimal clustering resolution
#-----------------------------------------------

# We'll use silhouette scores to assess cluster quality at each resolution
# Silhouette score measures how similar cells are to their own cluster vs other clusters
# Values range from -1 to 1 (higher is better)
library(cluster)

# Map UMAP reduction name to corresponding integrated reduction name
# We use integrated space (not UMAP) for silhouette calculation because:
# - Integrated space preserves more information (30+ dimensions vs UMAP's 2)
# - Silhouette scores are more meaningful in the space where clustering was performed
integrated_reduction <- switch(reduction_final,
  "umap.cca"= "integrated.cca",
  "umap.rpca"= "integrated.rpca",
  "umap.harmony"= "harmony",
  "umap.mnn"= "integrated.mnn",
  "umap"= "pca"# Fallback for naive merge
)

cat("\nEvaluating clustering resolutions using", integrated_reduction, "space\n")

# For large datasets (>10,000 cells), subsample for silhouette calculation (use 5000 cells in our case)
# Computing distance matrices for all cells is computationally prohibitive
n_cells <- ncol(integrated_final)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  cat("Dataset has", n_cells, "cells - subsampling", max_cells_for_silhouette, "cells for silhouette calculation\n")
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else{
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)

  # Get clusters for subsampled cells
  clusters <- as.numeric(integrated_final@meta.data[[cluster_col]][subsample_idx])

  # Get integrated coordinates for subsampled cells (use first 30 dimensions)
  coords <- Embeddings(integrated_final, reduction = integrated_reduction)[subsample_idx, ]
  coords <- coords[, 1:min(30, ncol(coords))]

  # Calculate silhouette (only if we have at least 2 clusters)
  if (length(unique(clusters)) > 1) {
    dist_matrix <- dist(coords)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  } else{
    NA# Return NA if only 1 cluster
  }
})

# Create resolution comparison table
resolution_comparison <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(res) {
    cluster_col <- paste0("clusters_res_", res)
    length(unique(integrated_final@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

cat("\nClustering resolution comparison:\n")
print(resolution_comparison)

##   resolution n_clusters silhouette_score
## 1        0.4         17        0.2040787
## 2        0.6         23        0.2013915
## 3        0.8         24        0.1951003
## 4        1.0         26        0.1962814
## 5        1.2         30        0.1718389
```

> **Note on Subsampling for Large Datasets**: When datasets have >5,000 cells, computing full distance matrices for silhouette scores becomes computationally prohibitive (would require >25 GB memory for this dataset). Subsampling 5,000 cells is standard practice and provides reliable estimates of cluster quality. The silhouette scores from the subsample accurately reflect the separation quality of the full clustering.

```
# Recommend optimal resolution
# Choose resolution with good silhouette score and interpretable cluster number
optimal_idx <- which.max(resolution_comparison$silhouette_score)
optimal_resolution <- resolution_comparison$resolution[optimal_idx]

cat("\nRecommended resolution:", optimal_resolution,
    "(", resolution_comparison$n_clusters[optimal_idx], "clusters )\n")

# Recommended resolution: 0.4 ( 17 clusters )

# Set the optimal clustering as the default
integrated_final$seurat_clusters <- integrated_final@meta.data[[paste0("clusters_res_", optimal_resolution)]]
Idents(integrated_final) <- "seurat_clusters"
```

> **Silhouette Score Interpretation**:
>
> - **\>0.5**: Strong cluster separation (good)
> - **0.3-0.5**: Moderate separation (acceptable)
> - **<0.3**: Weak separation (may be over-clustering)
>
> **Selection strategy**: Choose the resolution with the highest silhouette score while maintaining biologically interpretable cluster numbers. For PBMCs, 8-15 clusters is typical (major cell types + subtypes).

### 7.7 Cluster Quality Assessment

```
#-----------------------------------------------
# STEP 18: Assess cluster quality and stability
#-----------------------------------------------

# Check cluster sizes
cluster_sizes <- table(integrated_final$seurat_clusters)
cat("\nCluster sizes:\n")
print(cluster_sizes)

# Identify very small clusters (<1% of cells)
min_cluster_size <- 0.01 * ncol(integrated_final)
small_clusters <- which(cluster_sizes < min_cluster_size)

if (length(small_clusters) > 0) {
  cat("\n⚠  Warning: Small clusters detected (<1% of cells):\n")
  print(cluster_sizes[small_clusters])
  cat("   These may represent rare cell types or over-clustering\n")
}

# Check sample distribution across clusters
sample_cluster_table <- table(integrated_final$seurat_clusters, integrated_final$sample_id)
sample_cluster_pct <- prop.table(sample_cluster_table, margin = 1) * 100

# Identify sample-specific clusters (>70% from one sample)
sample_specific <- apply(sample_cluster_pct, 1, max) > 70
if (any(sample_specific)) {
  cat("\n⚠  Warning: Sample-dominated clusters detected (>70% from single sample):\n")
  print(names(which(sample_specific)))
  cat("   These may indicate incomplete batch correction\n")
}

# Visualize sample distribution in clusters
sample_dist_data <- as.data.frame.matrix(sample_cluster_pct)
sample_dist_data$cluster <- rownames(sample_dist_data)
sample_dist_long <- reshape2::melt(sample_dist_data, id.vars = "cluster",
                                    variable.name = "sample", value.name = "percentage")

p_sample_dist <- ggplot(sample_dist_long, aes(x = cluster, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Sample Distribution Across Clusters",
       subtitle = "Check for sample-dominated clusters (poor integration)",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave("plots/clustering/08_sample_distribution_clusters.png", p_sample_dist,
       width = 12, height = 6, dpi = 300)

# Check condition distribution
condition_dist <- table(integrated_final$seurat_clusters, integrated_final$condition)
condition_dist_pct <- prop.table(condition_dist, margin = 1) * 100

cat("\nCondition distribution across clusters:\n")
print(round(condition_dist_pct, 1))

# Visualize condition distribution
condition_dist_data <- as.data.frame.matrix(condition_dist_pct)
condition_dist_data$cluster <- rownames(condition_dist_data)
condition_dist_long <- reshape2::melt(condition_dist_data, id.vars = "cluster",
                                      variable.name = "condition", value.name = "percentage")

p_condition_dist <- ggplot(condition_dist_long, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Condition Distribution Across Clusters",
       subtitle = "Check if biological conditions are preserved",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Healthy"= "#2E86AB",
                                "Post_Treatment"= "#F18F01"))

ggsave("plots/clustering/09_condition_distribution_clusters.png", p_condition_dist,
       width = 12, height = 6, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/08_sample_distribution_clusters-scaled.png?resize=2048%2C1024&ssl=1)

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/09_condition_distribution_clusters-scaled.png?resize=2048%2C1024&ssl=1)

### 7.8 Final Integrated and Clustered Visualization

```
#-----------------------------------------------
# STEP 19: Create comprehensive visualization of final results
#-----------------------------------------------

# Main UMAP with clusters
p_final_clusters <- DimPlot(integrated_final,
                            reduction = reduction_final,
                            group.by = "seurat_clusters",
                            label = TRUE,
                            label.size = 5,
                            pt.size = 0.05) +
  ggtitle(paste0("Final Clustering (", best_method, " Integration)")) +
  theme(plot.title = element_text(face = "bold", size = 14))

# Split by condition to show preservation
p_by_condition <- DimPlot(integrated_final,
                          reduction = reduction_final,
                          group.by = "seurat_clusters",
                          split.by = "condition",
                          label = TRUE,
                          label.size = 4,
                          pt.size = 0.05,
                          ncol = 2) +
  ggtitle("Clusters by Treatment Condition") +
  theme(strip.text = element_text(face = "bold"))

# Colored by sample (integration quality check)
p_by_sample <- DimPlot(integrated_final,
                       reduction = reduction_final,
                       group.by = "sample_id",
                       pt.size = 0.05,
                       cols = sample_colors) +
  ggtitle("Sample Mixing (Integration Quality)") +
  theme(legend.text = element_text(size = 7))

# Colored by condition
p_by_condition_single <- DimPlot(integrated_final,
                                 reduction = reduction_final,
                                 group.by = "condition",
                                 pt.size = 0.05,
                                 cols = condition_colors) +
  ggtitle("Treatment Conditions")

# Combine final visualizations
combined_final <- (p_final_clusters | p_by_condition_single) /
                  (p_by_sample | plot_spacer())

ggsave("plots/clustering/10_final_integrated_clustered.png", combined_final,
       width = 16, height = 12, dpi = 300)

# Split by condition detailed view
ggsave("plots/clustering/11_clusters_by_condition_split.png", p_by_condition,
       width = 18, height = 6, dpi = 300)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/10_final_integrated_clustered-scaled.png?resize=2048%2C1536&ssl=1)

![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/11_clusters_by_condition_split-scaled.png?resize=2048%2C683&ssl=1)

## 8. Saving Integrated and Clustered Data

```
#-----------------------------------------------
# STEP 20: Save integrated and clustered data
#-----------------------------------------------

cat("\n=== Saving Final Results ===\n")

# Join layers for downstream analysis
integrated_final <- JoinLayers(integrated_final)

# Save the final integrated object
saveRDS(integrated_final, "integrated_data/integrated_clustered_seurat.rds")

# Save metadata
write.csv(integrated_final@meta.data,
          "metadata/integrated_cell_metadata.csv",
          row.names = TRUE)

# Create comprehensive summary
integration_summary <- data.frame(
  dataset = "GSE174609",
  n_samples = length(unique(integrated_final$sample_id)),
  n_cells_total = ncol(integrated_final),
  n_genes = nrow(integrated_final),
  integration_method = best_method,
  optimal_resolution = optimal_resolution,
  n_clusters = length(unique(integrated_final$seurat_clusters)),
  mixing_score = mixing_results$mixing_score[mixing_results$method == best_method],
  condition_separation = mixing_results$condition_separation[mixing_results$method == best_method]
)

write.csv(integration_summary, "metadata/integration_summary.csv", row.names = FALSE)
```
![](https://i0.wp.com/ngs101.com/wp-content/uploads/2025/12/scRNAseq_-integration_summary.png?w=1374&ssl=1)

## 9. Understanding Your Integration Results: Interpreting the Plots

### 9.1 What Good Integration Looks Like

After successful integration and clustering, you should observe:

**1\. Sample Mixing (plots/clustering/10\_final\_integrated\_clustered.png, bottom left)**

- ✓ All 8 samples intermingle within each cluster
- ✓ No obvious sample-specific sub-clusters
- ✓ Uniform color distribution across the UMAP

**Red flags:**

- ✗ Samples form separate islands
- ✗ One sample dominates a cluster (>70%)
- ✗ Clear batch structure still visible

**2\. Biological Signal Preservation (plots/clustering/10\_final\_integrated\_clustered.png, top right)**

- ✓ Healthy and Post-treatment conditions are detectable
- ✓ Some clusters show condition-specific enrichment (this is biology!)
- ✓ Treatment effects are visible but don’t drive overall structure

**Red flags:**

- ✗ All conditions perfectly mixed in every cluster (over-integration)
- ✗ No condition-specific patterns at all
- ✗ Loss of expected biological differences

**3\. Cluster Quality (plots/clustering/06\_multi\_resolution\_clustering.png)**

- ✓ Clear cluster boundaries with minimal bridging
- ✓ Clusters remain stable across similar resolutions
- ✓ Clusters have reasonable sizes (not 95% in one cluster, 5% scattered)

**Red flags:**

- ✗ Diffuse clouds with no clear boundaries
- ✗ Many tiny clusters (<50 cells each)
- ✗ Clusters disappear or merge dramatically at small resolution changes

### 9.2 Interpreting Multi-Resolution Clustering

**Low resolution (0.4-0.5):** Major cell types emerge

- Expect: 6-10 clusters for PBMCs
- Typical: “T cells”, “B cells”, “Myeloid cells”, “NK cells”

**Medium resolution (0.6-0.8):** Cell subtypes resolve

- Expect: 10-15 clusters
- Typical: “CD4+ T cells”, “CD8+ T cells”, “Classical monocytes”, “Non-classical monocytes”

**High resolution (1.0-1.2):** Fine-grained states appear

- Expect: 15-20+ clusters
- Typical: “Naive CD4+ T”, “Memory CD4+ T”, “Cytotoxic CD8+ T”, “CD14+ monocytes”, “CD16+ monocytes”

**Choosing resolution:**

- Start with medium (0.6-0.8) for standard analysis
- Use high resolution if focusing on specific cell types
- Use low resolution for broad overview or preliminary analysis

### 9.3 Common Integration Issues and Diagnoses

**Issue 1: Under-Integration**

- **Symptoms**: Samples still separate, mixing score <2.5
- **Diagnosis**: Check plots/integration\_comparison/02\_integration\_by\_sample.png – samples form distinct clusters
- **Solutions**:
- Try more aggressive method (CCA > RPCA, Harmony > FastMNN)
- Increase integration dimensions
- Check if samples are genuinely too different (different tissues, different protocols)

**Issue 2: Over-Integration**

- **Symptoms**: Perfect mixing but lost biological signal, separation score <1.5
- **Diagnosis**: Check plots/integration\_comparison/03\_integration\_by\_condition.png – conditions indistinguishable
- **Solutions**:
- Use less aggressive method (FastMNN instead of Harmony)
- Reduce integration strength (fewer dimensions)
- Check if biological differences are truly subtle

**Issue 3: Sample-Specific Clusters**

- **Symptoms**: One cluster is >70% from single sample
- **Diagnosis**: Check plots/clustering/08\_sample\_distribution\_clusters.png
- **Solutions**:
- Re-run integration with stronger correction
- Check if that sample has genuine unique biology (e.g., rare disease sample)
- Verify QC wasn’t too lenient for that sample

**Issue 4: Over-Clustering**

- **Symptoms**: Many tiny clusters, unstable across resolutions
- **Diagnosis**: Check plots/clustering/06\_multi\_resolution\_clustering.png
- **Solutions**:
- Lower resolution parameter
- Check for over-splitting of continuous populations
- May indicate heterogeneous cell states rather than discrete types

## 10. Best Practices for Integration and Clustering

### 10.1 General Principles

**1\. Always Compare Multiple Integration Methods**

- Don’t assume one method is universally best
- Different methods excel with different data characteristics
- Visual comparison is more informative than any single metric
- Save results from multiple methods for reviewer questions

**2\. Validate Integration Quality Systematically**

- Visual inspection (UMAPs colored by sample and condition)
- Quantitative metrics (mixing scores, silhouette scores)
- Biological validation (known marker genes still differential)
- Sample distribution across clusters

**3\. Choose Resolution Based on Biological Question**

- Broad comparison: Low resolution (major cell types)
- Cell type focus: Medium resolution (subtypes)
- Cell state analysis: High resolution (fine-grained states)
- Can always re-cluster at different resolution later

**4\. Document All Decisions**

- Record which integration method and why
- Note resolution choice and rationale
- Save all comparison plots
- Keep integration quality metrics

### 10.2 Method Selection Guidelines

**Use CCA when:**

- ✓ Standard dataset (10,000-100,000 cells)
- ✓ Different cell type compositions between samples
- ✓ You want robust, well-tested method
- ✓ Default choice for most experiments

**Use RPCA when:**

- ✓ Large datasets (>100,000 cells)
- ✓ Similar cell type composition across samples
- ✓ Computational resources are limited
- ✓ Speed is important

**Use Harmony when:**

- ✓ Very large datasets (>500,000 cells)
- ✓ Need fast integration for exploration
- ✓ Soft correction is preferred
- ✓ Multiple batch variables to correct

**Use FastMNN when:**

- ✓ Subtle biological differences need preservation
- ✓ Conservative correction preferred
- ✓ Complex batch structure
- ✓ Reference-based integration

> **Note**: scVI is another powerful integration method but requires complex Python/PyTorch setup. The four R-based methods above are sufficient for most analyses and much easier to configure.

### 10.3 Resolution Selection Strategy

**Step 1: Survey multiple resolutions**

```
# Test range of resolutions
resolutions <- seq(0.2, 1.5, by = 0.2)
```

**Step 2: Evaluate stability**

- Clusters should be similar at adjacent resolutions
- Major cell types stable across wide range
- Subclusters emerge gradually, not suddenly

**Step 3: Consider cluster sizes**

- Avoid <50 cell clusters unless biologically justified
- Aim for relatively balanced cluster sizes
- Very large clusters (>30%) may need higher resolution

**Step 4: Biological interpretation**

- Can you assign biological identity to most clusters?
- Do cluster numbers match expected cell type diversity?
- For PBMCs: 8-15 clusters typical

**Step 5: Downstream analysis needs**

- Differential expression: Medium resolution (statistical power)
- Trajectory analysis: Lower resolution (avoid fragmenting trajectories)
- Cell-cell communication: Medium-high resolution (detail needed)

### 10.4 Common Pitfalls to Avoid

**Pitfall 1: Using Single Integration Method Without Comparison**

- Problem: May not be optimal for your data
- Solution: Always test 2-3 methods, compare results

**Pitfall 2: Ignoring Biological Signal Loss**

- Problem: Over-aggressive integration removes true differences
- Solution: Check condition separation, validate marker genes

**Pitfall 3: Arbitrary Resolution Choice**

- Problem: Too high (over-clustering) or too low (under-clustering)
- Solution: Test multiple resolutions, use metrics and biological knowledge

**Pitfall 4: Not Validating Integration Quality**

- Problem: Proceed with poorly integrated data
- Solution: Check mixing scores, sample distribution, visual inspection

**Pitfall 5: Treating Clusters as Ground Truth**

- Problem: Clusters are computational constructs, not biological absolutes
- Solution: Validate with marker genes, biological knowledge, orthogonal methods

**Pitfall 6: Ignoring Small Clusters**

- Problem: Miss rare cell types or fail to detect over-clustering
- Solution: Investigate small clusters (<1%) – rare cell type or artifact?

### 10.5 When Integration May Not Be Necessary

**Skip integration if:**

- Single sample experiment
- No batch effects observed (rare)
- Samples genuinely represent different tissues/conditions you want to analyze separately
- Preliminary QC/exploration phase

**Use minimal integration if:**

- Mild batch effects
- Want to preserve maximum biological variation
- Samples very similar in composition

### 10.6 Issue-Specific Solutions

**Issue: Samples Still Separate After Integration**

Diagnosis:

```
# Check mixing score
mixing_results  # Should be >3.0

# Visual check
DimPlot(integrated_obj, group.by = "sample_id")
```

Solutions:

```
# Try more aggressive method
integrated <- IntegrateLayers(
  method = CCAIntegration,  # Instead of RPCAIntegration
  dims = 1:40,              # Increase dimensions
  k.anchor = 20             # Increase anchor strength
)
```

**Issue: Over-Integration (Lost Biological Signal)**

Diagnosis:

```
# Check condition separation
assess_condition_separation(integrated_obj)  # Should be >2.0

# Check marker genes
FeaturePlot(integrated_obj, features = c("CD3D", "CD14"))
```

Solutions:

```
# Use conservative method
integrated <- IntegrateLayers(
  method = FastMNNIntegration,  # Conservative
  dims = 1:20                    # Fewer dimensions
)

# Or reduce correction strength
integrated <- RunHarmony(
  theta = 1,  # Lower theta = less correction (default is 2)
  ...
)
```

**Issue: Very Unbalanced Cluster Sizes**

Diagnosis:

```
table(integrated_obj$seurat_clusters)
# One cluster with 80% of cells?
```

Solutions:

```
# Increase resolution
integrated_obj <- FindClusters(integrated_obj, resolution = 1.2)

# Check if PCA captures enough variation
ElbowPlot(integrated_obj, ndims = 50)

# Try more PCs
integrated_obj <- RunPCA(integrated_obj, npcs = 50)
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:40)
```

**Issue: Too Many Tiny Clusters**

Diagnosis:

```
# Many clusters with <50 cells
cluster_sizes <- table(integrated_obj$seurat_clusters)
sum(cluster_sizes < 50)
```

Solutions:

```
# Decrease resolution
integrated_obj <- FindClusters(integrated_obj, resolution = 0.4)

# Increase k.param (more neighbors)
integrated_obj <- FindNeighbors(integrated_obj, k.param = 30)
```

**Issue: Sample-Specific Clusters**

Diagnosis:

```
# Check sample distribution
table(integrated_obj$seurat_clusters, integrated_obj$sample_id)
```

Solutions:

```
# 1. Re-run integration with stronger correction
# 2. Check if that sample failed QC
# 3. Verify if truly unique biology (e.g., disease-specific cells)

# If technical, remove outlier sample:
integrated_obj <- subset(integrated_obj, subset = sample_id != "problematic_sample")
```

## 11. Conclusion: From Raw Data to Biological Insights

Congratulations! You’ve completed the integration and clustering phase of single-cell RNA-seq analysis.

### 11.1 Critical Insights

**1\. Integration is Not Optional**
For multi-sample experiments, proper integration is essential. Technical variation from sample preparation, sequencing, and processing will overwhelm biological signal without correction.

**2\. Method Choice Matters**
Different integration methods have different strengths. Always compare multiple approaches for your specific dataset rather than blindly applying one method.

**3\. Biology Should Guide Analysis**
While algorithms provide suggestions, biological knowledge should inform decisions about:

- Integration aggressiveness
- Clustering resolution
- Cluster merging or splitting
- Outlier handling

**4\. Validation is Continuous**
Integration and clustering aren’t one-time steps. As you proceed to cell type annotation and differential expression, you may need to revisit and refine these analyses.

### 11.2 What Makes Good Integration

Remember these principles:

- **Samples mix within cell types** (technical variation removed)
- **Conditions remain distinguishable** (biological variation preserved)
- **Known biology is recapitulated** (expected cell types present)
- **Rare populations are retained** (not over-corrected away)
- **Marker genes remain differential** (biological signal intact)

### 11.3 Resources for Continued Learning

**Official Documentation:**

- Seurat v5 Integration vignettes: https://satijalab.org/seurat/articles/seurat5\_integration
- Harmony documentation: https://github.com/immunogenomics/harmony
- FastMNN documentation: https://bioconductor.org/packages/batchelor/

**Recommended Reading:**

- Luecken & Theis (2019) – Current best practices in single-cell RNA-seq
- Tran et al. (2020) – Benchmark of batch-effect correction methods
- Korsunsky et al. (2019) – Fast, sensitive and accurate integration (Harmony)

## 12. References

1. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol. 2024;42(2):293-304. doi:10.1038/s41587-023-01767-y
2. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16(12):1289-1296. doi:10.1038/s41592-019-0619-0
3. Haghverdi L, Lun ATL, Morgan MD, Marioni JC. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol. 2018;36(5):421-427. doi:10.1038/nbt.4091
4. Stuart T, Butler A, Hoffman P, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019;177(7):1888-1902.e21. doi:10.1016/j.cell.2019.05.031
5. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15(12):1053-1058. doi:10.1038/s41592-018-0229-2
6. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol. 2019;15(6):e8746. doi:10.15252/msb.20188746
7. Tran HTN, Ang KS, Chevrier M, et al. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biol. 2020;21(1):12. doi:10.1186/s13059-019-1850-9
8. Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E. Fast unfolding of communities in large networks. J Stat Mech. 2008;2008(10):P10008. doi:10.1088/1742-5468/2008/10/P10008
9. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9(1):5233. doi:10.1038/s41598-019-41695-z
10. Büttner M, Miao Z, Wolf FA, Teichmann SA, Theis FJ. A test metric for assessing single-cell RNA-seq batch correction. Nat Methods. 2019;16(1):43-49. doi:10.1038/s41592-018-0254-1
11. Lee H, Joo JY, Sohn DH, et al. Single-cell RNA sequencing reveals rebalancing of immunological response in patients with periodontitis after non-surgical periodontal therapy. J Transl Med. 2022;20(1):504. doi:10.1186/s12967-022-03686-8
12. Amezquita RA, Lun ATL, Becht E, et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods. 2020;17(2):137-145. doi:10.1038/s41592-019-0654-x
13. Luecken MD et al. Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods*. 2022;19(1):41-50.
14. Gayoso A et al. A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*. 2022;40:163-166.

---

*This tutorial is part of the NGS101.com series on single cell sequencing analysis. If this tutorial helped advance your research, **please comment and share your experience** to help other researchers! **Subscribe to stay updated** with our latest bioinformatics tutorials and resources.*