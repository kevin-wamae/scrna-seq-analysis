---
title: "How to Analyze Single-Cell RNA-seq Data - Complete Beginner's Guide Part 1: From FASTQ to Count Matrix"
source: "https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-1-from-fastq-to-count-matrix/"
author:
  - "[[Lei]]"
published: 2025-12-11
created: 2026-07-12
description: "A comprehensive step-by-step tutorial for analyzing 10x Genomics single-cell RNA sequencing data using Cell Ranger Introduction: Understanding Single-Cell RNA Sequencing The revolution in molecular biology has been marked by our ability to zoom in from cell populations to individual cells. This shift reveals hidden heterogeneity that bulk measurements mask. Single-cell RNA sequencing (scRNA-seq) is one […]"
tags:
  - "clippings"
---
*A comprehensive step-by-step tutorial for analyzing 10x Genomics single-cell RNA sequencing data using Cell Ranger*

## Introduction: Understanding Single-Cell RNA Sequencing

The revolution in molecular biology has been marked by our ability to zoom in from cell populations to individual cells. This shift reveals hidden heterogeneity that bulk measurements mask. Single-cell RNA sequencing (scRNA-seq) is one of the most powerful technologies driving this revolution, allowing us to profile the transcriptome of thousands of individual cells simultaneously. This capability has transformed our understanding of cellular diversity, developmental trajectories, disease mechanisms, and therapeutic responses.

### What is Single-Cell RNA Sequencing?

Single-cell RNA sequencing measures gene expression at the individual cell level, rather than averaging across millions of cells as in traditional bulk RNA-seq. This distinction is crucial because:

**Cell Heterogeneity:** Even within seemingly uniform tissue samples or cell populations, individual cells can exhibit substantial differences in gene expression. These differences may reflect:

- Different cell types within mixed populations
- Cells at different stages of development or differentiation
- Variable responses to stimuli or environmental conditions
- Rare cell populations with specialized functions
- Cells transitioning between states

**Biological Insights:** By preserving single-cell resolution, scRNA-seq enables researchers to:

- Identify and characterize rare cell populations
- Reconstruct developmental trajectories and lineage relationships
- Discover novel cell types and states
- Understand cellular responses to perturbations at unprecedented detail
- Map disease-associated cell populations in complex tissues

> **Beginner’s Analogy:** If bulk RNA-seq is like measuring the average temperature of a city, single-cell RNA-seq is like measuring the temperature of every individual building. You might discover that while the average is 70°F, some buildings are sweltering at 85°F while others are freezing at 55°F—information completely lost in the average.

### The 10x Genomics Chromium Platform

Among various scRNA-seq technologies, the 10x Genomics Chromium platform has emerged as the gold standard for high-throughput single-cell analysis. Its success stems from several key innovations:

**Droplet-Based Cell Capture:**  
The Chromium system uses microfluidics to partition thousands of individual cells into nanoliter-scale oil droplets (called GEMs – Gel Bead-in-Emulsions). Each droplet contains:

- A single cell
- A gel bead coated with barcoded oligonucleotides
- Reverse transcription reagents

**Unique Cellular Barcoding:**  
Each gel bead carries millions of copies of oligonucleotides with:

- A cell barcode (unique to that gel bead, ~750,000 different barcodes available)
- A UMI (Unique Molecular Identifier) for counting individual mRNA molecules
- A poly(dT) sequence for capturing polyadenylated mRNA
- Primer sequences for PCR amplification

When the cell lyses inside the droplet, its mRNA molecules bind to the barcoded oligos. This means all mRNA from a single cell receives the same cell barcode, while each individual mRNA molecule gets a unique UMI. This elegant design allows us to:

- Pool thousands of cells together after capture
- Sequence everything in a single run
- Computationally assign reads back to their cell of origin using cell barcodes
- Count unique mRNA molecules (not PCR duplicates) using UMIs

**High Throughput and Sensitivity:**  
The 10x platform enables:

- Capture of 500-20,000 cells per sample in minutes
- Detection of 1,000-5,000 genes per cell on average
- Processing of up to 8 samples simultaneously

**Chemistry Versions:**  
10x Genomics has released several chemistry versions over the years:

- **v2 Chemistry (2017-2019):** The original widely-adopted version
- **v3 Chemistry (2019-2023):** Improved sensitivity, ~30% more genes detected per cell
- **v3.1 Chemistry (2021-2023):** Enhanced dual indexing, 98bp reads (vs 91bp in v3)
- **GEM-X Chemistry (2023-present):** Current version, higher throughput (up to 80% cell recovery)

> **Important:** The differences (read length, indexing) are automatically detected from your FASTQ files.

For this tutorial, we will work with **v3.1 chemistry data** which is one of the most widely used formats in current research.

### Understanding the Experimental Workflow

Before diving into computational analysis, it is essential to understand the experimental steps that generate your data:

**Step 1: Sample Preparation**

- Cells or nuclei are isolated from tissue or blood
- Viability should be >70%, ideally >85%
- Cell concentration adjusted to 700-1,200 cells/μL
- Dead cells can be removed using commercial kits if viability is low

**Step 2: GEM Generation and Barcoding (10x Chromium Instrument)**

- Cells are loaded into a microfluidic chip
- Cells + gel beads + oil combine to form droplets (GEMs)
- Takes ~6-7 minutes to process 8 samples
- Cells lyse and mRNA binds to barcoded beads

**Step 3: cDNA Synthesis and Amplification**

- Reverse transcription occurs inside droplets
- GEMs are broken, cDNA is pooled and amplified
- Results in full-length cDNA with cell barcodes and UMIs

**Step 4: Library Preparation**

- cDNA is enzymatically fragmented
- Illumina sequencing adapters are added
- Libraries are size-selected and quantified
- Typical library size: 400-600 bp

**Step 5: Sequencing**

- Libraries are sequenced on Illumina platforms (NovaSeq, NextSeq, HiSeq)
- Read structure for 3′ gene expression:
- **Read 1 (28 bp):** 16 bp cell barcode + 12 bp UMI
- **Read 2 (98 bp for v3.1, 91 bp for v3):** cDNA insert (mRNA sequence)
- **Index reads (8-10 bp):** Sample indices for demultiplexing
- Recommended depth: 20,000-50,000 reads per cell

**Step 6: Data Processing (What We’ll Do in This Tutorial)**

- Demultiplexing: Separate samples by index sequences (done by facility)
- Quality control: Assess read quality
- Alignment: Map reads to reference genome
- UMI counting: Generate gene expression matrices
- Cell calling: Distinguish real cells from empty droplets

### The Biological Context: Our Example Dataset

For this tutorial, we will analyze data from **GSE174609**, a study investigating immune cell responses in periodontitis patients before and after treatment.

**Scientific Background:**  
Periodontitis is a chronic inflammatory disease affecting the tissues surrounding teeth. While primarily a local oral disease, emerging evidence suggests systemic impacts through:

- Elevated circulating inflammatory mediators
- Bacteremia from periodontal pathogens
- Changes in systemic immune cell populations
- Associations with cardiovascular disease, diabetes, and rheumatoid arthritis

**Research Question:**  
How does non-surgical periodontal therapy (scaling and root planing) affect the systemic immune cell landscape, and can these changes explain the link between oral and systemic health?

**Experimental Design:**

- **Samples:** 12 total
- 4 healthy donors (controls)
- 4 periodontitis patients before treatment (pre-treatment)
- 4 periodontitis patients after treatment (post-treatment, same individuals)
- **Cell Type:** Peripheral blood mononuclear cells (PBMCs)
- **Technology:** 10x Genomics Chromium 3′ v3.1
- **Goal:** Identify immune cell populations affected by periodontitis and treatment

**Expected Cell Types in PBMCs:**

- T cells (CD4+, CD8+, regulatory, memory, naive)
- B cells (naive, memory, plasma cells)
- Monocytes (classical, intermediate, non-classical)
- NK cells (natural killer cells)
- Dendritic cells
- Potentially rare populations

## Setting Up Your Analysis Environment

Before beginning analysis, we need to install Cell Ranger and prepare the necessary software environment.

### Understanding System Requirements

Single-cell RNA-seq analysis is computationally intensive. Before proceeding, ensure your system meets these requirements:

**Optimal Setup:**

- **HPC Access:** University cluster or cloud computing (AWS, Google Cloud)
- **RAM:** 128 GB or more
- **Storage:** 1-2 TB for multiple datasets
- **CPU:** 16+ cores for faster processing

> **Note for Beginners:** If you do not have access to an HPC system, contact your institution’s IT department about computational resources. Most universities provide HPC access for researchers, though cloud computing is another option that can be expensive for large datasets.

### Installing Cell Ranger

Cell Ranger is 10x Genomics’ official software suite for processing single-cell data, handling everything from demultiplexing to gene expression matrix generation.

**Installation Steps:**

```bash
#-----------------------------------------------
# STEP 0: Install Cell Ranger
#-----------------------------------------------
 
# Create a directory for Cell Ranger installation
mkdir-p ~/software
cd~/software
 
# Download Cell Ranger (version 10.0.0)
# Note: Check https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links for the most current version
wget -O cellranger-10.0.0.tar.gz "YOUR_CELLRANGER_DOWNLOAD_LINK_HERE"
 
# Extract the tarball
tar-xzvf cellranger-10.0.0.tar.gz
 
# Add Cell Ranger to PATH (optional: add to ~/.bashrc for permanent access)
exportPATH=$HOME/software/cellranger-10.0.0:$PATH
 
# To make this permanent, add to your ~/.bashrc:
# echo 'export PATH=$HOME/software/cellranger-10.0.0:$PATH' >> ~/.bashrc
# source ~/.bashrc
 
# Check Cell Ranger components
cellranger -h
```

**Understanding Cell Ranger Components:**

Cell Ranger includes several key pipelines:

1. **cellranger mkfastq:** Demultiplexes raw BCL files into FASTQ files (rarely needed and is deprecated – facilities typically provide FASTQs)
2. **cellranger count:** Primary pipeline – aligns reads, counts genes, detects cells
3. **cellranger aggr:** Combines multiple samples for integrated analysis
4. **cellranger reanalyze:** Re-runs secondary analysis with different parameters
5. **cellranger mat2csv:** Converts output matrices to CSV format

For this tutorial, we will primarily use `cellranger count`.

### Installing Additional Software

We’ll also need several supporting tools for quality control and data handling:

```bash
#-----------------------------------------------
# STEP 1: Install supporting software using Conda
#-----------------------------------------------
 
# Create a dedicated conda environment for single-cell analysis
conda create -p ~/Env_scRNApython=3.10
 
# Activate the environment
conda activate ~/Env_scRNA
 
# Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --setchannel_priority strict
 
# Install required tools
conda install-y \
    sra-tools \          # For downloading data from SRA
    fastqc \             # For quality control
    multiqc \            # For aggregating QC reports
    wget \               # For downloading files
    samtools             # For BAM file manipulation
 
# Verify installations
fastqc --version         # FastQC v0.12.1
multiqc --version        # multiqc, version 1.19
prefetch --version       # prefetch : 3.0.10
```

> **Why These Tools?**
> 
> - **SRA-tools:** Downloads data from NCBI’s Sequence Read Archive
> - **FastQC:** Generates quality control reports for FASTQ files
> - **MultiQC:** Aggregates multiple FastQC reports into a single summary
> - **samtools:** Views and manipulates BAM alignment files

### Downloading Reference Genome

Cell Ranger requires a pre-built reference genome for read alignment. 10x Genomics provides optimized reference packages:

```bash
#-----------------------------------------------
# STEP 2: Download Cell Ranger reference genome
#-----------------------------------------------
 
# Create a directory for reference genomes
mkdir-p ~/references/cellranger
cd~/references/cellranger
 
# Download human reference genome (GRCh38/hg38)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
 
# Download mouse reference genome (GRCm39)
# wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
 
# Download rat reference genome (mRatBN7.2)
# wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mRatBN7-2-2024-A.tar.gz"
 
# Extract the reference
tar-xzvf refdata-gex-GRCh38-2024-A.tar.gz
 
# The extracted directory contains:
# - genome.fa: Reference genome FASTA
# - genes.gtf: Gene annotations
# - star/: Pre-built STAR aligner index
```

> **Important:** Always use the same reference genome version for samples you will compare later. Mixing reference versions will cause problems in downstream analysis.

With Cell Ranger installed and the reference genome downloaded, we can now examine the structure of single-cell sequencing data.

## Understanding Your Data: FASTQ Files and 10x File Structure

Before processing data, let us understand what we are working with.

### The 10x FASTQ File Naming Convention

When you receive data from a sequencing facility (or download from public repositories), you will encounter specific file naming patterns:

**Standard 10x File Naming:**

```bash
SampleName_S1_L001_R1_001.fastq.gz
SampleName_S1_L001_R2_001.fastq.gz
SampleName_S1_L001_I1_001.fastq.gz  # (optional - i7 sample index)
SampleName_S1_L001_I2_001.fastq.gz  # (optional - i5 sample index for dual indexing)
```

**Breaking Down the Filename:**

- **SampleName:** Your sample identifier
- **S1:** Sample number (S1, S2, S3, etc.)
- **L001:** Lane number (L001, L002, etc. for multiple lanes)
- **R1/R2/I1/I2:** Read type
- R1: Read 1 (cell barcode + UMI)
- R2: Read 2 (cDNA sequence)
- I1: Index read 1 (i7 sample index for demultiplexing)
- I2: Index read 2 (i5 sample index for dual-indexed libraries)
- **001:** File segment (for very large files split into chunks)

**What Each Read Contains:**

**Read 1 (R1) – 28 bases:**

```bash
Structure: [16 bp cell barcode][12 bp UMI]
Example:  AAACCTGAGCGATATA GCACCAACGA
          └──cell barcode─┘ └───UMI───┘
```
- Cell barcode identifies which cell the read came from
- UMI identifies the specific mRNA molecule

**Read 2 (R2) – 98 bases (v3.1 chemistry):**

```bash
Structure: [cDNA sequence from mRNA]
Example:  ATGAAGATCCTGACCTCAGATGATCCATGTCCTGACAACACCTTGAAGGACAAGTTTGACCCTAAGGTGAAGAAGATCCAGAACATTGTACGA...
```
- This is the actual gene sequence we will align and count
- **Note:** v3 chemistry uses 91bp R2 reads; v3.1 uses 98bp for better gene coverage

**Index Read 1 (I1) – 8-10 bases:**

```bash
Structure: [i7 sample index for demultiplexing]
Example:  ATCACG
```
- Used to separate different samples sequenced together (single indexing)
- Typically handled by the sequencing facility
- Cell Ranger does not need this for already-demultiplexed data

**Index Read 2 (I2) – 8-10 bases (optional):**

```bash
Structure: [i5 sample index for dual-indexed libraries]
Example:  GTAACG
```
- Second sample index used in dual-indexed sequencing runs
- Provides additional specificity when multiplexing many samples (>96)
- Reduces index hopping on patterned flow cells (NovaSeq, NextSeq 2000)
- Like I1, this is also not needed by Cell Ranger for demultiplexed data

> **Why Dual Indexing (I1 + I2)?**  
> Dual indexing uses both i7 and i5 sample indices, providing:
> 
> - **Increased multiplexing capacity:** Can pool more samples per lane
> - **Reduced index hopping:** On patterned flow cells, free-floating indices can attach to wrong templates. With dual indexing, a read needs BOTH indices to match for assignment to a sample, drastically reducing misassignment
> - **Better sample identification:** Unique i7+i5 combinations provide more specific sample barcodes
> 
> For already-demultiplexed FASTQ files (from SRA or your sequencing facility), you can safely ignore or remove I1 and I2 files (if present) – Cell Ranger only uses R1 (cell barcode + UMI) and R2 (cDNA sequence).

### FASTQ File Format

Each read in a FASTQ file is represented by 4 lines:

```bash
@SRR14575500.1 1 length=28
AAACCTGAGCGATATAAGAGTGGCAA
+
FFFFFFFFFFFFFFFFFFFFFFFF55
```

**Line-by-Line Breakdown:**

1. **Header:** Starts with `@`, contains read ID and metadata
2. **Sequence:** The actual DNA bases (A, C, G, T, N)
3. **Separator:** Just a `+` symbol or same as the first line
4. **Quality Scores:** Phred quality scores (encoded as ASCII characters)

**Quality Score Interpretation:**

- Each character represents a quality score (Q-score)
- Higher ASCII value = better quality
- Common encoding: Phred+33 (FastQ-Sanger)
- `F` = Q score 37 (error rate 0.0002, 99.98% accuracy)
- `5` = Q score 20 (error rate 0.01, 99% accuracy)
- `!` = Q score 0 (lowest quality)

**Rule of Thumb:**

- Q30 (99.9% accuracy) or higher is considered good quality
- Q20 (99% accuracy) is acceptable
- Below Q20 may need quality trimming

## Data Demultiplexing: Understanding the First Step

Before we can analyze single-cell data, we need to understand how sequencing facilities process raw instrument output into usable FASTQ files.

### What is Demultiplexing?

Modern Illumina sequencers generate **BCL files** (base call files) – a raw binary format containing sequences from ALL samples sequenced together in a single run. Since multiple samples (often 8-96) are pooled on the same flow cell lane to maximize efficiency, we need to **demultiplex** – separate the pooled data into individual sample files based on sample indices.

> **Important Note:** This demultiplexing step is **almost always performed by your sequencing core facility**. They will provide you with demultiplexed FASTQ files ready for Cell Ranger analysis. This section explains what happens behind the scenes to help you understand the data you receive.

### Understanding Sample Indices vs. Cell Barcodes

It’s crucial to distinguish between two types of barcodes in 10x data:

**Sample Indices (for Demultiplexing):**

- Short sequences (8-10 bp) added during library preparation
- Identify which *sample* a read belongs to (e.g., Patient 1 vs Patient 2)
- Used to separate multiplexed samples after sequencing
- Present in the I1 (and I2 for dual indexing) read
- Removed after demultiplexing – NOT in your final FASTQ files

**Cell Barcodes (for Single-Cell Identity):**

- Longer sequences (16 bp) from 10x gel beads
- Identify which *cell* within a sample a read belongs to
- Present in Read 1 (R1)
- Remain in FASTQ files and are used by Cell Ranger
- One of ~750,000 possible barcodes per 10x kit

**Analogy:**

- Sample indices are like apartment building addresses (123 Main St vs 456 Oak Ave)
- Cell barcodes are like apartment numbers within a building (Apt 1, Apt 2, etc.)

### What Happens During Demultiplexing

**Step 1: Sample Sheet Preparation**

Your sequencing facility needs information about which indices correspond to which samples. This is provided in a **sample sheet** (CSV format):

```bash
Lane,Sample,Index
1,Healthy_Donor_1,SI-TT-A1
1,Healthy_Donor_2,SI-TT-A2
1,Patient_Pre_1,SI-TT-A3
1,Patient_Pre_2,SI-TT-A4
```

For 10x samples, the indices look like `SI-TT-A1`, which actually encodes 4 different Illumina index sequences (10x uses a specific indexing scheme for improved performance).

**Step 2: BCL to FASTQ Conversion**

The facility runs Illumina’s `bcl-convert` (or older `bcl2fastq`) software, which:

1. Reads raw BCL files from the sequencer
2. Matches index sequences to samples in the sample sheet
3. Creates separate FASTQ files for each sample
4. Generates quality metrics (Q30 scores, read counts, etc.)

**Step 3: Quality Control**

The facility checks:

- Total reads per sample (typically 200-500 million reads for a 10x run)
- Index read quality (should be >95% Q30)
- Per-sample read distribution (should be relatively even)
- Undetermined reads (reads that do not match any index, should be <10%)

### What You Receive from the Sequencing Facility

After demultiplexing, you will receive:

**For Each Sample:**

- `SampleName_S1_L001_R1_001.fastq.gz` – Cell barcodes + UMIs (required)
- `SampleName_S1_L001_R2_001.fastq.gz` – cDNA sequences (required)
- `SampleName_S1_L001_I1_001.fastq.gz` – i7 sample index (optional, Cell Ranger ignores)
- `SampleName_S1_L001_I2_001.fastq.gz` – i5 sample index (optional, only in dual-indexed runs)

**Quality Control Reports:**

- `Stats.json` – Detailed sequencing statistics
- `Top_Unknown_Barcodes.html` – Most common unmatched indices (diagnostics)
- `Demultiplex_Stats.csv` – Per-sample read counts and quality metrics

**Expected Metrics:**

- **Read depth:** 20,000-50,000 reads per cell × estimated cell number
- For 5,000 cells: 100-250 million reads per sample
- Our PBMC samples: ~150-200 million reads each
- **Q30 Score:** >75% of bases should be Q30 or higher
- **Read length:**
- R1: 28 bp (fixed – cell barcode + UMI)
- R2: 98 bp (v3.1 chemistry)
- I1/I2: 10 bp each (if present)

> **When to Contact the Facility:**
> 
> - Drastically uneven read distribution (some samples have 10x more than others)
> - Very high undetermined reads (>15%)
> - Low Q30 scores (<70%)
> - Missing files

### For Advanced Users: Running Demultiplexing Yourself

In rare cases, you might need to perform demultiplexing yourself (e.g., you have raw BCL files from a collaborator). The standard workflow would be:

**Option 1: Using Cell Ranger mkfastq (10x-specific):**

> **DEPRECATION WARNING:** The `cellranger mkfastq` pipeline is deprecated and will be removed in a future Cell Ranger release. 10x Genomics now recommends using Illumina’s `bcl-convert` directly for demultiplexing. If you are starting a new workflow, use Option 2 below.

```bash
# Requires: BCL files, sample sheet, bcl2fastq, Cell Ranger installation
cellranger mkfastq --run=/path/to/bcl\
                   --csv=sample_sheet.csv \
                   --output-dir=fastqs
```

**Option 2: Using Illumina bcl-convert (Recommended):**

```bash
# Note: bcl-convert requires Illumina customer registration
# Not freely available for public download
bcl-convert --bcl-input-directory /path/to/bcl\
            --output-directory /path/to/fastqs\
            --sample-sheet sample_sheet.csv
```

**We will not cover this in detail** because:

1. Facilities almost always handle this step
2. BCL files are very large (often 1-2 TB per run)
3. The process is facility-specific
4. Most researchers never see BCL files
5. `bcl-convert` requires Illumina registration and is not freely available

If you absolutely need to demultiplex yourself, consult your facility’s documentation or 10x Genomics support.

## Downloading the Tutorial Dataset

Now let us download our example dataset: GSE174609 from the Sequence Read Archive (SRA).

### Understanding the Dataset

**Study:** Peripheral blood mononuclear cells (PBMCs) from periodontitis patients

**Samples (12 total):**

- Healthy donors: 4 samples
- Pre-treatment patients: 4 samples
- Post-treatment patients: 4 samples (same individuals as pre-treatment)

**SRA Accessions:** SRR14575500 – SRR14575511

**Important Note on SRA 10x Data:**  
Public 10x Genomics data in SRA can be stored in three different formats:

1. **Properly formatted:** R1 and R2 FASTQs (easiest to use)
2. **BAM format:** Requires conversion with 10x’s bamtofastq tool
3. **Three-file format:** R1, R2, and I1 (need to identify which is which by file size)

We will show you how to handle the most common case and provide troubleshooting for the others.

### Creating Project Directory Structure

First, let us set up an organized directory structure:

```bash
#-----------------------------------------------
# STEP 3: Create project directory structure
#-----------------------------------------------
 
# Create main project directory
mkdir-p ~/GSE174609_scRNA
cd~/GSE174609_scRNA
 
# Create subdirectories for organization
mkdir-p raw_data          # For downloaded FASTQ files
mkdir-p fastqc_reports    # For quality control reports
mkdir-p cellranger_output # For Cell Ranger results
mkdir-p scripts           # For analysis scripts
mkdir-p logs              # For log files
 
# Verify structure
tree -L 1 .
# Expected output:
# .
# ├── raw_data
# ├── fastqc_reports
# ├── cellranger_output
# ├── scripts
# └── logs
```

### Downloading Data from SRA

For this tutorial, we will demonstrate the workflow using **one sample: SRR14575500** (Healthy\_1). The exact same process applies to all other 11 samples in the dataset.

> **Note:** We are processing a single sample here for educational purposes. In a real analysis, you would download and process all 12 samples using a loop or parallel processing. The commands shown here can be easily adapted to process multiple samples – simply repeat the workflow for each SRR accession (SRR14575500 through SRR14575511).

```bash
#-----------------------------------------------
# STEP 4: Download one sample from GSE174609
#-----------------------------------------------
 
cd~/GSE174609_scRNA/raw_data
 
# We'll download SRR14575500 (Healthy donor #1)
# This is a representative PBMC sample from the study
SRR="SRR14575500"
 
# Step 1: Download SRA file
# This creates a directory ~/ncbi/public/sra/SRR14575500.sra
prefetch $SRR
 
# Step 2: Convert to FASTQ with split files
# --split-files: Create separate R1 and R2 files
# --include-technical: Include index reads if present (I1, I2)
# --threads: Use multiple cores for faster conversion
# --progress: Show progress bar
fasterq-dump --split-files \
             --include-technical \
             --threads 8 \
             --progress \
             --outdir . \
             $SRR/${SRR}.sra
 
# Step 3: Compress FASTQ files to save space
# gzip reduces file size by ~70%
gzip${SRR}*.fastq
 
# Step 4: Clean up SRA file to save disk space
rm-rf $SRR
```

After successful download, you should have for each sample:

```
SRR14575500_1.fastq.gz   # R1: Cell barcodes + UMIs
SRR14575500_2.fastq.gz   # R2: cDNA sequences
```

### Renaming Files for Cell Ranger Compatibility

Cell Ranger expects specific file naming conventions. Let’s rename our downloaded file:

```bash
#-----------------------------------------------
# STEP 5: Rename files with biological sample name
#-----------------------------------------------
 
cd~/GSE174609_scRNA/raw_data
 
# SRR14575500 corresponds to Healthy donor #1
# Cell Ranger expects format: SampleName_S1_L001_R1_001.fastq.gz
 
# Rename R1 (cell barcodes + UMIs)
mvSRR14575500_1.fastq.gz Healthy_1_S1_L001_R1_001.fastq.gz
 
# Rename R2 (cDNA sequences)
mvSRR14575500_2.fastq.gz Healthy_1_S1_L001_R2_001.fastq.gz
```

**Understanding the Naming Convention:**

- **SampleName:** Biological sample name (e.g., Healthy\_1)
- **S1:** Sample index (S1 for first sample, S2 for second, etc.)
- **L001:** Lane 1 (if data was from multiple lanes, you’d have L001, L002, etc.)
- **R1/R2:** Read 1 or Read 2
- **001:** Segment number (for split files)

> **Important:** If your original data came from multiple lanes (L001, L002), Cell Ranger will automatically merge them during processing. Just make sure all lanes for the same sample have the same sample name prefix.

## Quality Control with FastQC

Before running Cell Ranger, let us assess the quality of our FASTQ files.

### Running FastQC

FastQC analyzes sequence quality, GC content, adapter contamination, and other metrics:

```bash
#-----------------------------------------------
# STEP 7: Quality control with FastQC
#-----------------------------------------------
 
cd~/GSE174609_scRNA
 
# Create output directory
mkdir-p fastqc_reports
 
# Run FastQC on all FASTQ files
# We'll run on R2 only (cDNA reads)
fastqc raw_data/*_R2_001.fastq.gz \
    --outdir fastqc_reports \
    --threads 8 \
    --quiet
 
# Aggregate results with MultiQC
multiqc fastqc_reports \
    --outdir fastqc_reports \
    --filename multiqc_report \
    --title "GSE174609 Quality Control"
```

**Why We Focus on R2:**

- **R1 (Cell Barcode + UMI):** Expected to have low sequence complexity
- Will “fail” several FastQC modules (this is normal!)
- Per-base sequence content will be uneven (expected)
- Sequence duplication will be high (expected – only 750K possible barcodes)
- **R2 (cDNA):** Should look like good RNA-seq data
- High per-base quality (>Q30)
- Diverse sequence content
- Should pass most FastQC modules

### Interpreting FastQC Results

**Key Metrics to Evaluate:**

**1\. Per Base Sequence Quality**

- **Good:** Green across all bases, Q scores >28
- **Acceptable:** Yellow in the first few bases or at the end
- **Poor:** Red zones, Q scores <20
- *For our data:* Should see high quality (>Q30) across most of the read

**2\. Per Sequence Quality Scores**

- **Good:** Sharp peak at Q35-40
- **Poor:** Broad distribution or peak at low quality
- *For our data:* Expect peak around Q36-38

**3\. Adapter Content**

- **Good:** No adapters detected
- **Poor:** Adapters present (>5%)
- *For our data:* Should be minimal (<1%) as facilities remove adapters

## Running Cell Ranger Count

Now for the main event: processing our single-cell data with Cell Ranger.

### Understanding Cell Ranger Count

`cellranger count` is the core Cell Ranger pipeline. It performs:

1. **Read alignment:** Maps R2 (cDNA) reads to reference genome using STAR aligner
2. **Barcode processing:** Corrects cell barcode errors in R1 reads
3. **UMI counting:** Collapses PCR duplicates using UMIs, counts unique molecules
4. **Cell calling:** Distinguishes real cells from empty droplets
5. **Gene quantification:** Creates gene × cell expression matrix
6. **Secondary analysis:** Performs clustering, t-SNE, differential expression

**Key Parameters:**

- `--id`: Output folder name (required)
- `--fastqs`: Directory containing FASTQ files
- `--sample`: Sample name prefix in FASTQ filenames
- `--transcriptome`: Path to reference genome
- `--expect-cells`: Expected cell recovery (optional but helpful)
- `--chemistry`: Chemistry version (SC3Pv3 for both v3 and v3.1)
- `--create-bam`: Set to `false` to skip BAM file generation (saves 10-50 GB per sample)
- `--output-dir`: Optional – specify output directory (if not set, uses current directory)
- `--localcores`: Number of CPU cores to use
- `--localmem`: Memory limit in GB

> **Chemistry Parameter Note:** You can also use `--chemistry=auto` for automatic detection, but explicit specification is recommended for reproducibility.
> 
> **BAM File Note:** By default, Cell Ranger creates large BAM files (10-50 GB) containing aligned reads. These are only needed for genome browser visualization or custom alignment analysis. For standard single-cell analysis, set `--create-bam=false` to save significant disk space.

### Running Cell Ranger Count

Let’s process our sample with Cell Ranger:

```bash
#-----------------------------------------------
# STEP 8: Run Cell Ranger count on single sample
#-----------------------------------------------
 
cd~/GSE174609_scRNA/cellranger_output
 
# Configuration
FASTQ_DIR=~/GSE174609_scRNA/raw_data
TRANSCRIPTOME=~/references/cellranger/refdata-gex-GRCh38-2024-A
SAMPLE="Healthy_1"
 
# System resources (adjust based on your system)
CORES=16
MEM=64
 
# Expected cells for PBMC samples (typically 3,000-8,000 cells)
EXPECTED_CELLS=5000
 
# Run Cell Ranger count
cellranger count \
    --id=${SAMPLE} \
    --fastqs=${FASTQ_DIR} \
    --sample=${SAMPLE} \
    --transcriptome=${TRANSCRIPTOME} \
    --expect-cells=${EXPECTED_CELLS} \
    --localcores=${CORES} \
    --localmem=${MEM} \
    --chemistry=SC3Pv3 \
    --create-bam=false
```
``
> **Processing Multiple Samples:**  
> To process all 12 samples in the dataset, use a loop or create an array job on HPC:

```bash
#!/bin/bash
Create sample list
SAMPLES=("Healthy_1""Healthy_2""Healthy_3""Healthy_4"\
         "Pre_Patient_1""Pre_Patient_2""Pre_Patient_3""Pre_Patient_4"\
         "Post_Patient_1""Post_Patient_2""Post_Patient_3""Post_Patient_4")

# Loop through all samples
forSAMPLE in"${SAMPLES[@]}"; do
    cellranger count \
        --id=${SAMPLE} \
        --output-dir~/GSE174609_scRNA/cellranger_output\
        --fastqs=${FASTQ_DIR} \
        --sample=${SAMPLE} \
        --transcriptome=${TRANSCRIPTOME} \
        --expect-cells=5000 \
        --localcores=16 \
        --localmem=64 \
        --chemistry=SC3Pv3 \
        --create-bam=false
done
```


> **For HPC Users (SLURM):**

```bash
#!/bin/bash
#SBATCH --partition=partition_name              # Adjust based on your HPC partitions
#SBATCH --cpus-per-task=16             # Number of CPU cores
#SBATCH --mem=120G                     # Memory allocation (120GB recommended for human genome)
#SBATCH --job-name=cellranger_array
#SBATCH --array=1-12                   # Process samples 1-12 in parallel
#SBATCH --output=logs/cellranger_%A_%a.out    # %A = job ID, %a = array task ID
#SBATCH --error=logs/cellranger_%A_%a.err
#SBATCH --mail-user=your.email@institution.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --time=0-4:00:00               # Time limit per sample

#-----------------------------------------------
# Cell Ranger Count - Array Job for Multiple Samples
#-----------------------------------------------

# Activate conda environment with Cell Ranger
source~/miniforge3/etc/profile.d/conda.sh
conda activate ~/Env_scRNA

# Define all sample names
SAMPLES=("Healthy_1""Healthy_2""Healthy_3""Healthy_4"\
         "Pre_Patient_1""Pre_Patient_2""Pre_Patient_3""Pre_Patient_4"\
         "Post_Patient_1""Post_Patient_2""Post_Patient_3""Post_Patient_4")

# Map SLURM array task ID to sample name
# SLURM arrays are 1-indexed, bash arrays are 0-indexed
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define paths
PROJECT_DIR=~/GSE174609# Main project directory
FASTQ_DIR=${PROJECT_DIR}/fastq# Directory containing FASTQ files
OUTPUT_DIR=${PROJECT_DIR}/cellranger_output
REFERENCE=~/references/cellranger/refdata-gex-GRCh38-2024-A

# Cell Ranger parameters
EXPECTED_CELLS=5000                    # Expected number of cells (adjust based on experiment)
LOCALCORES=16                          # Must match --cpus-per-task
LOCALMEM=110                           # Leave ~10GB buffer for system (120 - 10 = 110)

# Create necessary directories
mkdir-p ${OUTPUT_DIR}
mkdir-p ${PROJECT_DIR}/logs

# Change to output directory (Cell Ranger creates subdirectories here)
cd${OUTPUT_DIR}

# Run Cell Ranger Count
cellranger count \
    --id=${SAMPLE_ID} \
    --fastqs=${FASTQ_DIR} \
    --sample=${SAMPLE_ID} \
    --transcriptome=${REFERENCE} \
    --expect-cells=${EXPECTED_CELLS} \
    --localcores=${LOCALCORES} \
    --localmem=${LOCALMEM} \
    --chemistry=auto \
    --create-bam=false```


### Understanding Cell Ranger Output

Cell Ranger creates an extensive directory structure:

```bash
cellranger_output/
└── Healthy_1/
    ├── outs/
    │   ├── web_summary.html                        # Main QC report (VIEW THIS FIRST)
    │   ├── metrics_summary.csv                     # Key metrics in CSV format
    │   │
    │   ├── filtered_feature_bc_matrix/             # **CELLS ONLY** - use for analysis
    │   │   ├── barcodes.tsv.gz                    # Cell barcodes
    │   │   ├── features.tsv.gz                    # Gene identifiers
    │   │   └── matrix.mtx.gz                      # Expression counts (sparse)
    │   ├── filtered_feature_bc_matrix.h5           # Same data in HDF5 format (single file)
    │   ├── filtered_feature_bc_matrix.tar.gz       # Tarball of matrix directory
    │   ├── Healthy_1_filtered_barcodes.csv         # List of cell barcodes (same as above)
    │   │
    │   ├── raw_feature_bc_matrix/                  # **ALL BARCODES** - includes empty droplets
    │   │   ├── barcodes.tsv.gz                    # All barcodes (cells + empty droplets)
    │   │   ├── features.tsv.gz                    # Gene identifiers
    │   │   └── matrix.mtx.gz                      # Expression counts (sparse)
    │   ├── raw_feature_bc_matrix.h5                # Same data in HDF5 format
    │   ├── raw_feature_bc_matrix.tar.gz            # Tarball of matrix directory
    │   ├── raw_probe_bc_matrix.h5                  # For probe-based assays (Fixed RNA, not standard)
    │   │
    │   ├── Healthy_1_alignments.bam                # Aligned reads with cell barcode tags
    │   ├── Healthy_1_alignments.bam.bai            # BAM index file
    │   ├── Healthy_1_molecule_info.h5              # Molecule-level information (UMI counts per gene per cell)
    │   │
    │   ├── Healthy_1_cloupe.cloupe                 # For Loupe Browser visualization
    │   ├── probe_set.csv                           # Probe sequences (Fixed RNA assays only)
    │   │
    │   ├── analysis.tar.gz                         # Compressed automated analysis results
    │   │   └── (unpacked contains:)
    │   │       ├── clustering/                     # Graph-based clustering results
    │   │       ├── diffexp/                        # Differential expression between clusters
    │   │       ├── pca/                            # Principal component analysis
    │   │       ├── tsne/                           # t-SNE dimensionality reduction
    │   │       └── umap/                           # UMAP dimensionality reduction
    │   │
    │   └── cell_types.tar.gz                       # Automated cell type annotation (if enabled)
    │
    └── _log                                        # Detailed processing logs
```

**Detailed Explanation of Each Output File:**

**Primary Analysis Files (Use These):**

**1\. `web_summary.html` START HERE**

- Interactive HTML report with all QC metrics
- Open in any web browser
- Contains: barcode rank plots, t-SNE, sequencing saturation, mapping statistics
- **Action:** Always review this first to assess data quality

**2\. `filtered_feature_bc_matrix/` USE FOR DOWNSTREAM ANALYSIS**

- Gene expression matrix containing ONLY high-quality cells
- Cell Ranger’s cell calling algorithm removes empty droplets
- Three files work together:
- `barcodes.tsv.gz`: Cell identifiers (one per line)
- `features.tsv.gz`: Gene annotations (ENSEMBL ID, symbol, type)
- `matrix.mtx.gz`: Sparse matrix in Market Matrix format
- Typical size: 3,000-8,000 cells × 30,000-35,000 genes
- **When to use:** All downstream analysis (Seurat, Scanpy, Monocle, etc.)

**3\. `filtered_feature_bc_matrix.h5`**

- Same content as `filtered_feature_bc_matrix/` but in HDF5 format
- Single binary file (easier to transfer)
- Can be read directly by Seurat (`Read10X_h5()`) or Scanpy
- **When to use:** More convenient than three separate files

**Raw Data Files (Advanced Users):**

**4\. `raw_feature_bc_matrix/`**

- Contains ALL barcodes captured, including empty droplets
- Useful for custom cell calling algorithms
- Much larger than filtered matrix (100,000-500,000 barcodes)
- **When to use:**
- Developing custom QC pipelines
- Comparing different cell calling methods
- Troubleshooting low cell recovery

**5\. `raw_feature_bc_matrix.h5`**

- HDF5 version of raw matrix
- Same content, single file format

**Alignment Files:**

**6\. `Healthy_1_alignments.bam`**

- Position-sorted BAM file with cell barcode tags
- Each read tagged with `CB` (corrected barcode) and `UB` (corrected UMI)
- **Size warning:** Can be 10-50 GB per sample
- **Note:** NOT created if you used `--create-bam=false` (recommended)
- **When to use:**
- Visualizing specific loci in genome browsers (IGV, UCSC)
- Custom transcript quantification
- Investigating alignment artifacts
- **When NOT needed:** Standard single-cell analysis (clustering, differential expression, etc.)

**7\. `Healthy_1_alignments.bam.bai`**

- Index file for BAM (required for random access)
- Generated automatically by Cell Ranger (only if BAM was created)
- **Note:** NOT created if you used `--create-bam=false`

**Molecule-Level Data:**

**8\. `Healthy_1_molecule_info.h5`**

- UMI counts for each gene in each cell
- More detailed than the count matrix
- Contains per-molecule information including gem group, barcode, UMI, gene, reads
- **When to use:**
- Aggregating multiple samples
- Recalling cells with different parameters
- Detailed QC investigations

**Visualization Files:**

**9\. `Healthy_1_cloupe.cloupe`**

- For 10x Genomics’ Loupe Browser software
- Interactive exploration of clusters, gene expression, cell types
- Free software download from 10x Genomics website
- **When to use:** Quick interactive data exploration without coding

**Automated Analysis Results:**

**10\. `analysis`**

- Contains Cell Ranger’s automated analysis
- Unpacks to multiple subdirectories:
- `clustering/`: Graph-based and k-means clustering results
- `diffexp/`: Differentially expressed genes between clusters
- `pca/`: Principal component analysis (top 10 PCs)
- `tsne/`: 2D t-SNE projection (deprecated in favor of UMAP)
- `umap/`: 2D UMAP projection
- **When to use:** Quick exploratory analysis before custom analysis
- **Limitation:** Uses default parameters; custom analysis usually preferred

**11\. `cell_types.tar.gz`**

- Automated cell type annotation
- Compares cells to reference transcriptomes
- **Availability:** Only in newer Cell Ranger versions with cell typing enabled
- **Reliability:** Use as starting point; always validate with markers

**Metadata Files:**

**12\. `metrics_summary.csv`**

- All metrics from web\_summary.html in CSV format
- Easy to parse programmatically
- Use for comparing metrics across many samples

**Which Files Should You Keep?**

**Essential (always keep):**

- `web_summary.html` – Quick QC reference
- `filtered_feature_bc_matrix/` – Gene expression matrix
- `filtered_feature_bc_matrix.h5` – For all downstream analysis
- `metrics_summary.csv` – Programmatic access to QC metrics

**Important (keep if space allows):**

- `molecule_info.h5` – Allows reanalysis with different parameters
- `cloupe.cloupe` – Interactive exploration

**Optional (can delete after QC):**

- `alignments.bam` (.bai) – Very large; delete unless needed for visualization
- **Note:** Won’t exist if you used `--create-bam=false` (recommended)
- `raw_feature_bc_matrix/` – Only needed for custom cell calling

### Interpreting Cell Ranger Metrics

#### Quality Assessment Summary for Sample Healthy\_1

<iframe src="https://ngs101.com/wp-content/uploads/2025/12/web_summary.html" width="600" height="800"></iframe>

#### Detailed Metric Interpretation

Let’s break down what each metric tells us about data quality:

**1\. Estimated Number of Cells: 9,726**

- **Interpretation:** Excellent! We expected ~5,000 cells and recovered 9,726
- **Why this matters:** More cells = better statistical power for rare cell type detection
- **Red flag if:** <50% of expected (e.g., <2,500 if expecting 5,000)
- **This sample:** Nearly 2x expected – high-quality cell capture

**2\. Mean Reads per Cell: 39,180**

- **Interpretation:** Excellent sequencing depth
- **Why this matters:** More reads = more genes detected per cell
- **Minimum acceptable:** 10,000 reads/cell
- **Recommended:** 20,000-50,000 reads/cell
- **This sample:** In the sweet spot for comprehensive gene detection

**3\. Median Genes per Cell: 2,689**

- **Interpretation:** Good gene capture for PBMCs
- **Why this matters:** Reflects cell quality and sequencing depth
- **Expected for PBMCs:** 1,500-3,000 genes
- **Expected for other tissues:**
	- Brain: 2,000-4,000 genes
		- Liver: 1,000-2,500 genes
		- Tumor: Highly variable (1,000-5,000)
- **Red flag if:** <500 genes (poor quality cells or low depth)
- **This sample:** Right in the expected range for healthy PBMCs

**4\. Median UMI Counts per Cell: 7,652**

- **Interpretation:** Healthy UMI count
- **Why this matters:** UMIs represent unique molecules (not PCR duplicates)
- **Relationship to genes:** UMI/Gene ratio ≈ 3:1 is typical (this sample: 7,652/2,689 = 2.8)
- **Red flag if:** UMI count similar to gene count (indicates saturation issues)

**5\. Total Genes Detected: 30,617**

- **Interpretation:** Detected nearly all expressed genes in the transcriptome
- **Why this matters:** Shows comprehensive transcriptome coverage
- **Human genome:** ~33,000 protein-coding genes + lncRNAs
- **This sample:** Captured ~93% of human genes across all cells

**6\. Fraction Reads in Cells: 94.4%**

- **Interpretation:** Excellent – very little background
- **Why this matters:** High percentage = clean GEM formation with minimal empty droplets
- **Good:** >70%
- **Acceptable:** 50-70%
- **Poor:** <50% (indicates GEM formation issues or cell loss)
- **This sample:** Outstanding purity – 94.4% of sequencing went to real cells

**7\. Sequencing Saturation: 70.5%**

- **Interpretation:** Good saturation level
- **Why this matters:** Indicates whether sequencing depth is sufficient
- **What it means:** 70.5% of reads are duplicates of molecules already sequenced
- **Interpretation guide:**
	- <40%: Under-sequenced (sequence deeper to detect more genes)
		- 40-60%: Adequate (more sequencing helps moderately)
		- 60-80%: Good (diminishing returns from more sequencing)
		- \> 80%: Saturated (more sequencing adds little value)
- **This sample:** Well-sequenced; additional sequencing would provide minimal benefit

**8\. Valid Barcodes: 98.4%**

- **Interpretation:** Excellent barcode quality
- **Why this matters:** Invalid barcodes are wasted reads
- **Good:** >90%
- **Red flag if:** <80% (sequencing quality issues or wrong chemistry specified)
- **This sample:** Nearly perfect barcode reads

**9\. Q30 Bases in Barcode: 96.8%**

- **Interpretation:** High-quality barcode sequencing
- **Why this matters:** Low-quality barcodes → cell assignment errors
- **Q30 means:** 99.9% base call accuracy (1 error per 1,000 bases)
- **Required:** >75% for reliable cell calling
- **Good:** >90%
- **This sample:** Excellent – very few sequencing errors in barcodes

**10\. Q30 Bases in UMI: 96.7%**

- **Interpretation:** High-quality UMI sequencing
- **Why this matters:** Low-quality UMIs → PCR duplicate detection errors
- **This sample:** Nearly perfect UMI reads

**11\. Q30 Bases in RNA Read: 90.7%**

- **Interpretation:** Excellent read quality
- **Why this matters:** Low-quality reads → mapping errors
- **Required:** >65%
- **Good:** >75%
- **This sample:** Outstanding read quality

**12\. Reads Mapped Confidently to Transcriptome: 84.0%**

- **Interpretation:** Very good transcriptome mapping
- **Why this matters:** High mapping = correct species/genome + good RNA quality
- **Good:** >70%
- **Acceptable:** 50-70%
- **Red flag if:** <50% (wrong genome, contamination, or degraded RNA)
- **This sample:** Excellent mapping rate

**13\. Exonic vs Intronic Reads: 65.3% vs 26.8%**

- **Interpretation:** Normal distribution for 10x v3 with `--include-introns=true`
- **Why this matters:**
	- **Exonic reads:** From mature mRNA (standard gene expression)
		- **Intronic reads:** From pre-mRNA (nuclear RNA, useful for RNA velocity)
- **Expected ratio:** 2-3:1 exonic:intronic
- **This sample:** 65.3/26.8 = 2.4:1 (typical)
- **Note:** Cell Ranger counts both by default; intronic reads help with lowly-expressed genes

**14\. Reads Mapped Antisense: 7.5%**

- **Interpretation:** Acceptable level
- **Why this matters:** Antisense reads may indicate:
	- Natural antisense transcripts (normal biology)
		- Non-strand-specific artifacts
		- Index hopping on patterned flow cells
- **Good:** <10%
- **This sample:** Within normal range

**Conclusion for Healthy\_1:** This is **high-quality single-cell RNA-seq data**. All metrics exceed the recommended thresholds. The sample is suitable for downstream analysis including clustering, differential expression, and cell type annotation.

Use this table to evaluate your own samples:

| Metric | Good | Acceptable | Poor | Healthy\_1 |
| --- | --- | --- | --- | --- |
| **Estimated Cells** | Match expectation ±30% | ±50% | \>2x off | 9,726 (expected ~5,000) |
| **Mean Reads per Cell** | \>20,000 | 10,000-20,000 | <10,000 | 39,180 |
| **Median Genes per Cell** | \>1,500 | 1,000-1,500 | <1,000 | 2,689 |
| **Sequencing Saturation** | \>60% | 40-60% | <40% | 70.5% |
| **Reads in Cells** | \>70% | 50-70% | <50% | 94.4% |
| **Valid Barcodes** | \>90% | 80-90% | <80% | 98.4% |
| **Q30 in Barcode** | \>90% | 80-90% | <80% | 96.8% |
| **Q30 in UMI** | \>90% | 80-90% | <80% | 96.7% |
| **Q30 in RNA Read** | \>75% | 65-75% | <65% | 90.7% |
| **Mapped to Transcriptome** | \>70% | 50-70% | <50% | 84.0% |

#### When to Be Concerned

**Red Flags that Require Investigation:**

1. **<50% of expected cells**
	- Possible causes: Cell clumping, dead cells, GEM generation failure
		- Action: Review experimental protocol, check cell viability
2. **<10,000 reads per cell**
	- Possible causes: Under-sequencing, too many cells loaded
		- Action: Sequence deeper or load fewer cells next time
3. **<1,000 genes per cell**
	- Possible causes: Low sequencing depth, poor RNA quality, wrong cell type
		- Action: Check RNA quality metrics, consider re-sequencing
4. **<50% reads in cells**
	- Possible causes: Empty droplets, ambient RNA contamination
		- Action: Review GEM generation protocol, check cell concentration
5. **<80% valid barcodes**
	- Possible causes: Wrong chemistry specified, sequencing quality issues
		- Action: Verify chemistry parameter, check sequencing quality
6. **<70% mapped to transcriptome**
	- Possible causes: Wrong genome, contamination, degraded RNA
		- Action: Verify species/genome version, check for contamination

## Understanding the Output: Gene Expression Matrices

Let’s explore the count matrices Cell Ranger generated.

### Matrix Format: Market Matrix (MTX)

Cell Ranger outputs sparse matrices in Market Matrix format to save disk space:

```bash
#-----------------------------------------------
# STEP 11: Explore count matrices
#-----------------------------------------------
 
cd~/GSE174609_scRNA/cellranger_output/Healthy_1/outs/filtered_feature_bc_matrix
 
# Look at barcodes (cell IDs)
zcat barcodes.tsv.gz | head-5
# Example output:
# AAACCTGAGAAACCAT-1
# AAACCTGAGAAACCTA-1
# AAACCTGAGAAACGAG-1
 
# Look at features (genes)
zcat features.tsv.gz | head-5
# Example output:
# ENSG00000243485  MIR1302-2HG  Gene Expression
# ENSG00000237613  FAM138A      Gene Expression
# ENSG00000186092  OR4F5        Gene Expression
 
# Look at matrix structure
zcat matrix.mtx.gz | head-10
# Format: gene_index cell_index count
# Example:
# %%MatrixMarket matrix coordinate integer general
# 33538 5234 18726410
# 1 1 4
# 2 1 2
# ...
```

**Understanding the Files:**

**barcodes.tsv.gz:**

- One line per cell
- Format: `BARCODE-1` (the -1 is the GEM well)
- These are the “real” cells Cell Ranger identified

**features.tsv.gz:**

- One line per gene
- Columns: ENSEMBL\_ID, Gene\_Symbol, Feature\_Type
- Typically 30,000-35,000 genes (all annotated genes)

**matrix.mtx.gz:**

- Sparse format: Only non-zero values stored
- Three columns: gene\_index, cell\_index, count
- Saves massive amounts of disk space
- Dense matrix (all zeros included): ~4 GB
- Sparse matrix: ~100-200 MB

## Best Practices for Single-Cell RNA-seq Analysis

Based on extensive experience with 10x data, here are key recommendations:

### Experimental Design

**Sample Size:**

- **Minimum:** 2-3 biological replicates per condition
- **Recommended:** 3-4 replicates for robust statistics
- **Our study:** 4 replicates per group (excellent design)

**Sequencing Depth:**

- **Target:** 20,000-50,000 reads per cell
- **Minimum:** 10,000 reads per cell
- **More is not always better:** Diminishing returns above 50K
- **Saturation matters:** Check sequencing saturation in web summary

**Cell Numbers:**

- **PBMCs:** Target 5,000-10,000 cells per sample
- **Tissue:** May need more (10,000-20,000) due to heterogeneity
- **Rare populations:** Increase cell numbers to capture rare types

### Quality Control Checkpoints

**Before Sequencing:**

- Cell viability >80% (measure with Trypan blue or flow cytometry)
- Remove dead cells if viability <75%
- Accurate cell counting (automated counters recommended)
- Appropriate dilution (700-1,200 cells/μL for 10x)

**After Sequencing:**

- Check index read quality (Q30 >95%)
- Verify even distribution across samples
- Check for index hopping (<5% in dual-indexed runs)

**After Cell Ranger:**

- Review all web summaries carefully
- Compare across samples for consistency
- Flag any samples with <50% reads in cells
- Check for excessive mitochondrial reads (>10% may indicate stress)

### Common Issues and Solutions

**Issue 1: Low Cell Recovery (<50% of expected)**

**Possible Causes:**

- Cell clumping (aggregate cells)
- Too many dead cells
- Incorrect cell concentration
- GEM generation failure

**Solutions:**

- Filter through 40 μm strainer before loading
- Use dead cell removal kit
- Re-count cells and adjust concentration
- Check 10x instrument for clogs

**Issue 2: High Mitochondrial Content (>10%)**

**Possible Causes:**

- Cell stress during isolation
- Apoptotic cells
- Tissue type (some tissues naturally higher)

**Solutions:**

- Optimize isolation protocol (faster, colder)
- Process samples immediately after collection
- Filter mitochondrial genes in downstream analysis
- Consider nuclei isolation instead of whole cells

**Issue 3: Low Genes per Cell (<500)**

**Possible Causes:**

- Poor sequencing depth
- Low-quality RNA
- Cell type with naturally low complexity

**Solutions:**

- Increase sequencing depth
- Improve RNA isolation protocol
- Check if expected for cell type (RBCs have low complexity)

**Issue 4: Doublet Rate Higher Than Expected (>10%)**

**Possible Causes:**

- Loaded too many cells (>10,000 target)
- Cell clumping

**Solutions:**

- Reduce target cell recovery
- Filter samples more stringently
- Use computational doublet detection (DoubletFinder, Scrublet)

### Data Management

**Storage Requirements:**

- **Raw FASTQs:** Keep compressed, ~5-10 GB per sample
- **Cell Ranger output (with `--create-bam=false`):** ~1-4 GB per sample ✓ **Recommended**
- **Cell Ranger output (with BAM files):** ~11-54 GB per sample
- **BAM files alone:** ~10-50 GB per sample (only needed for genome browser visualization)
- **Downstream analysis:** ~1-5 GB per project

**Best Practice:** Use `--create-bam=false` in your Cell Ranger command to save **90% disk space**. BAM files are rarely needed for standard single-cell analysis and can always be regenerated later if required for specific visualization tasks.

**File Organization:**

```bash
project/
├── raw_data/           # Original FASTQs (archive after processing)
├── cellranger_output/  # Cell Ranger results (keep indefinitely)
├── downstream_analysis/# Seurat/Scanpy objects (primary analysis)
├── figures/            # Publication-quality figures
├── scripts/            # All analysis scripts (version control!)
└── metadata/           # Sample information, experimental design
```

**Backup Strategy:**

- Keep raw FASTQs backed up (can regenerate everything else)
- Archive Cell Ranger output after QC
- Version control all scripts (use Git)
- Document parameter choices

## Additional Resources

- **10x Genomics Support:** https://support.10xgenomics.com/
- **Cell Ranger Documentation:** https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- **Loupe Browser:** https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser

## References

1. Zheng GX, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun. 2017;8:14049.
2. Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-3587.e29.
3. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15.
4. Lee H, Joo JY, Sohn DH, et al. Single-cell RNA sequencing reveals rebalancing of immunological response in patients with periodontitis after non-surgical periodontal therapy. J Transl Med. 2022;20(1):504.
5. Stuart T, Butler A, Hoffman P, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019;177(7):1888-1902.e21.
6. Wang S, Sun S-T, Zhang X-Y, et al. The Evolution of Single-Cell RNA Sequencing Technology and Application: Progress and Perspectives. *International Journal of Molecular Sciences*. 2023; 24(3):2943.
7. Swaminath S, Russell AB (2024) The use of single-cell RNA-seq to study heterogeneity at varying levels of virus–host interactions. PLoS Pathog 20(1): e1011898.

---

**This tutorial is part of the NGS101.com series on single cell sequencing analysis. If this tutorial helped advance your research, **please comment and share** **your experience** to help other researchers! **Subscribe to stay updated** with our latest bioinformatics tutorials and resources.**