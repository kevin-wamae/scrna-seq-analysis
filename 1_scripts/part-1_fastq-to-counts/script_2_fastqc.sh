#!/bin/bash
#*********************************************************

# ********************************************************
#--- Start of slurm commands ---
# ********************************************************

#SBATCH -J fastqc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kevin_wamae@brown.edu

# ********************************************************
#--- End of slurm commands ---
# ********************************************************

# variable names
PIXI_MANIFEST_PATH="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/pixi.toml"
RAW_DATA_PATH="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/resources/data/GSE174609/raw_data/"
FASTQC_REPORTS_PATH="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/resources/data/GSE174609/fastqc_reports/"

# Quality control with FastQC
#----------------------------------------------- 
cd $RAW_DATA_PATH
 
# # Create output directory
# mkdir -p fastqc_reports
 
# Run FastQC on all FASTQ files
# We'll run on R2 only (cDNA reads)
#----------------------------------------------- 
pixi run \
    --manifest-path /nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/pixi.toml \
    -e part1 \
    fastqc *R2_001.fastq.gz \
    --outdir $FASTQC_REPORTS_PATH \
    --threads 8 \
    --quiet
 
# Aggregate results with MultiQC
#----------------------------------------------- 
pixi run \
    --manifest-path /nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/pixi.toml \
    -e part1 \
    multiqc $FASTQC_REPORTS_PATH \
    --outdir $FASTQC_REPORTS_PATH \
    --filename multiqc_report \
    --title "GSE174609 Quality Control"