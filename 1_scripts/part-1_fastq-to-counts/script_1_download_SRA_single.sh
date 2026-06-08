#!/bin/bash
#*********************************************************

# ********************************************************
#--- Start of slurm commands ---
# ********************************************************

#SBATCH -J SRA-download
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wamaekevin@gmail.com

# ********************************************************
#--- End of slurm commands ---
# ********************************************************

# Activate the pixi environment shell hook
eval "$(pixi shell-hook)"

# We'll download SRR14575500 (Healthy donor #1)
# This is a representative PBMC sample from the study
SRR="SRR14575500"

# Step 1: Download SRA file
# This creates a directory ~/ncbi/public/sra/SRR14575500.sra
pixi run \
    -e part1 \
    prefetch $SRR

# Step 2: Convert to FASTQ with split files
# --split-files: Create separate R1 and R2 files
# --include-technical: Include index reads if present (I1, I2)
# --threads: Use multiple cores for faster conversion
# --progress: Show progress bar
pixi run \
    -e part1 \
    fasterq-dump \
    --split-files \
    --include-technical \
    --threads 8 \
    --progress \
    --outdir . \
    $SRR/${SRR}.sra

# Step 3: Compress FASTQ files using ALL cores (reduces file size by ~70%)
# 🟢 FIX: Check if 'pigz' is inside your pixi environment, or use system pigz
if command -v pigz >/dev/null 2>&1; then
    pigz -p $SLURM_CPUS_PER_TASK ${SRR}*.fastq  # matches cores in slurm header
else
    # Fallback to standard gzip if pigz isn't installed on the cluster
    gzip ${SRR}*.fastq
fi

# Step 4: Clean up SRA file to save disk space
rm -rf $SRR