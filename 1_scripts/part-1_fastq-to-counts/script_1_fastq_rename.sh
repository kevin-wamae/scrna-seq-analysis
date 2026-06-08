#!/bin/bash
#*********************************************************
#SBATCH -J rename-fastqs
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -o logs/job.%j.out
#SBATCH -e logs/job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kevin_wamae@brown.edu
#*********************************************************

FASTQ_DIR="/oscar/scratch/kwamae/scRNA/raw_data/fastqs"

cd "$FASTQ_DIR" || { echo "ERROR: Cannot cd to $FASTQ_DIR"; exit 1; }

# --- Healthy donors ---
mv SRR14575500_1.fastq.gz Healthy_1_S1_L001_R1_001.fastq.gz
mv SRR14575500_2.fastq.gz Healthy_1_S1_L001_R2_001.fastq.gz

mv SRR14575501_1.fastq.gz Healthy_2_S2_L001_R1_001.fastq.gz
mv SRR14575501_2.fastq.gz Healthy_2_S2_L001_R2_001.fastq.gz

mv SRR14575502_1.fastq.gz Healthy_3_S3_L001_R1_001.fastq.gz
mv SRR14575502_2.fastq.gz Healthy_3_S3_L001_R2_001.fastq.gz

mv SRR14575503_1.fastq.gz Healthy_4_S4_L001_R1_001.fastq.gz
mv SRR14575503_2.fastq.gz Healthy_4_S4_L001_R2_001.fastq.gz

# --- Periodontitis pre-treatment ---
mv SRR14575504_1.fastq.gz Periodontitis_Pre_1_S5_L001_R1_001.fastq.gz
mv SRR14575504_2.fastq.gz Periodontitis_Pre_1_S5_L001_R2_001.fastq.gz

mv SRR14575505_1.fastq.gz Periodontitis_Pre_2_S6_L001_R1_001.fastq.gz
mv SRR14575505_2.fastq.gz Periodontitis_Pre_2_S6_L001_R2_001.fastq.gz

mv SRR14575506_1.fastq.gz Periodontitis_Pre_3_S7_L001_R1_001.fastq.gz
mv SRR14575506_2.fastq.gz Periodontitis_Pre_3_S7_L001_R2_001.fastq.gz

mv SRR14575507_1.fastq.gz Periodontitis_Pre_4_S8_L001_R1_001.fastq.gz
mv SRR14575507_2.fastq.gz Periodontitis_Pre_4_S8_L001_R2_001.fastq.gz

# --- Periodontitis post-treatment ---
mv SRR14575508_1.fastq.gz Periodontitis_Post_1_S9_L001_R1_001.fastq.gz
mv SRR14575508_2.fastq.gz Periodontitis_Post_1_S9_L001_R2_001.fastq.gz

mv SRR14575509_1.fastq.gz Periodontitis_Post_2_S10_L001_R1_001.fastq.gz
mv SRR14575509_2.fastq.gz Periodontitis_Post_2_S10_L001_R2_001.fastq.gz

mv SRR14575510_1.fastq.gz Periodontitis_Post_3_S11_L001_R1_001.fastq.gz
mv SRR14575510_2.fastq.gz Periodontitis_Post_3_S11_L001_R2_001.fastq.gz

mv SRR14575511_1.fastq.gz Periodontitis_Post_4_S12_L001_R1_001.fastq.gz
mv SRR14575511_2.fastq.gz Periodontitis_Post_4_S12_L001_R2_001.fastq.gz

echo "Renaming complete. Verifying..."
ls -lh *.fastq.gz