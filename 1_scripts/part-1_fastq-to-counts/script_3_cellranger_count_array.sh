#!/bin/bash
#*********************************************************
# Cell Ranger count — array job for all 12 samples
# Runs alignment and gene expression quantification
# against the GRCh38 human reference genome
#*********************************************************
#SBATCH -J cellranger-count
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH -o logs/job.%A_%a.out
#SBATCH -e logs/job.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kevin_wamae@brown.edu
#SBATCH --array=1-12%3
#*********************************************************

# --- Path Configuration ---
# PROJECT_DIR: where pixi.toml, sample_list.txt, and logs live
# FASTQ_DIR: raw FASTQs
# TRANSCRIPTOME: pre-built Cell Ranger GRCh38 reference from 10x genomics
# SCRATCH_OUTPUT_DIR: all Cell Ranger output goes to scratch;
#   named with the shared array job ID so all samples land
#   in one folder (e.g. cellranger_runs/job_3054746/Healthy_1/outs/)
PROJECT_DIR="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis"
FASTQ_DIR="/users/kwamae/scratch/scRNA/raw_data"
TRANSCRIPTOME="/users/kwamae/scratch/scRNA/refdata-gex-GRCh38-2024-A"
SCRATCH_OUTPUT_DIR="/users/kwamae/scratch/scRNA/cellranger_runs/job_${SLURM_ARRAY_JOB_ID}"
SAMPLE_LIST="$PROJECT_DIR/scripts/part-1_fastq-to-counts/sample_list.txt"

# --- Slurm Resource Auto-Mapping ---
# Pulling directly from SLURM variables keeps these in sync
CORES=$SLURM_CPUS_PER_TASK
MEM=154  # stay a few GB under the ceiling of $SLURM_MEM_PER_NODE

# --- Directory Setup ---
# logs/ is relative to the project directory (where you run sbatch)
# SCRATCH_OUTPUT_DIR is created fresh for this job if it doesn't exist
mkdir -p logs "$SCRATCH_OUTPUT_DIR"

# --- Sample Extraction ---
# Pulls the sample name at this task's line number from sample_list.txt
# e.g. task 1 -> Healthy_1, task 5 -> Periodontitis_Pre_1
SAMPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID {print \$1}" "$SAMPLE_LIST")

# shows an error if no sample is found
if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found at line $SLURM_ARRAY_TASK_ID in $SAMPLE_LIST"
    exit 1
fi

echo "========================================"
echo "Task:    $SLURM_ARRAY_TASK_ID / 12"
echo "Sample:  $SAMPLE"
echo "Cores:   $CORES"
echo "Memory:  ${MEM}G"
echo "Output:  $SCRATCH_OUTPUT_DIR/$SAMPLE/outs/"
echo "========================================"

# --- Move to Scratch ---
# Cell Ranger writes thousands of temporary files during alignment.
# Running from scratch avoids overwhelming shared project filesystem.
cd "$SCRATCH_OUTPUT_DIR" || { echo "ERROR: Cannot cd to $SCRATCH_OUTPUT_DIR"; exit 1; }

# --- Run Cell Ranger count ---
# --id:            output folder name, one per sample
# --fastqs:        directory containing all renamed FASTQ files
# --sample:        prefix to match the correct FASTQs for this sample
#                  (e.g. Healthy_1 matches Healthy_1_S1_L001_R1_001.fastq.gz)
# --transcriptome: pre-built GRCh38 Cell Ranger reference from 10x genomics
# --expect-cells:  expected cell count; guides Cell Ranger's cell calling
# --localcores:    cap CPU usage to our SLURM allocation
# --localmem:      cap memory usage to our SLURM allocation
# --chemistry:     3' v3 chemistry used in this study
# --create-bam:    skipped to save disk space; enable if you need
#                  per-read alignment for downstream tools
cellranger count \
    --id="$SAMPLE" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --transcriptome="$TRANSCRIPTOME" \
    --expect-cells=5000 \
    --localcores="$CORES" \
    --localmem="$MEM" \
    --chemistry=SC3Pv3 \
    --create-bam=false

# --- Validate Output ---
# Checks that the key output file exists and is non-empty before
# declaring success. Catches silent Cell Ranger failures that still
# exit with code 0.
if [ ! -f "$SCRATCH_OUTPUT_DIR/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" ]; then
    echo "ERROR: Expected Cell Ranger output missing for $SAMPLE"
    echo "Check: $SCRATCH_OUTPUT_DIR/$SAMPLE/_log"
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID ($SAMPLE) completed successfully."
echo "Output: $SCRATCH_OUTPUT_DIR/$SAMPLE/outs/"