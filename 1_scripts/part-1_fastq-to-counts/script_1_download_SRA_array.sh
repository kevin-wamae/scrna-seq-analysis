#!/bin/bash
#*********************************************************
#SBATCH -J SRA-bulk-download
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH -o logs/job.%A_%a.out
#SBATCH -e logs/job.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=''
#SBATCH --array=1-12%4
#*********************************************************

# Make pixi findable on compute nodes
export PATH="$HOME/.pixi/bin:$PATH"
eval "$(pixi shell-hook)"

# --- Path Configuration ---
PROJECT_DIR="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis"
SCRATCH_DIR="/users/kwamae/scratch/scRNA/raw_data/sra_cache"
ACC_LIST="/users/kwamae/scratch/scRNA/raw_data/SRR_list.txt"

# Create output folders
mkdir -p logs "$SCRATCH_DIR/fastqs" "$SCRATCH_DIR/sra_cache"

# Extract the accession for this array task
SRR=$(awk "NR==$SLURM_ARRAY_TASK_ID {print \$1}" "$ACC_LIST")

if [ -z "$SRR" ]; then
    echo "Error: No accession found at line $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing task $SLURM_ARRAY_TASK_ID: $SRR"

# --- Step 1: Download ---
pixi run \
    --manifest-path "$PROJECT_DIR/pixi.toml" \
    -e part1 \
    prefetch "$SRR" \
    --output-directory "$SCRATCH_DIR/sra_cache"

# --- Step 2: Convert ---
pixi run \
    --manifest-path "$PROJECT_DIR/pixi.toml" \
    -e part1 \
    fasterq-dump \
    --split-3 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --progress \
    --outdir "$SCRATCH_DIR/fastqs" \
    "$SCRATCH_DIR/sra_cache/${SRR}/${SRR}.sra"

# --- Step 3: Compress ---
if command -v pigz >/dev/null 2>&1; then
    pigz -p "$SLURM_CPUS_PER_TASK" "$SCRATCH_DIR/fastqs/${SRR}"*.fastq
else
    gzip "$SCRATCH_DIR/fastqs/${SRR}"*.fastq
fi

# --- Step 4: Validate ---
failed=0
for f in "$SCRATCH_DIR/fastqs/${SRR}_1.fastq.gz" "$SCRATCH_DIR/fastqs/${SRR}_2.fastq.gz"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Expected output $f is missing or empty"
        failed=1
    fi
done

if [[ $failed -eq 1 ]]; then
    echo "Skipping cleanup for $SRR — verify outputs manually"
    exit 1
fi

# --- Step 5: Cleanup ---
rm -rf "$SCRATCH_DIR/sra_cache/${SRR}"

echo "Task $SLURM_ARRAY_TASK_ID ($SRR) completed successfully."
