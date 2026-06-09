#!/bin/bash
#*********************************************************
#SBATCH -J copy-cellranger-output
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH -o logs/job.%j.out
#SBATCH -e logs/job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kevin_wamae@brown.edu
#*********************************************************

# --- Path Configuration ---
# Source: Cell Ranger output on scratch
SOURCE_DIR="/users/kwamae/scratch/scRNA/cellranger_runs/job_3058993"

# Destination: project input directory on home
DEST_DIR="/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/2_input/cellranger-matrix-counts/GSE174609_All_Participants"

# Create destination if it doesn't exist
mkdir -p "$DEST_DIR"

# --- Copy ---
# Loops over every sample folder in the job directory,
# creates a matching named folder in the destination,
# and copies only the four files/directories needed
# for downstream QC — skipping extras/ and analysis/
# which are large and not used in the pipeline.
for SAMPLE_DIR in "$SOURCE_DIR"/*/; do

    SAMPLE=$(basename "$SAMPLE_DIR")
    OUTS="$SAMPLE_DIR/outs"

    echo "Copying $SAMPLE..."

    mkdir -p "$DEST_DIR/$SAMPLE"

    cp -r "$OUTS/filtered_feature_bc_matrix" "$DEST_DIR/$SAMPLE/"
    cp -r "$OUTS/raw_feature_bc_matrix"      "$DEST_DIR/$SAMPLE/"
    cp    "$OUTS/metrics_summary.csv"         "$DEST_DIR/$SAMPLE/"
    cp    "$OUTS/web_summary.html"            "$DEST_DIR/$SAMPLE/"

    echo "Done: $DEST_DIR/$SAMPLE"

done

echo "========================================"
echo "All samples copied successfully."
echo "Destination: $DEST_DIR"
ls "$DEST_DIR"
echo "========================================"