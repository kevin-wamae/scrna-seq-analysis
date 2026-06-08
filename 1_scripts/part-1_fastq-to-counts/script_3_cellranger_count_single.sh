#!/bin/bash

# ********************************************************
#--- Start of slurm commands ---
# ********************************************************

#SBATCH -J cellranger_count
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G # recommended for human genome
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kevin_wamae@brown.edu

# ********************************************************
#--- End of slurm commands ---
# ********************************************************

# Input Configurations & Path Variables
#----------------------------------------------- 
FASTQ_DIR=/users/kwamae/scratch/scRNA/raw_data/
TRANSCRIPTOME=/users/kwamae/scratch/scRNA/refdata-gex-GRCh38-2024-A
SAMPLE="Healthy_1"
EXPECTED_CELLS=5000

# Slurm Resource Auto-Mapping
#----------------------------------------------- 
CORES=16 # pulls cores from --cpus-per-task
MEM=120  # pulls from --mem and converts to integer GB for Cell Ranger

# High-Speed Output Scratch Space
#----------------------------------------------- 
# Cell Ranger writes thousands of temporary files; scratch is critical here.
SCRATCH_OUTPUT_DIR="/users/kwamae/scratch/scRNA/cellranger_runs/job_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_OUTPUT_DIR"

# Move directly inside high-speed scratch space to run the analysis
cd "$SCRATCH_OUTPUT_DIR"

#-----------------------------------------------
# STEP 8: Run Cell Ranger count on single sample
#-----------------------------------------------
echo "Starting Cell Ranger count inside high-speed scratch..."
echo "Allocated Cores: $CORES, Allocated Memory: ${MEM}G"

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

echo "Cell Ranger run complete."
echo "Your final output directory is located here: ${SCRATCH_OUTPUT_DIR}/${SAMPLE}/outs/"
