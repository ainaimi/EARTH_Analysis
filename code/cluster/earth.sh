#!/usr/bin/env bash
#SBATCH --job-name=afc_sens
#SBATCH --partition=naimi
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=30-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anaimi@emory.edu

# ==============================================================================
# SLURM Submission Script for AFC Sensitivity Analysis
# ==============================================================================
# This script runs 24 sensitivity analysis scenarios across 6 jobs.
# Each job processes 4 scenarios sequentially.
#
# Usage:
#   sbatch code/earth.sh
#
# To check job status:
#   squeue -u $USER
#
# To cancel jobs:
#   scancel <job_id>
# ==============================================================================

# Print job information
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Started at: $(date)"
echo "========================================="

# Load required modules (adjust for your cluster)
module purge
module load R/4.4.0

# Set number of threads for OpenBLAS/MKL (if applicable)
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to project directory
cd /home/anaimi/EARTH  

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the R script with the array task ID
# This tells R which subset of scenarios to run
Rscript --no-save --no-restore --verbose ./code/5_cluster.R $SLURM_ARRAY_TASK_ID > ./output/cluster_output_${SLURM_ARRAY_TASK_ID}.Rout 2>&1

# Print completion information
echo "========================================="
echo "Finished at: $(date)"
echo "========================================="
