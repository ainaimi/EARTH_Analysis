#!/bin/bash
#SBATCH --job-name=afc_sens          # Job name
#SBATCH --output=logs/sens_%A_%a.out # Standard output (%A=job ID, %a=array task ID)
#SBATCH --error=logs/sens_%A_%a.err  # Standard error
#SBATCH --partition=naimi            # My Partition on the cluster
#SBATCH --array=1-6                  # Job array: 6 jobs for 24 scenarios
#SBATCH --ntasks=1                   # Number of tasks (1 per job)
#SBATCH --cpus-per-task=8            # CPUs for parallel processing within each scenario
#SBATCH --mem=64G                    # Memory per job
#SBATCH --time=6:00:00               # Max runtime (6 hours per job)

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
module load gcc/11.2.0  # Some R packages need a compiler

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
Rscript --no-save --no-restore --verbose code/5_cluster.R $SLURM_ARRAY_TASK_ID

# Print completion information
echo "========================================="
echo "Finished at: $(date)"
echo "========================================="
