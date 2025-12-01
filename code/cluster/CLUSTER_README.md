# Running Sensitivity Analysis on HPC Cluster

This guide explains how to run the AFC sensitivity analysis on an HPC cluster using SLURM.

## Overview

The sensitivity analysis tests 24 scenarios (6 age thresholds × 2 AFC imputation options × 2 EDC imputation options). The cluster workflow splits these across 6 jobs, with each job processing 4 scenarios sequentially.

## Files

- `earth.sh` - SLURM submission script
- `5_cluster.R` - Main R script for cluster execution
- `5_combine_results.R` - Script to combine results after jobs complete

## Setup Instructions

### 1. Prepare Your Environment

Before submitting jobs, you need to update the cluster-specific settings:

**Edit `earth.sh`:**
- Line 10: Update `--partition=standard` to match your cluster's partition name
- Line 41: Update `/path/to/your/project` to your actual project path
- Lines 35-36: Uncomment and adjust module loads if needed (e.g., `module load R/4.3.0`)

### 2. Create Required Directories

On the cluster, navigate to your project directory and create:

```bash
mkdir -p logs
```

This directory will store job output and error logs.

### 3. Make Scripts Executable

```bash
chmod +x code/earth.sh
chmod +x code/5_cluster.R
chmod +x code/5_combine_results.R
```

### 4. Transfer Files to Cluster

Transfer the following to your cluster:
- `data/imputed_EARTH.Rdata`
- All R scripts in `code/`
- Any required package installations

## Running the Analysis

### Step 1: Submit Jobs

From your project directory on the cluster:

```bash
sbatch code/earth.sh
```

This will submit 6 jobs as an array (job IDs 1-6).

### Step 2: Monitor Jobs

Check job status:
```bash
squeue -u $USER
```

Check specific job details:
```bash
scontrol show job <job_id>
```

View live output (replace JOBID and TASKID):
```bash
tail -f logs/sens_JOBID_TASKID.out
```

### Step 3: Check for Errors

After jobs complete, check for errors:
```bash
grep -i error logs/sens_*.err
```

### Step 4: Combine Results

Once all 6 jobs complete successfully:

```bash
Rscript code/5_combine_results.R
```

This will:
- Find all `sensitivity_task_*.rds` files in `output/`
- Combine them into a single file: `output/sensitivity_analysis_results.rds`
- Report summary statistics
- Verify all 24 scenarios completed

### Step 5: Generate Plots

After combining results, generate plots:

```bash
Rscript code/5c. AFC_sensitivity_plot.R
```

Or transfer `output/sensitivity_analysis_results.rds` back to your local machine and run the plotting script there.

## Troubleshooting

### Jobs Fail with Memory Errors
- Increase `--mem=64G` in `earth.sh` to `--mem=128G`
- Or reduce `--cpus-per-task` from 8 to 4

### Jobs Timeout
- Increase `--time=3:00:00` in `earth.sh` to `--time=6:00:00`

### Missing Scenarios After Completion
- Check error logs: `cat logs/sens_*.err`
- Identify which task failed
- Resubmit specific array task: `sbatch --array=<failed_task_id> code/earth.sh`

### Module Load Errors
- Check available R modules: `module avail R`
- Update module load commands in `earth.sh`

### Package Installation Issues
On the cluster, you may need to install R packages in your home directory:

```r
# In R on cluster:
install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load_gh("ecpolley/SuperLearner")
pacman::p_load(tidyverse, broom, lmtest, sandwich, ranger, xgboost, missForest)
```

## Resource Usage

Per job (4 scenarios each):
- **Runtime:** ~1.5-2.5 hours (depends on data size and SuperLearner complexity)
- **Memory:** 64GB (adjust if needed)
- **CPUs:** 8 cores
- **Total cluster usage:** 6 jobs × 2 hours × 8 cores = 96 core-hours

## Output Files

After successful completion:
- `output/sensitivity_task_1.rds` through `output/sensitivity_task_6.rds` (individual task results)
- `output/sensitivity_analysis_results.rds` (combined results)
- `logs/sens_JOBID_*.out` (standard output logs)
- `logs/sens_JOBID_*.err` (error logs)

## Canceling Jobs

Cancel all array jobs:
```bash
scancel <job_id>
```

Cancel specific array task:
```bash
scancel <job_id>_<array_task_id>
```

## Cluster-Specific Notes

### For Different Schedulers

If your cluster uses PBS/Torque or SGE instead of SLURM, you'll need to adapt `earth.sh`:

**PBS/Torque:**
```bash
#PBS -N afc_sens
#PBS -t 1-6
#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=3:00:00
```

**SGE:**
```bash
#$ -N afc_sens
#$ -t 1-6
#$ -pe smp 8
#$ -l h_vmem=64G
#$ -l h_rt=3:00:00
```

### For Different Numbers of Jobs

To change from 6 jobs to a different number, update:
1. `earth.sh`: Line 5 `--array=1-N` (where N = number of jobs)
2. `5_cluster.R`: Line 54 `scenarios_per_job` (should be 24/N)

Example for 8 jobs (3 scenarios each):
- `earth.sh`: `--array=1-8`
- `5_cluster.R`: `scenarios_per_job <- 3`
