# EARTH Treatment Heterogeneity Analysis

This repository contains the complete analysis pipeline for examining whether
and how environmental endocrine disrupting chemicals modify the well established
association between age and antral follicle count (a marker of fertility) in 775 women
aged 21 to 46 years.

The analysis uses data from the EARTH (Environment and Reproductive Health) Study and
employs the doubly robust (DR) learner, a method to estimate conditional
average treatment effects. We use this learner to estimate how the association between age and antral follicle count changes as a function of environmental endocrine disrupting chemicals.

## Analysis Pipeline

The analysis consists of two main components: **local analyses** (Steps 0-4) and **cluster-based sensitivity analyses** (Step 5).

### Local Analyses (Steps 0-4)

These core analyses are designed to run on a local machine and can be executed via the included Makefile.

0. **`0_data_gen.R`** - Data generation and random forest imputation
   - Loads raw EARTH study data
   - Analyzes missing data patterns
   - Performs random forest imputation for variables with <40% missingness
   - Excludes variables with >40% missingness (BPS, BPF, BP3, TCS, OPFR metabolites)
   - Generates imputation diagnostics and comparison plots

1. **`1_data_man.R`** - Data management and preprocessing
   - Loads imputed data from Step 0
   - Creates cleaned datasets with and without truncation of extreme values
   - Prepares analysis-ready datasets

2. **`2_chem_analysis.R`** - Chemical modifier exploration and outlier detection
   - Examines distributions of environmental modifiers
   - Compares truncated vs. non-truncated modifier data
   - Generates descriptive statistics
   - Performs multivariate outlier detection (Mahalanobis, PCA, LOF, Cook's D)

3. **`3_IF_scores_gen.R`** - Influence function score generation
   - Implements the DR learner models using cross-validated SuperLearner
   - Generates influence function (IF) scores for statistical inference
   - Produces variable importance measures and propensity score diagnostics

4. **`4_IF_scores_analysis.R`** - Heterogeneity analysis and visualization
   - Analyzes estimated treatment effect heterogeneity using best linear projection
   - Performs hypothesis tests for effect modification (unconditional and conditional models)
   - Creates CATE function plots (linear and SuperLearner-based)
   - Generates manuscript figures

**Quarto Report** - `EARTH_Analysis_Report.qmd`
   - Compiles all results into a comprehensive HTML report
   - Includes missing data analysis, descriptive statistics, model diagnostics, and CATE results

### Cluster-Based Sensitivity Analyses (Step 5)

Computationally intensive sensitivity analyses testing robustness across 24 scenarios:
- Age thresholds: 35, 36, 37, 38, 39, 40 years
- AFC>30 imputation (yes/no)
- EDC>97.5th percentile imputation (yes/no)

**Cluster Scripts:**
- **`cluster/5_cluster.R`** - Main sensitivity analysis script for HPC execution
- **`cluster/earth.sh`** - SLURM job submission script (6 parallel jobs, 4 scenarios each)
- **`cluster/5_combine_results.R`** - Combines results from all cluster jobs
- **`cluster/5_sensitivity_plot.R`** - Generates sensitivity analysis visualizations

See `code/cluster/CLUSTER_README.md` for detailed instructions on running sensitivity analyses on an HPC cluster with SLURM.

## Requirements

### R Version
- R >= 4.0.0 recommended

## Project Structure

```
EARTH_Analysis/
├── code/                      # Analysis scripts
│   ├── 0_data_gen.R           # Data generation and imputation
│   ├── 1_data_man.R           # Data management
│   ├── 2_chem_analysis.R      # Chemical analysis and outliers
│   ├── 3_IF_scores_gen.R      # IF score generation
│   ├── 4_IF_scores_analysis.R # Heterogeneity analysis
│   ├── cluster/               # Cluster-based sensitivity analysis
│   │   ├── 5_cluster.R        # Main cluster analysis script
│   │   ├── 5_combine_results.R # Combine cluster results
│   │   ├── 5_sensitivity_plot.R # Sensitivity plots
│   │   ├── earth.sh           # SLURM submission script
│   │   └── CLUSTER_README.md  # Cluster usage guide
│   ├── EARTH_Analysis_Report.qmd # Quarto report
│   └── .gitignore             # Ignores Quarto cache and temp files
├── data/                      # Raw and processed data files (not tracked)
├── figures/                   # Generated figures (not tracked)
├── output/                    # Analysis results as .rds files (not tracked)
├── logs/                      # Cluster job logs (not tracked)
├── manuscript/                # Manuscript drafts (not tracked)
├── sandbox/                   # Exploratory analyses (not tracked)
├── misc/                      # Miscellaneous files (not tracked)
├── Makefile                   # Automated local pipeline orchestration
└── README.md                  # This file
```

## Usage

### Setup
1. Ensure all required R packages are installed
2. Place raw EARTH data files in the `data/` directory
3. Open the R project: `EARTH_Analysis.Rproj`

### Running Local Analyses

**Using Make (Recommended)**
```bash
# Run entire local pipeline (Steps 0-4 + report generation)
make all

# Or run individual steps
make step0  # Data generation and imputation
make step1  # Data management
make step2  # Chemical analysis
make step3  # IF score generation
make step4  # Heterogeneity analysis
make report # Generate HTML report only

# Clean outputs
make clean        # Remove intermediate marker files
make clean-report # Remove report only
make clean-all    # Remove ALL generated outputs (use with caution!)
```

**Manual Execution**
```bash
# Run scripts sequentially in R
Rscript code/0_data_gen.R
Rscript code/1_data_man.R
Rscript code/2_chem_analysis.R
Rscript code/3_IF_scores_gen.R
Rscript code/4_IF_scores_analysis.R

# Generate report
quarto render code/EARTH_Analysis_Report.qmd --output-dir output
```

### Running Cluster-Based Sensitivity Analyses

For detailed instructions on running sensitivity analyses on an HPC cluster, see `code/cluster/CLUSTER_README.md`.

**Quick Start:**
```bash
# On cluster: Submit sensitivity analysis jobs
sbatch code/cluster/earth.sh

# After jobs complete: Combine results
Rscript code/cluster/5_combine_results.R

# Generate sensitivity plots
Rscript code/cluster/5_sensitivity_plot.R
```

## Citation

coming soon.

## License

This project is licensed under the MIT License - see below for details:

```
MIT License

Copyright (c) 2025 [Your Name/Institution]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

**Note:** The data from the EARTH Study is not publicly available due to 
privacy restrictions and data use agreements. Access to the data 
requires approval from the EARTH Study team.
