# EARTH Treatment Heterogeneity Analysis

This repository contains the complete analysis pipeline for examining whether 
and how environmental endocrine disrupting chemicals modify the well established 
association between age and antral follicle count (a marker of fertility) in 775 women 
aged 21 to 46 years.
 
The analysis uses data from the EARTH (Environment and Reproductive Health) Study and 
employs two "learners" commonly used to estimate conditional average treatment effects.
Note, however, that we use these learners to estimate associations.

## Analysis Pipeline

The analysis consists of five sequential scripts that should be run in order:

1. **`1. AFC_data_man.R`** - Data management and preprocessing
   - Loads and cleans raw EARTH study data
   - Handles missing data
   - Creates analysis-ready datasets

2. **`2. AFC_chem_analysis.R`** - Chemical modifier exploration
   - Examines distributions of environmental modifiers
   - Compares trimmed vs. untrimmed modifier data
   - Generates descriptive statistics

3. **`3. AFC_IF_scores_gen.R`** - Influence function score generation
   - Implements causal forest models using GRF
   - Implements the DR learner models using cross validated SuperLearner
   - Generates influence function (IF) scores for statistical inference for both

4. **`4. AFC_IF_scores_analysis.R`** - Heterogeneity analysis
   - Analyzes estimated treatment effect heterogeneity
   - Conducts subgroup analyses
   - Performs hypothesis tests for effect modification

5. **`5. AFC_CATE_plots.R`** - Visualization
   - Creates plots of modified associations
   - Generates figures for manuscript

## Requirements

### R Version
- R >= 4.0.0 recommended

### Required R Packages
The analysis uses `pacman` for package management and requires:

**Core packages:**
- `tidyverse` - Data manipulation and visualization
- `here` - Path management

**Causal inference:**
- `grf` - Generalized random forests for CATE estimation
- `SuperLearner` - Ensemble machine learning

**Statistical analysis:**
- `lmtest`, `sandwich` - Robust standard errors
- `pROC` - ROC curve analysis
- `broom` - Model tidying

**Machine learning:**
- `ranger` - Random forest implementation
- `xgboost` - Gradient boosting
- `polspline` - Polynomial splines

**Visualization:**
- `ggplot2` (via tidyverse)
- `gridExtra` - Multi-panel plots
- `scales` - Plot scaling
- `vip` - Variable importance plots

**Data handling:**
- `rio` - Import/export
- `haven` - SAS/Stata file reading
- `skimr` - Data summaries

**Other:**
- `xtable` - Table formatting
- `parallel` - Parallel computing

## Project Structure

```
EARTH_Analysis/
├── code/               # Analysis scripts (run in numbered order)
├── data/              # Raw and processed data files (not tracked in git)
├── figures/           # Generated figures (not tracked in git)
├── manuscript/        # Manuscript drafts (not tracked in git)
├── sandbox/           # Exploratory analyses (not tracked in git)
├── misc/              # Miscellaneous files (not tracked in git)
└── README.md          # This file
```

## Usage

### Setup
1. Ensure all required R packages are installed
2. Place raw EARTH data files in the `data/` directory
3. Open the R project: `EARTH_Analysis.Rproj`

### Running the Analysis
Execute scripts sequentially:

```r
# Run from R console with project open
source(here::here("code", "1. AFC_data_man.R"))
source(here::here("code", "2. AFC_chem_analysis.R"))
source(here::here("code", "3. AFC_IF_scores_gen.R"))
source(here::here("code", "4. AFC_IF_scores_analysis.R"))
source(here::here("code", "5. AFC_CATE_plots.R"))
```

### Output
- Cleaned datasets: `data/` directory
- Figures: `figures/` directory
- Intermediate model objects: `misc/` directory (large .RDS files)

## Data

This analysis uses data from the EARTH Study. Data files are not included in 
this repository due to privacy and data use agreements. Authorized users should 
place data files in the `data/` directory.

## Reproducibility

- All file paths use the `here` package for platform-independent path management
- Random seed setting is implemented in scripts for reproducibility
- Session info is captured in analysis outputs

## Citation

coming soon.

## Acknowledgments

This work was supported by [funding sources]. We thank the EARTH Study participants 
and staff for their contributions to this research.

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
