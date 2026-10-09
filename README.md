# Environmental Chemicals as Modifiers of the Association between Chronological Age and Ovarian Reserve

Analysis code for:

> Naimi AI, Kennedy EH, Yu Y-H, Hauser R, Mínguez-Alarcón L, Gaskins A.
> *Environmental Chemicals as Modifiers of the Association between Chronological Age and Ovarian Reserve.*
> Preprint: medRxiv, [doi:10.64898/2025.12.09.25341902](https://doi.org/10.64898/2025.12.09.25341902).
> Peer-reviewed version in revision at *Environmental Epidemiology*.

The repository holds the full pipeline from the raw EARTH Study extract to the
figures and tables in the manuscript, together with a rendered HTML report,
[`output/EARTH_Analysis_Report.html`](output/EARTH_Analysis_Report.html), that
serves as the online supporting material referenced in the paper. GitHub shows
the raw source of that file; to view it rendered, download it or open it through
[htmlpreview](https://htmlpreview.github.io/?https://github.com/ainaimi/EARTH_Analysis/blob/main/output/EARTH_Analysis_Report.html).

## Study question

Antral follicle count (AFC), a marker of ovarian reserve, declines with age. We
asked whether the size of the age–AFC association differs across concentrations
of 16 endocrine-disrupting chemicals (EDCs) in 775 women enrolled in the
Environment and Reproductive Health (EARTH) Study at the Massachusetts General
Hospital Fertility Center between 2004 and 2019.

- **Population.** Women aged 18–45 at enrollment who had one AFC measurement
  (day 3 of an unstimulated cycle or of a progesterone-withdrawal bleed) and at
  least one urinary or hair chemical measurement. After excluding scans done on
  Lupron, incomplete scans, scans with polycystic ovaries, and repeated scans,
  775 women contributed one scan each. Observed age at the scan ranged from 21
  to 46 years.
- **Exposure contrast.** Age at the AFC scan, dichotomized at 35 years (the
  ACOG definition of advanced maternal age): ≥35 vs <35.
- **Outcome.** AFC, the sum of antral follicles in both ovaries, used as
  recorded. AFC values were not truncated or trimmed.
- **Candidate modifiers.** Sixteen chemicals: 11 urinary phthalate metabolites
  (MBP, MiBP, MCNP, MCOP, MECPP, MEHHP, MEHP, MEOHP, MCPP, MEP, MBzP), 3 urinary
  parabens (methyl-, propyl-, and butylparaben), urinary bisphenol A, and hair
  mercury (ng/g). Urinary concentrations are specific-gravity corrected;
  values below the limit of detection were replaced by LOD/√2 in the source
  data before this pipeline. Seven other chemicals measured in EARTH (BPS,
  BPF, benzophenone-3, triclosan, and three organophosphate flame-retardant
  metabolites) were missing for more than 40% of women and were not analyzed.
- **Other covariates.** Year and month of the AFC scan, smoking status, race,
  education, BMI, prior intrauterine insemination, prior in vitro
  fertilization, gravidity, and indicator variables for imputed values.

## What is estimated, and what is not

The analysis targets two covariate-adjusted *associations*, written in the
manuscript as

- ψ = μ₁ − μ₀, with μₐ = ∫ E(Y | A = a, X) dP(X): the mean difference in AFC
  between women aged ≥35 and <35 years, standardized to the covariate
  distribution of the cohort; and
- τ(x) = E(Y | A = 1, X = x) − E(Y | A = 0, X = x): the same difference among
  women with covariate values x. The association is modified when τ(x) varies
  with the chemical concentrations in x.

These quantities are descriptive. The manuscript does not interpret them as
causal effects of age or of the chemicals, and the language in this repository
follows suit, for the reasons given in the paper:

1. Age is not a manipulable exposure, so the identification conditions usually
   invoked for average treatment effects are not justified here.
2. The data are cross-sectional, and most chemical measurements postdate the
   AFC scan: 80% of phthalate/BPA samples and 72% of hair mercury samples were
   collected after the scan, a median of about four weeks later (Step 0 of the
   report). A chemical measured after the outcome cannot be a cause of it.
3. Several covariates (BMI, gravidity, IUI/IVF history, and the chemicals
   themselves) may lie on pathways between age and AFC, so the standardized
   contrasts may be conditioned on intermediates.

The DR-learner comes from the heterogeneous-treatment-effect literature, and
the code inherits that vocabulary in a few object names (`aipw_ate`, `res_ate`,
`create_cate_plot`) and in the names of the methods it calls. Read these as
labels for the standardized associational contrasts above, not as claims about
effects.

## Methods in brief

1. **Missing data.** Random forest imputation (`missForest`, 50 iterations,
   2,000 trees, seed 123) for all retained variables, none of which was missing
   for more than 18% of women. Indicator variables flagging imputed values enter
   the nuisance models.
2. **Nuisance models.** Ten-fold cross-validated Super Learner (non-negative
   least squares meta-learner) for the outcome regression E(Y | A, X) and the
   propensity score P(A = 1 | X), using the same folds for both so that all
   predictions are cross-fitted. Library: mean, GLM, random forests (`ranger`,
   6 tuning settings), gradient boosting (`xgboost`, 4 settings), elastic net
   (`glmnet`, 6 values of α), and a GLM with all two-way interactions preceded
   by a correlation screen.
3. **Pseudo-outcome.** The augmented inverse probability weighted (AIPW) score,
   the efficient influence function for ψ, computed for each woman from the
   cross-fitted predictions. Its sample mean estimates ψ.
4. **Best linear projection.** OLS regression of the AIPW scores on the
   log-transformed chemicals, unconditionally (one chemical at a time) and
   conditionally (all 16 together), with HC3 robust standard errors. The Wald
   statistics from these models are Figure 1 of the manuscript.
5. **Smoothed functions.** For chemicals whose Wald statistic fell below −1.28
   in either projection model (MiBP, MCOP, MEHP), the AIPW scores are smoothed
   against each log-chemical with a 10-fold cross-validated LOESS (Super
   Learner over a grid of 21 span values). These curves are Figure 2. For the
   smoother only, the sparse tails are excluded (MEHP above its 95th
   percentile, MiBP below its 5th percentile); the rug plots show the data
   used.
6. **Descriptive supplement.** Age-stratified OLS associations between each
   log-chemical and AFC with HC3 standard errors (Table S2), summary
   statistics underlying Table 1 and Table S1, multivariate outlier
   diagnostics, and the timing of chemical sampling relative to the scan.

## Main results

| Quantity | Result |
|---|---|
| Standardized mean difference in AFC, ≥35 vs <35 years (ψ) | −5.1 antral follicles (95% CI −6.2, −3.9) |
| Chemicals with the largest Wald statistics in the linear projections | MiBP, MCOP, MEHP (negative); hair Hg (positive) |
| Shape across concentration (cross-validated LOESS) | MiBP and MCOP curvilinear, with the negative age–AFC difference growing beyond ≈4.5 μg/L and ≈50 μg/L respectively; MEHP approximately linear with no apparent threshold |

Higher MiBP, MCOP, and MEHP concentrations accompanied a more negative age–AFC
difference; higher hair mercury accompanied a less negative one. The manuscript
discusses hair mercury as a likely marker of seafood and omega-3 fatty acid
intake rather than as a mercury effect. The remaining 12 chemicals showed little
evidence of modification.

## Pipeline

Five numbered R scripts plus a Quarto report, orchestrated by the `Makefile`.
Each script reads the previous script's output and uses `here::here()` paths,
so they can be run from any working directory inside the project.

0. **`code/0_data_gen.R`** — data preparation and imputation
   - Reads the raw EARTH extract (`data/afc_mixtures.sas7bdat`).
   - Summarizes the timing of phthalate/BPA and hair-mercury sampling relative
     to the AFC scan (`figures/date_diff_combined.png`,
     `output/date_comparison_summary.rds`, `output/date_timing_counts.rds`).
   - Tabulates missingness and drops variables missing for more than 40% of
     women (BPS, BPF, BP3, TCS, and the OPFR metabolites BDCIPP, DPHP, and
     ipPPP, whose specific-gravity ratio was 75% missing).
   - Imputes the remaining variables with `missForest` and saves out-of-bag
     error and before/after distribution comparisons.
   - Creates non-redundant indicators of which values were imputed.
   - Builds specific-gravity-corrected concentrations (phthalates by the
     phthalate SG ratio; BPA and parabens by the phenol SG ratio). A molar
     ΣDEHP variable is also built but is not used downstream.
   - Writes `data/imputed_EARTH.Rdata` and an ID crosswalk.

1. **`code/1_data_man.R`** — analysis data sets
   - Selects the analysis variables (AFC, scan year and month, age, BMI, race,
     education, smoking, prior IVF, prior IUI, gravidity, the 16 chemicals,
     and the imputation indicators) and saves summary statistics.
   - Writes `data/afc_clean_notrunc.Rdata`, used for all analyses, and
     `data/afc_clean_trunc.Rdata`, in which each chemical is capped at its
     92.5th percentile. The capped version is used only for the descriptive
     comparison in Step 2.

2. **`code/2_chem_analysis.R`** — chemical distributions, outlier diagnostics,
   and age-stratified associations
   - Plots the distributions of AFC, age, and the chemicals with and without
     capping.
   - Runs multivariate outlier diagnostics (Mahalanobis distance, PCA
     distance, local outlier factor, Cook's distance) on AFC, age, and the
     chemicals jointly. These are descriptive; no observations were removed.
   - Fits OLS models of AFC on each log-chemical, unconditionally and
     conditionally on the other chemicals, separately for women <35 and ≥35,
     with HC3 standard errors (`output/association_table_by_age.rds` and
     `.xlsx`; Table S2).

3. **`code/3_IF_scores_gen.R`** — nuisance models and AIPW scores
   - Defines the exposure indicator (age ≥35) and fits the cross-validated
     Super Learner outcome and propensity score models described above. Fits
     are cached in `misc/fit_mu.RDS` and `misc/fit_pi.RDS`; delete them to
     refit.
   - Saves the Super Learner summaries, a permutation-based variable
     importance plot, and the propensity score overlap plot.
   - Computes the AIPW score for each woman and the overall standardized mean
     difference ψ (`output/ate_results.rds`).
   - Writes `data/afc_clean_notrunc_IF.Rdata`, which appends the AIPW scores
     (`dr_scores`) to the analysis data.

4. **`code/4_IF_scores_analysis.R`** — modification of the age–AFC association
   - Pairwise correlations among the log-chemicals.
   - Best linear projections of the AIPW scores on the log-chemicals,
     unconditional and conditional, with HC3 standard errors
     (`output/linear_projection_output.rds`; `figures/teststat_scatter.png`
     is Figure 1).
   - Linear association functions for all 16 chemicals
     (`figures/cate_functions_full_linear.png`).
   - Cross-validated LOESS functions for the chemicals selected by the
     projection models (`figures/cate_functions_paper.png` is Figure 2;
     `output/significant_chemicals_blp.rds` lists the selected chemicals).

**Report.** `code/EARTH_Analysis_Report.qmd` renders every saved table and
figure into `output/EARTH_Analysis_Report.html`, the only generated file
tracked in this repository.

## Requirements

- R (the results were produced with R 4.6.1; R ≥ 4.1 should work).
- Quarto (rendered with 1.10) for the report.
- R packages, loaded with `pacman::p_load()`, which installs any that are
  missing: rio, here, skimr, tidyverse, scales, haven, VIM, naniar,
  missForest, doParallel, foreach, patchwork, lmtest, sandwich, broom,
  reshape2, dbscan, factoextra, MVN, writexl, SuperLearner, vip, ranger,
  xgboost, glmnet, polspline, earth, clubSandwich, mgcv, mgcViz, ggh4x,
  ggrepel, gridExtra, GGally, knitr.

Steps 0 and 3 are the slow steps (random forest imputation and the
cross-validated Super Learner fits). Both use `parallel::detectCores() - 2`
cores. Random seeds are fixed (`set.seed(123)`) in Steps 0, 3, and 4.

### Reproducibility notes

- The manuscript's results come from a single run of the pipeline on
  2025-12-01. The random forest imputation in Step 0 is parallelized across
  variables with `doParallel`, and `set.seed()` does not reach the worker
  processes, so re-running Step 0 returns a different imputed data set:
  observed values are unchanged, but every imputed value differs. Re-running
  Steps 0 through 3 therefore reproduces the manuscript's estimates only up
  to imputation noise. The 2025-12-01 realization is preserved locally in
  `data/afc_clean_notrunc_IF.Rdata`, and `data/afc_clean_notrunc.Rdata` and
  `data/afc_clean_trunc.Rdata` hold the same realization, so Steps 2 and 4
  run on the data the paper used. For reproducible parallel imputation in
  future runs, register `doRNG::registerDoRNG(123)` after
  `registerDoParallel()`.
- The 2025-12-01 nuisance models also included an "ever smoker" indicator
  that is a deterministic function of smoking status (never, former,
  current). It was dropped from the analysis data on 2025-12-05 as
  redundant, after the main results had been generated. The cached Super
  Learner fits in `misc/` were trained with it and fail when applied to the
  current covariate matrix; delete them before re-running Step 3.
- Steps 2 and 4 are deterministic given their inputs. On the preserved data
  they reproduce Table S2, Figure 1, and Figure 2 exactly (checked
  2026-10-09).

## Usage

```bash
# Full pipeline: Steps 0-4, then the report
make all

# Individual steps (each also runs any earlier step whose outputs are stale)
make step0   # data preparation and imputation (slow)
make step1   # analysis data sets
make step2   # chemical distributions, outliers, age-stratified associations
make step3   # Super Learner fits and AIPW scores (slow)
make step4   # linear projections and smoothed association functions
make report  # render the HTML report only

make help    # list targets
make clean-all   # remove ALL generated data, figures, and outputs (destructive)
```

`make` decides what to rerun from file timestamps. To rerun a single script
without touching upstream steps, call it directly:

```bash
cd code && Rscript 2_chem_analysis.R
quarto render code/EARTH_Analysis_Report.qmd --output-dir ../output
```

## Repository structure

```
EARTH_Analysis/
├── code/
│   ├── 0_data_gen.R               # data preparation and imputation
│   ├── 1_data_man.R               # analysis data sets
│   ├── 2_chem_analysis.R          # distributions, outliers, Table S2
│   ├── 3_IF_scores_gen.R          # Super Learner fits, AIPW scores, psi
│   ├── 4_IF_scores_analysis.R     # linear projections, Figures 1 and 2
│   └── EARTH_Analysis_Report.qmd  # Quarto source of the report
├── output/
│   └── EARTH_Analysis_Report.html # rendered report (tracked)
├── Makefile
├── EARTH_Analysis.Rproj
└── README.md

Not tracked (generated locally or restricted):
├── data/        # raw EARTH extract and derived data sets
├── figures/     # generated figures
├── output/*     # generated tables (.rds, .xlsx), apart from the report
├── misc/        # cached Super Learner fits, outlier CSVs
├── manuscript/  # manuscript drafts
└── sandbox/     # exploratory code
```

## Data availability

The EARTH Study data are not publicly available because of privacy and
confidentiality restrictions and the terms of the data use agreement. Requests
to access the data should be directed to Russ Hauser
(rhauser@hsph.harvard.edu), as stated in the manuscript. The pipeline expects
the raw extract at `data/afc_mixtures.sas7bdat`.

## Funding

National Institute of Environmental Health Sciences of the National Institutes
of Health, award numbers P30ES019776 and ES009718. The content is solely the
responsibility of the authors and does not necessarily represent the official
views of the National Institutes of Health.

## Citation

Until the journal version appears, please cite the preprint:

```bibtex
@article{naimi2025environmental,
  title   = {Environmental Chemicals as Modifiers of the Association between
             Chronological Age and Ovarian Reserve},
  author  = {Naimi, Ashley I and Kennedy, Edward H and Yu, Ya-Hui and
             Hauser, Russ and M{\'i}nguez-Alarc{\'o}n, Lidia and Gaskins, Audrey},
  journal = {medRxiv},
  year    = {2025},
  doi     = {10.64898/2025.12.09.25341902}
}
```

Methodological reference for the DR-learner: Kennedy EH. Towards optimal
doubly robust estimation of heterogeneous causal effects. *Electronic Journal
of Statistics*. 2023;17(2):3008–3049.
[doi:10.1214/23-EJS2157](https://doi.org/10.1214/23-EJS2157).

## Contact

Ashley I. Naimi, Department of Epidemiology, Rollins School of Public Health,
Emory University. ashley.naimi@emory.edu

## License

MIT License

Copyright (c) 2025 Ashley I. Naimi

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
