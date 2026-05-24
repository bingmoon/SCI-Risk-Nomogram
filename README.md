# Development and Internal Validation of a Multivariable Nomogram for Systemic Chronic Inflammation: A Population-Based Machine Learning Study

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.3.0-blue)](https://www.r-project.org/)
[![Data](https://img.shields.io/badge/Data-NHANES%202017--2018-green)](https://wwwn.cdc.gov/nchs/nhanes/)

## Overview
This repository contains the complete analytical pipeline, feature engineering datasets, and master R script for the manuscript: **"A Cost-Effective Multivariable Nomogram for Predicting Systemic Chronic Inflammation: Implications for Outpatient Screening and Surgical Prehabilitation"**. 

The study leverages interpretable machine learning (SHAP) and restricted cubic splines (RCS) to construct a non-invasive, high-performance clinical nomogram integrating the dietary stearic-to-palmitic acid (18:0/16:0) ratio with established immune-metabolic baselines (TyG and SII indices).

## Repository Structure
- `analysis_master_script.R`: The core executable script. It encompasses the entire workflow from NHANES data fetching, table generation, machine learning (Random Forest & SHAP), RCS dose-response analysis, nomogram construction, internal validation (Bootstrap), to Decision Curve Analysis (DCA).
- `data/`: Contains the processed datasets required to replicate the statistical models.
  - `Stage1_NHANES_Cleaned_Raw.csv`: Cleaned foundational cohort data.
  - `Stage2_Features_Engineered.csv`: Engineered dataset including calculated composite indices (TyG, SII) and lipid ratios.

## Computational Reproducibility
To ensure absolute methodological transparency, a fixed random seed (`set.seed(2026)`) is uniformly applied across all stochastic algorithms (including Random Forest tree generation, SHAP Monte Carlo simulations, and bootstrap resampling). Reviewers and researchers can execute the master script to reproduce identical figures, C-indices, and P-values reported in the manuscript.

### Prerequisites
The analysis was performed using **R software (version 4.3.x or higher)**. The script automatically checks and installs the necessary packages via `pacman`. Key dependencies include:
- `nhanesA` (Data retrieval)
- `tableone` (Baseline characteristics)
- `randomForest`, `fastshap`, `shapviz` (Machine learning & Interpretability)
- `rms` (RCS integration & Nomogram construction)
- `dcurves` (Decision Curve Analysis)
- `ggplot2`, `dplyr`, `tidyr` (Data manipulation & Visualization)

### How to Run
1. Clone this repository to your local machine:
   ```bash
   git clone [https://github.com/bingmoon/SCI-Risk-Nomogram.git](https://github.com/bingmoon/SCI-Risk-Nomogram.git)
