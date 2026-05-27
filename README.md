# SCI-Risk-Nomogram: Machine Learning for Systemic Chronic Inflammation

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-4.2.2-green.svg)](https://www.r-project.org/)

## Overview
This repository contains the complete, reproducible analytical pipeline (R code and feature engineering scripts) for the manuscript: **"Development and Internal Validation of a Multivariable Nomogram for Systemic Chronic Inflammation: A Population-Based Machine Learning Study."** Our study integrates machine learning (Random Forest & SHAP) with traditional epidemiological modeling (Restricted Cubic Splines & DCA) to evaluate the threshold-dependent protective effect of the dietary stearic-to-palmitic acid (18:0/16:0) ratio on systemic chronic inflammation (hs-CRP > 3 mg/L).

## Directory Structure & Working Environment
To ensure absolute methodological transparency and local computational reproducibility, all R scripts are structurally designed to be executed within a fixed absolute root directory.

**Primary Working Directory:** `/Users/bing/AA`

Please ensure that your local R Studio working directory is set to this absolute path before executing the master script. The script will automatically generate the following sub-directories upon initialization:
* `/results/data/`: Outputs of all cleaned CSVs, engineered features, and statistical tables (e.g., Table 1).
* `/results/plots/`: High-resolution PDFs for all clinical charts (RCS curves, Nomograms, Calibration Curves, DCA).

## Analytical Pipeline (Master Script)
The entire analysis is consolidated into a single, seamless master script. The pipeline is structured into 7 sequential stages:
1. **Stage 1:** Data retrieval from the NHANES (2017-2018) database, extraction of complex covariates (handling Skip Patterns), and baseline Table 1 generation.
2. **Stage 2:** Feature engineering (TyG and SII indices) and initial Random Forest modeling.
3. **Stage 3:** Black-box interpretability analysis using SHapley Additive exPlanations (SHAP) to ascertain feature directionality.
4. **Stage 4:** Non-linear dose-response assessment utilizing Restricted Cubic Splines (RCS).
5. **Stage 5:** Construction and bootstrap-based internal validation (Calibration & C-index) of the clinical nomogram.
6. **Stage 6:** Clinical net benefit evaluation via Decision Curve Analysis (DCA).
7. **Stage 7:** Complete-case sensitivity analyses (Energy adjustment, stringent hs-CRP thresholding, and massive confounder adjustment).

## Dependencies
* Core modeling: `randomForest`, `fastshap`, `shapviz`, `rms`, `dcurves`
* Data processing & Visualization: `nhanesA`, `dplyr`, `tidyr`, `ggplot2`, `tableone`

## Contact & Correspondence
For any inquiries regarding the code, data engineering logic, or statistical methodology, please contact the corresponding author via official institutional academic email channels. The use of personal email addresses for academic correspondence regarding this repository is strictly prohibited.
