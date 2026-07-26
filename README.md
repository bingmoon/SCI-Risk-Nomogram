# Dietary Stearic-to-Palmitic Acid Ratio and Systemic Inflammation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21388142.svg)](https://doi.org/10.5281/zenodo.21388142)

This repository contains the complete analytical code for the manuscript:

**Dietary stearic-to-palmitic acid ratio and systemic inflammation: a cross-sectional dual-track analysis with independent external validation**

*Submitted to SAGE Open Medicine* (previously under review at other journals)

---

## 📄 About the Study

This study employs a dual‑track framework to investigate the association between the dietary 18:0/16:0 fatty acid ratio and systemic chronic inflammation, and to develop a pragmatic clinical prediction model.

- **Track I** – Observational association and mediation analysis (NHANES 2015–2018, *N* = 8,382)
- **Track II** – Machine learning‑based prediction model with independent external validation (CHARLS 2015, *N* = 12,613)

> **Note on the exposure**: The dietary 18:0/16:0 ratio was energy‑adjusted by first calculating the raw ratio of stearic to palmitic acid, then applying the residual method to the ratio itself. This is the corrected approach used in the final manuscript, and differs from earlier versions that adjusted the individual fatty acids before calculating the ratio. The code in this repository reflects the final, peer‑reviewed method.

> **Note on data availability**: Blood biomarker data from more recent CHARLS waves (2018, 2020) have not yet been publicly released, therefore CHARLS 2015 remains the most recent wave with complete biochemical data necessary for model validation.

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `analysis_master_script.R` | Complete analytical pipeline – data extraction, statistical modeling, figure generation |
| `README.md` | This file |
| `LICENSE` | MIT License |

All results (tables, figures, intermediate data) are generated automatically when you run the script – no additional files are required.

---

## 🚀 How to Run

### 1. Prerequisites

- **R** version ≥ 4.2.2
- **RStudio** (recommended)
- An internet connection (NHANES data are downloaded automatically via the `nhanesA` package)

### 2. Install Required Packages

The script will automatically install any missing packages using `pacman`. You can also install them manually:

```r
install.packages(c(
  "nhanesA", "dplyr", "tidyr", "tableone", "rms", "ggplot2", "broom",
  "MatchIt", "survey", "psych", "mediation", "glmnet", "caret",
  "pROC", "Hmisc", "ResourceSelection", "car", "haven", "dcurves",
  "survival", "cobalt", "splines"
))

3. Set Up Local Paths

You must edit the following paths inside analysis_master_script.R before running:

Working directory – change setwd() to your local project folder
CHARLS 2015 data – update the path in the read_dta() calls to point to your downloaded CHARLS .dta files
⚠️ The script contains absolute paths (e.g., /Users/bing/...). Replace them with your own local paths.
4. Run the Script

Open analysis_master_script.R in RStudio, select all code (Ctrl/Cmd + A), and click Run. The script will execute all analyses in the correct order and save the outputs.

📊 Outputs

After the script completes, you will find two new folders in your working directory:

data/ – Intermediate R objects (.rds files)
results/ –

Publication‑quality figures (PDF format)
Tables (CSV format)
No pre‑computed results are included in this repository – everything is generated on‑the‑fly for full reproducibility.

🔒 Reproducibility

A fixed random seed (2026) is used for all stochastic processes (data splitting, bootstrap, matching).
All package versions are managed by pacman.
LASSO‑selected dummy variables are explicitly mapped back to original factor names to ensure transparent feature selection.
The corrected energy‑adjustment method (raw ratio → residual adjustment) is applied.
📦 Archived Version

The version of the code that accompanies the final manuscript is permanently archived on Zenodo with the DOI:

10.5281/zenodo.21388142

The Zenodo archive contains a snapshot of this repository as of the date of submission, ensuring complete reproducibility of the published results.

📧 Contact

For questions about the code or the manuscript, please contact:

Hongbo He
Division of Surgery, Institute of Integrated Traditional Chinese and Western Medicine,
West China Hospital, Sichuan University
Chengdu, China
Email: hhb89008684@163.com

📜 License

This project is licensed under the MIT License – see the LICENSE file for details.
