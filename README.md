# Dietary Stearic-to-Palmitic Acid Ratio and Systemic Inflammation]{Dietary Stearic-to-Palmitic Acid Ratio and Systemic Inflammation: A Dual-Track Analysis Integrating Observational Associations with a Machine Learning-Derived Clinical Prediction Nomogram

Welcome to the official repository for the **SCI-Risk-Nomogram** project. This repository contains the complete replication code, data schemas, and analytical pipelines for our study investigating the observational association between micro-level dietary saturated fatty acid imbalances (stearic acid 18:0 to palmitic acid 16:0 ratio) and systemic chronic inflammation (SCI), as well as the development of a parsimonious clinical screening tool.

---

## 📌 Project Architecture

This study implements a novel **Dual-Track Analytical Framework** utilizing population-based data from the National Health and Nutrition Examination Survey (NHANES) 2015–2018:

```text
                  ┌─────────────────────────────────────────┐
                  │       NHANES 2015–2018 Cohort          │
                  │             (N = 8,382)                 │
                  └────────────────────┬────────────────────┘
                                       │
                    ┌──────────────────┴──────────────────┐
                    ▼                                     ▼
     【Track I: Etiological Track】        【Track II: Clinical Translation】
        Full Cleaned Cohort                     Fasting Serological Sub-cohort
            (N = 8,382)                                  (N = 3,887)
                    │                                     │
    ┌───────────────┴───────────────┐                     ├──────────────────────┐
    ▼                               ▼                     ▼                      ▼
Observational Models         Mechanism Discovery    70% Training Set       30% Testing Set
- Energy-Adj (Residual)      - PCA (Fatty Acids)       (n = 2,721)            (n = 1,166)
- RCS Linear Fit             - ACME Mediation             │                      │
- PSM (1:1 Matching)           (TyG & SII paths)          ▼                      ▼
- Sensitivity Analyses                                 LASSO Selection        Blinded Validation
                                                      4-Var Nomogram Fit     - C-index (0.769)
                                                      (Gender, BMI, TyG,     - Brier Score (0.177)
                                                           and SII)          - DCA Benefit Curve
```

- **Track I: Robust Association & Mechanistic Exploration ($N = 8,382$)**
  Establishes the independent linear dose-response relationship between the energy-adjusted dietary 18:0/16:0 ratio and elevated high-sensitivity C-reactive protein (hs-CRP $> 3$ mg/L). Features rigorous confounding control via Propensity Score Matching (PSM), exploratory Principal Component Analysis (PCA), and association-based mediation modeling (Average Causal Mediation Effects - ACME).
- **Track II: Machine Learning & Clinical Risk Prediction ($N = 3,887$)**
  Translates upstream nutritional etiology into downstream clinical utility within a fasting sub-cohort. Implements an outcome-stratified 70/30 physical split, Least Absolute Shrinkage and Selection Operator (LASSO) feature contraction, multivariable nomogram deployment, and fully independent blinded validation alongside Decision Curve Analysis (DCA).

---

## 📂 Repository Layout

```text
├── data/
│   ├── Stage_一_2_Cleaned_Cohort.rds     # Cleaned full analytical cohort database (N=8,382)
│   └── Stage_五_Locked_Test_Set.rds      # Strictly isolated 30% blind test subset (N=1,166)
├── scripts/
│   └── Master_Analytical_Pipeline.R      # End-to-head consolidated compilation R script
├── results/
│   ├── Stage_二_3_RCS_Curve.pdf          # Multivariable restricted cubic spline plot
│   ├── Stage_三_5_Subgroup_Forest.pdf    # Clinical stratified stratification forest plot
│   ├── Stage_五_1_LASSO_Path.pdf         # Cross-validation coefficient penalty path
│   ├── Stage_五_2_Training_Nomogram.pdf  # Final dynamic visual predictive nomogram
│   ├── Stage_六_1_Testing_Calibration.pdf# Tester calibration continuous curve
│   └── Stage_六_2_Testing_DCA.pdf        # Clinical decision curve analysis benefit output
└── README.md                             # Repository introductory roadmap document
```

---

## 💻 Environment & Dependencies

All data processing, multivariable modeling, and visualization were executed under **R Environment (v4.2.2)**. The package infrastructure is strictly anchored to the following versions for full mathematical reproducibility:

```R
# Core library setup and version pinning
install.packages("pacman")
library(pacman)
p_load(
  nhanesA,    # v0.7.1  - Programmatic data extraction from CDC servers
  MatchIt,    # v4.5.5  - Nearest-neighbor propensity score matching
  mediation,  # v4.5.0  - Nonparametric bootstrap mediation modeling
  psych,      # v2.3.3  - Principal component analysis with varimax rotation
  glmnet,     # v4.1-7  - High-dimensional LASSO regularization
  rms,        # v6.7-1  - Predictive nomogram fitting and internal validation
  dcurves,    # v0.4.0  - Clinical Net Benefit Decision Curve Analysis
  broom,      # v1.0.4  - Clean model parameter tidying
  dplyr       # v1.1.1  - Data manipulation base engine
)
```

---

## 🚀 Step-by-Step Analytical Pipeline

The file `scripts/Master_Analytical_Pipeline.R` contains the complete end-to-end execution path. Below is a breakdown of the core methodological steps:

### 1. Energy Adjustment via the Residual Method
To eliminate confounding from total caloric volume, absolute dietary intakes of stearic acid (18:0) and palmitic acid (16:0) are individually regressed on total energy intake (`Energy_Kcal`). The math follows:
$$\text{Intake}_{\text{adjusted}} = \text{Residual}_{\text{lm}} + \text{Mean}(\text{Intake}_{\text{population}})$$
The final exposure variable is synthesized as: $\text{Ratio}_{\text{adj}} = \frac{18:0_{\text{adj}}}{16:0_{\text{adj}} + 0.0001}$.

### 2. Track I Association & Sensitivity Backstops
- **RCS Curve:** Fits a 4-knot restricted cubic spline to test non-linear deviation ($P_{\text{non-linear}} = 0.7971$, confirming strict linearity).
- **Confounding Matching (PSM):** Implements a 1:1 nearest-neighbor matching algorithm under a tight caliper (0.05). Achieves superb balance with all 13 baseline covariates reaching a Standardized Mean Difference (SMD $< 0.05$, maximum absolute SMD $= 0.0458$).
- **Extreme Multicollinearity Defense:** Verifies that simultaneously forcing absolute SFA intakes into the model triggers severe variance inflation (VIF $> 7$), statistically vindicating our ratio-based model.
- **Multidimensional Sensitivity Checks:** Sequentially isolates the independent association by fitting cohorts strictly excluding history of malignancy (OR $= 1.20, P = 0.015$), diabetes (OR $= 1.178, P = 0.0458$), and hypertension (OR $= 1.230, P = 0.0331$), as well as testing an elevated diagnostic threshold (hs-CRP $> 5$ mg/L: OR $= 1.258, P = 0.0154$).

### 3. Track I Mechanism Deconvolution
- **Dietary Deconvolution (PCA):** Extracts 3 components with eigenvalues $> 1$, explaining $63.72\%$ of total fatty acid variance. RC1 isolates the overarching "long-chain SFA pattern" ($34.86\%$ variance explained), which shows zero independent link to inflammation ($P = 0.7446$), highlighting that the micro-level 18:0/16:0 ratio is uniquely associated with SCI.
- **Biomarker Mediation (ACME):** Deploys 500-resample nonparametric bootstrap mediation models. Confirms a highly dominant direct pathway (ADE $= 0.015, P = 0.016$) completely unmediated by metabolic insulin resistance (TyG Index ACME $= -0.0003, 95\%$ CI: $-0.0015$ to $0.0013$) and minimally impacted by the systemic immune microenvironment (SII Index ACME $= 0.0011, 95\%$ CI: $-0.0000$ to $0.0022$).

### 4. Track II Machine Learning & Prediction Translation
- **Overfitting Prevention:** Verifies an adequate Events Per Variable ratio within the training cohort ($n=2,721, 841$ events, EPV $\approx 52$) to guard against model over-parameterization.
- **LASSO Contraction:** Compresses a 16-variable candidate pool down to 4 stable host-milieu signatures under the 1-SE parsimony criterion: Gender, BMI, TyG Index, and SII Index. The upstream dietary exposure is objectively compressed out due to high intra-individual variance.
- **Blinded Test Set Validation:** Evaluates the model directly on the physically isolated testing subset ($n=1,166, 360$ events), capturing stable discrimination (C-index $= 0.769$) and superb calibration (Brier Score $= 0.177$).
- **Value-Added Verification:** Confirms via formal DeLong and Net Reclassification tests that forcing dietary recalls into the clinical screening tool yields zero statistical gain ($\Delta$AUC $P = 0.570$; NRI $= 0.042, P = 0.500$; IDI $= 0.002, P = 0.106$).
- **DCA Decision Curves:** Clinically justifies the nomogram by showing strong, stable net benefit gains across a wide continuum of threshold probabilities.

---

## 📈 Summary of Main Scientific Findings

| Metric / Parameter | Experimental Statistical Output | Scientific Contextual Verdict |
| :--- | :--- | :--- |
| **Track I: Q4 vs Q1 OR** | **1.20** ($95\%$ CI: $1.04 - 1.40, P = 0.013$) | Robust, linear independent risk upregulation. |
| **Post-Match Covariates**| Max Absolute **SMD = 0.0458** | Perfect balance across extreme strata. |
| **TyG Path ACME** | **-0.0003** ($95\%$ CI: $-0.0015$ to $0.0013, P = 0.736$) | Zero mediation through metabolic dysregulation. |
| **SII Path ACME** | **0.0011** ($95\%$ CI: $-0.0000$ to $0.0022, P = 0.080$) | Suggestive but non-significant immune pathway. |
| **Validation C-index** | **0.769** in physically isolated testing cohort | Superb discriminative capacity and generalizability. |
| **Validation Brier** | **0.177** in physically isolated testing cohort | High accuracy, zero predictive calibration drift. |
| **Dietary Addition Gain**| $\Delta$AUC $P=0.570$, NRI $P=0.500$, IDI $P=0.106$| Vindicates clinical screening variable parsimony. |

---

## ⚖️ License & Terms of Use
This repository is open-sourced under the MIT License. Data vectors are synthesized directly from public CDC/NHANES records. The R codebase is entirely original and free for scholastic adaptation and clinical translation purposes with proper attribution.

## ✉️ Contact & Collaboration
For scientific inquiries regarding data architectures, multi-omics translation, or clinical prehabilitation screening protocols, please contact:
- **Corresponding Author:** Dr. Hongbo He (`hhb89008684@163.com`)
- **Lead Analyst:** Bing Tang, Department of Integrated Traditional Chinese and Western Medicine Surgery, West China Hospital, Sichuan University.
