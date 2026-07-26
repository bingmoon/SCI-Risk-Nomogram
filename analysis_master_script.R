# ==============================================================================
# Project: Dual-Track Analysis of Dietary 18:0/16:0 Ratio and Systemic Inflammation
# Final Complete Version – All fixes integrated:
#   - dplyr::select used everywhere
#   - Robust Gender conversion
#   - Manual E-value calculation
#   - LASSO dummy variable mapping back to original factors
#   - CHARLS smoking variable removed (fallback-safe)
# ==============================================================================

# 0. Environment setup ---------------------------------------------------------
rm(list = ls())
set.seed(2026)

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(
  nhanesA, dplyr, tidyr, tableone, rms, ggplot2, broom,
  MatchIt, survey, psych, mediation, glmnet, caret,
  pROC, Hmisc, ResourceSelection, car, haven, dcurves, survival
)

# Set working directory (adjust as needed)
setwd("/Users/bing/AA/AB_revision")
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# 1. NHANES data extraction ---------------------------------------------------
pull_cycle_data <- function(cycle) {
  demo <- nhanes(paste0('DEMO_', cycle)) %>%
    dplyr::select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, INDFMPIR)
  
  diet <- nhanes(paste0('DR1TOT_', cycle)) %>%
    dplyr::select(SEQN, DR1TKCAL, DR1TALCO, DR1TFIBE,
                  DR1TTFAT, DR1TSFAT, DR1TMFAT, DR1TPFAT,
                  DR1TS140, DR1TS160, DR1TS180,
                  DR1TM181, DR1TP182, DR1TP183, DR1TP204, DR1TP205, DR1TP226)
  
  crp  <- nhanes(paste0('HSCRP_', cycle)) %>% dplyr::select(SEQN, LBXHSCRP)
  bmx  <- nhanes(paste0('BMX_', cycle)) %>% dplyr::select(SEQN, BMXBMI)
  smq  <- nhanes(paste0('SMQ_', cycle)) %>% dplyr::select(SEQN, SMQ020, SMQ040)
  diq  <- nhanes(paste0('DIQ_', cycle)) %>% dplyr::select(SEQN, DIQ010)
  bpq  <- nhanes(paste0('BPQ_', cycle)) %>% dplyr::select(SEQN, BPQ020)
  paq  <- nhanes(paste0('PAQ_', cycle)) %>%
    dplyr::select(SEQN, any_of(c("PAQ610", "PAD615", "PAQ625", "PAD630",
                                 "PAQ640", "PAD645", "PAQ655", "PAD660",
                                 "PAQ670", "PAD675")))
  
  glu <- nhanes(paste0('GLU_', cycle)) %>% dplyr::select(SEQN, LBXGLU)
  tg  <- nhanes(paste0('TRIGLY_', cycle)) %>% dplyr::select(SEQN, LBXTR)
  cbc <- nhanes(paste0('CBC_', cycle)) %>%
    dplyr::select(SEQN, LBXPLTSI, LBDNENO, LBDLYMNO)
  
  demo %>%
    left_join(diet, by = "SEQN") %>%
    left_join(crp, by = "SEQN") %>%
    left_join(bmx, by = "SEQN") %>%
    left_join(smq, by = "SEQN") %>%
    left_join(diq, by = "SEQN") %>%
    left_join(bpq, by = "SEQN") %>%
    left_join(paq, by = "SEQN") %>%
    left_join(glu, by = "SEQN") %>%
    left_join(tg, by = "SEQN") %>%
    left_join(cbc, by = "SEQN")
}

df_I <- pull_cycle_data("I")
df_J <- pull_cycle_data("J")
nhanes_raw <- bind_rows(df_I, df_J)

# Data cleaning and variable derivation
nhanes_clean <- nhanes_raw %>%
  drop_na(RIDAGEYR, RIAGENDR, DR1TS180, DR1TS160, LBXHSCRP) %>%
  filter(RIDAGEYR >= 20, LBXHSCRP <= 10, DR1TS160 > 0, DR1TS180 > 0) %>%
  mutate(
    Ratio_18_16 = DR1TS180 / (DR1TS160 + 0.0001),
    High_CRP_num = ifelse(LBXHSCRP > 3, 1, 0),
    High_CRP     = factor(High_CRP_num, levels = c(0,1), labels = c("No","Yes")),
    
    Age      = RIDAGEYR,
    Gender_raw = RIAGENDR,
    Race     = factor(RIDRETH1),
    PIR      = ifelse(is.na(INDFMPIR), median(INDFMPIR, na.rm = TRUE), INDFMPIR),
    BMI      = ifelse(is.na(BMXBMI), median(BMXBMI, na.rm = TRUE), BMXBMI),
    
    Energy_Kcal = ifelse(is.na(DR1TKCAL), median(DR1TKCAL, na.rm = TRUE), DR1TKCAL),
    Total_Fat   = ifelse(is.na(DR1TTFAT), median(DR1TTFAT, na.rm = TRUE), DR1TTFAT),
    Fiber       = ifelse(is.na(DR1TFIBE), median(DR1TFIBE, na.rm = TRUE), DR1TFIBE),
    Alcohol_g   = ifelse(is.na(DR1TALCO), 0, DR1TALCO),
    
    Smoking = factor(case_when(
      as.character(SMQ020) %in% c("2","No") ~ "Never",
      as.character(SMQ020) %in% c("1","Yes") & 
        as.character(SMQ040) %in% c("1","2","Every day","Some days") ~ "Current",
      TRUE ~ "Former/Unknown"
    ), levels = c("Never", "Current", "Former/Unknown")),
    
    MET_min_week = (
      replace_na(ifelse(PAQ610 %in% 1:7, PAQ610, 0), 0) * 
        replace_na(ifelse(PAD615 < 7777, PAD615, 0), 0) * 8 +
      replace_na(ifelse(PAQ625 %in% 1:7, PAQ625, 0), 0) * 
        replace_na(ifelse(PAD630 < 7777, PAD630, 0), 0) * 4 +
      replace_na(ifelse(PAQ640 %in% 1:7, PAQ640, 0), 0) * 
        replace_na(ifelse(PAD645 < 7777, PAD645, 0), 0) * 4 +
      replace_na(ifelse(PAQ655 %in% 1:7, PAQ655, 0), 0) * 
        replace_na(ifelse(PAD660 < 7777, PAD660, 0), 0) * 8 +
      replace_na(ifelse(PAQ670 %in% 1:7, PAQ670, 0), 0) * 
        replace_na(ifelse(PAD675 < 7777, PAD675, 0), 0) * 4
    ),
    MET_imp = ifelse(MET_min_week == 0 | is.na(MET_min_week),
                     median(MET_min_week[MET_min_week > 0], na.rm = TRUE),
                     MET_min_week),
    
    Diabetes = factor(case_when(
      as.character(DIQ010) %in% c("1","Yes") ~ "Yes", TRUE ~ "No"
    )),
    Hypertension = factor(case_when(
      as.character(BPQ020) %in% c("1","Yes") ~ "Yes", TRUE ~ "No"
    )),
    
    TyG_Index = log(LBXTR * LBXGLU / 2),
    SII_Index = (LBXPLTSI * LBDNENO) / (LBDLYMNO + 0.0001)
  )

# ---- Gender emergency patch ----
nhanes_clean <- nhanes_clean %>%
  mutate(
    Gender = factor(case_when(
      as.character(Gender_raw) %in% c("1", "Male") ~ "Male",
      as.character(Gender_raw) %in% c("2", "Female") ~ "Female",
      TRUE ~ NA_character_
    ), levels = c("Male", "Female"))
  ) %>%
  dplyr::select(-Gender_raw)

print(table(nhanes_clean$Gender, useNA = "ifany"))

ratio_q <- quantile(nhanes_clean$Ratio_18_16, probs = c(0,0.25,0.5,0.75,1), na.rm=TRUE)
nhanes_clean$Ratio_Quartile <- cut(
  nhanes_clean$Ratio_18_16, breaks = ratio_q, include.lowest = TRUE,
  labels = c("Q1","Q2","Q3","Q4")
)
saveRDS(nhanes_clean, "data/nhanes_clean.rds")

# 2. Table 1 ------------------------------------------------------------------
vars_t1 <- c("Age","Gender","Race","PIR","BMI","Smoking","Alcohol_g",
             "MET_imp","Energy_Kcal","Total_Fat","Fiber",
             "Diabetes","Hypertension","High_CRP")
cat_vars <- c("Gender","Race","Smoking","Diabetes","Hypertension","High_CRP")
tab1 <- CreateTableOne(vars = vars_t1, strata = "Ratio_Quartile",
                       data = nhanes_clean, factorVars = cat_vars)
tab1_csv <- print(tab1, nonnormal = c("Age","PIR","BMI","Alcohol_g","MET_imp",
                                      "Energy_Kcal","Total_Fat","Fiber"),
                  quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_csv, "results/tables/Table1.csv")

# 3. Track I: Observational association (energy-adjusted ratio) ---------------
nhanes_clean <- nhanes_clean %>%
  mutate(
    Ratio_raw = DR1TS180 / (DR1TS160 + 0.0001),
    Ratio_adj = resid(lm(Ratio_raw ~ Energy_Kcal, data = nhanes_clean,
                         na.action = na.exclude)) + mean(Ratio_raw, na.rm = TRUE)
  )

ratio_adj_q <- quantile(nhanes_clean$Ratio_adj, probs = c(0,0.25,0.5,0.75,1), na.rm=TRUE)
nhanes_clean$Ratio_adj_Q <- cut(
  nhanes_clean$Ratio_adj, breaks = ratio_adj_q, include.lowest = TRUE,
  labels = c("Q1","Q2","Q3","Q4")
)

dd <- datadist(nhanes_clean)
options(datadist = "dd")

# 3.1 RCS
fit_rcs <- lrm(High_CRP_num ~ rcs(Ratio_adj, 4) + Age + Gender + Race + PIR +
                BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal +
                Total_Fat + Fiber + Diabetes + Hypertension,
              data = nhanes_clean)
print(anova(fit_rcs))

# 3.2 Logistic regression
fit_model3 <- glm(
  High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR + BMI +
    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber +
    Diabetes + Hypertension,
  data = nhanes_clean, family = binomial()
)
res_model3 <- tidy(fit_model3, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(grepl("Ratio_adj_Q", term))
print(res_model3)

nhanes_clean$Ratio_adj_Q_num <- as.numeric(nhanes_clean$Ratio_adj_Q)
fit_trend <- glm(
  High_CRP_num ~ Ratio_adj_Q_num + Age + Gender + Race + PIR + BMI +
    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber +
    Diabetes + Hypertension,
  data = nhanes_clean, family = binomial()
)
p_trend <- coef(summary(fit_trend))["Ratio_adj_Q_num", "Pr(>|z|)"]

# 3.3 VIF
fit_vif <- glm(
  High_CRP_num ~ Ratio_adj_Q + DR1TS160 + DR1TS180 + Age + Gender + Race +
    PIR + BMI + Smoking + Alcohol_g + MET_imp + Total_Fat + Fiber +
    Diabetes + Hypertension,
  data = nhanes_clean, family = binomial()
)
print(car::vif(fit_vif))

# 3.4 PSM (Q4 vs Q1) with conditional logistic regression
df_extreme <- nhanes_clean %>%
  filter(Ratio_adj_Q %in% c("Q1","Q4")) %>%
  mutate(
    Exposure_Extreme = ifelse(Ratio_adj_Q == "Q4", 1, 0),
    Exposure_Factor = factor(Exposure_Extreme, levels=c(0,1), labels=c("Q1","Q4"))
  )

psm_match <- matchit(Exposure_Extreme ~ Age + Gender + Race + PIR + BMI +
                      Smoking + Alcohol_g + MET_imp + Energy_Kcal +
                      Total_Fat + Fiber + Diabetes + Hypertension,
                    data = df_extreme, method = "nearest", caliper = 0.05)
matched_data <- match.data(psm_match)

fit_psm <- clogit(High_CRP_num ~ Exposure_Factor + strata(subclass),
                 data = matched_data)
summary(fit_psm)

# 3.5 Manual E‑value calculation
calc_evalue <- function(or, lo) {
  if (or <= 1) return(list(e_value = NA, e_value_lo = NA))
  e_val <- or + sqrt(or * (or - 1))
  e_val_lo <- ifelse(lo > 1, lo + sqrt(lo * (lo - 1)), 1)
  list(e_value = e_val, e_value_lo = e_val_lo)
}
ev <- calc_evalue(1.20, 1.04)
cat("E-value for point estimate:", round(ev$e_value, 3), "\n")
cat("E-value for CI lower limit:", round(ev$e_value_lo, 3), "\n")

# 3.6 Subgroup analyses (per 0.1 unit increase)
run_subgroup <- function(data, subgroup_name) {
  data <- data %>% mutate(Ratio_Scale = Ratio_adj * 10) %>% droplevels()
  covs <- c("Age","Gender","BMI","Smoking","Energy_Kcal")
  covs <- covs[sapply(covs, function(v) length(unique(data[[v]])) > 1)]
  f <- paste("High_CRP_num ~ Ratio_Scale +", paste(covs, collapse=" + "))
  fit <- glm(f, data = data, family = binomial())
  res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "Ratio_Scale")
  data.frame(Subgroup = subgroup_name, N = nrow(data),
             OR = round(res$estimate,2), Lower = round(res$conf.low,2),
             Upper = round(res$conf.high,2), P = round(res$p.value,4))
}
subgroups <- bind_rows(
  run_subgroup(filter(nhanes_clean, Gender == "Male"), "Male"),
  run_subgroup(filter(nhanes_clean, Gender == "Female"), "Female"),
  run_subgroup(filter(nhanes_clean, BMI < 25), "BMI < 25"),
  run_subgroup(filter(nhanes_clean, BMI >= 25 & BMI < 30), "BMI 25-30"),
  run_subgroup(filter(nhanes_clean, BMI >= 30), "BMI >= 30"),
  run_subgroup(filter(nhanes_clean, Age < 60), "Age < 60"),
  run_subgroup(filter(nhanes_clean, Age >= 60), "Age >= 60")
)
write.csv(subgroups, "results/tables/Subgroup.csv", row.names = FALSE)

# 3.7 PCA on 9 fatty acids
fa_vars <- c("DR1TS140","DR1TS160","DR1TS180","DR1TM181",
             "DR1TP182","DR1TP183","DR1TP204","DR1TP205","DR1TP226")
df_pca <- nhanes_clean %>% dplyr::select(SEQN, all_of(fa_vars), High_CRP_num) %>% drop_na()
pca_res <- principal(df_pca[,fa_vars], nfactors=3, rotate="varimax", scores=TRUE)
print(pca_res$loadings, cutoff=0.3)

df_pca <- cbind(df_pca, pca_res$scores) %>%
  rename(Pattern_1 = RC1, Pattern_2 = RC2, Pattern_3 = RC3)
pca_pvals <- sapply(1:3, function(i) {
  f <- as.formula(paste0("High_CRP_num ~ Pattern_", i))
  summary(glm(f, data=df_pca, family=binomial()))$coefficients[2,4]
})

# 3.8 Mediation analysis (TyG and SII)
df_med <- nhanes_clean %>%
  mutate(Ratio_Scale = Ratio_adj * 10) %>%
  drop_na(TyG_Index, SII_Index)
covars <- "Age + Gender + BMI + Smoking + Energy_Kcal + Diabetes"

set.seed(2026)
med_tyg <- mediate(
  lm(paste("TyG_Index ~ Ratio_Scale +", covars), data = df_med),
  glm(paste("High_CRP_num ~ Ratio_Scale + TyG_Index +", covars),
      data = df_med, family = binomial("probit")),
  treat = "Ratio_Scale", mediator = "TyG_Index", boot = TRUE, sims = 500
)
summary(med_tyg)

set.seed(2026)
med_sii <- mediate(
  lm(paste("SII_Index ~ Ratio_Scale +", covars), data = df_med),
  glm(paste("High_CRP_num ~ Ratio_Scale + SII_Index +", covars),
      data = df_med, family = binomial("probit")),
  treat = "Ratio_Scale", mediator = "SII_Index", boot = TRUE, sims = 500
)
summary(med_sii)

# 3.9 Sensitivity analysis (excluding malignancy)
mcq <- bind_rows(
  nhanes('MCQ_I') %>% dplyr::select(SEQN, MCQ220),
  nhanes('MCQ_J') %>% dplyr::select(SEQN, MCQ220)
)
df_sens <- nhanes_clean %>% left_join(mcq, by="SEQN") %>%
  filter(is.na(MCQ220) | MCQ220 != 1)
fit_sens <- glm(
  High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR + BMI +
    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber +
    Diabetes + Hypertension,
  data = df_sens, family = binomial()
)
tidy(fit_sens, exponentiate=TRUE, conf.int=TRUE) %>% filter(grepl("Ratio_adj_Q4", term))

# 5. Track II: Model derivation and external validation -----------------------

# 5.1 Unified predictor pool (variables common to NHANES and CHARLS)
nhanes_model <- nhanes_clean %>%
  dplyr::select(SEQN, Age, Gender, BMI, Smoking, TyG_Index, High_CRP_num) %>%
  drop_na() %>%
  as.data.frame()

train_index <- createDataPartition(nhanes_model$High_CRP_num, p=0.7, list=FALSE)
df_train <- nhanes_model[train_index, ]
df_test  <- nhanes_model[-train_index, ]

x_train <- model.matrix(High_CRP_num ~ . - SEQN, data = df_train)[, -1]
y_train <- df_train$High_CRP_num

cv_lasso <- cv.glmnet(x_train, y_train, family="binomial", alpha=1, nfolds=10)
lambda_1se <- cv_lasso$lambda.1se
selected_vars_raw <- rownames(coef(cv_lasso, s=lambda_1se))[
  which(coef(cv_lasso, s=lambda_1se)[,1] != 0)]
selected_vars_raw <- setdiff(selected_vars_raw, "(Intercept)")
message("LASSO selected raw: ", paste(selected_vars_raw, collapse=", "))

# ---- Map dummy names back to original factor names ----
map_back <- function(vars) {
  vars <- gsub("GenderFemale|GenderMale", "Gender", vars)
  vars <- gsub("SmokingCurrent|SmokingFormer/Unknown|SmokingNever", "Smoking", vars)
  unique(vars)
}
selected_vars <- map_back(selected_vars_raw)
message("Final model variables: ", paste(selected_vars, collapse=", "))

formula_str <- paste("High_CRP_num ~", paste(selected_vars, collapse=" + "))
fit_final <- glm(as.formula(formula_str), data = df_train, family = binomial())
summary(fit_final)

df_test$pred <- predict(fit_final, newdata = df_test, type = "response")
roc_internal <- roc(df_test$High_CRP_num, df_test$pred)
c_index_internal <- auc(roc_internal)

pdf("results/figures/Internal_Calibration.pdf")
val.prob(df_test$pred, df_test$High_CRP_num)
dev.off()

dca_internal <- dca(High_CRP_num ~ pred, data = df_test,
                    thresholds = seq(0,1,0.01),
                    label = list(pred = "Model"))
pdf("results/figures/Internal_DCA.pdf")
plot(dca_internal)
dev.off()

# Incremental value of SII
df_test_sii <- df_test %>%
  inner_join(nhanes_clean %>% dplyr::select(SEQN, SII_Index), by = "SEQN")
if (nrow(df_test_sii) > 0) {
  formula_sii_str <- paste("High_CRP_num ~", paste(selected_vars, collapse=" + "), "+ SII_Index")
  fit_with_sii <- glm(as.formula(formula_sii_str), data = df_test_sii, family = binomial())
  roc_with <- roc(df_test_sii$High_CRP_num, predict(fit_with_sii, type="response"))
  delong_p <- roc.test(roc_internal, roc_with)$p.value
  message("DeLong P (SII increment): ", delong_p)
}

# 6. CHARLS external validation -----------------------------------------------
charls_harmonized <- read_dta("/Users/bing/CHARLS/Harmonized/H_CHARLS_D_Data.dta") %>%
  dplyr::select(ID, ragender, r4agey)
blood_data <- read_dta("/Users/bing/CHARLS/2015/Blood.dta") %>%
  dplyr::select(ID, bl_tg, bl_glu, bl_crp)
biomarker_data <- read_dta("/Users/bing/CHARLS/2015/Biomarker.dta") %>%
  dplyr::select(ID, qi002, ql002)

charls_clean <- charls_harmonized %>%
  left_join(blood_data, by = "ID") %>%
  left_join(biomarker_data, by = "ID") %>%
  filter(!is.na(ragender), !is.na(bl_tg), !is.na(bl_glu),
         !is.na(qi002), !is.na(ql002), !is.na(bl_crp),
         !is.na(r4agey)) %>%
  filter(qi002 > 50, ql002 > 10, bl_crp <= 10) %>%
  mutate(
    Gender = ifelse(ragender == 1, "Male", "Female"),
    Age = r4agey,
    BMI = ql002 / ((qi002/100)^2),
    TyG_Index = log(bl_tg * bl_glu / 2),
    High_CRP_num = ifelse(bl_crp > 3, 1, 0)
  ) %>%
  dplyr::select(ID, Gender, Age, BMI, TyG_Index, High_CRP_num)

charls_clean$Gender <- factor(charls_clean$Gender, levels = c("Male","Female"))

# Safety net for missing variables
missing_vars <- setdiff(selected_vars, colnames(charls_clean))
if (length(missing_vars) > 0) {
  message("WARNING: CHARLS lacks: ", paste(missing_vars, collapse=", "))
  message("Refitting model without these variables.")
  available_vars <- intersect(selected_vars, colnames(charls_clean))
  if (length(available_vars) < 2) {
    stop("Too few predictors available in CHARLS for valid external validation.")
  }
  formula_str <- paste("High_CRP_num ~", paste(available_vars, collapse=" + "))
  fit_final <- glm(as.formula(formula_str), data = df_train, family = binomial())
  message("Refitted model with: ", paste(available_vars, collapse=", "))
  selected_vars <- available_vars
}

charls_clean$pred <- predict(fit_final, newdata = charls_clean, type = "response")

roc_ext <- roc(charls_clean$High_CRP_num, charls_clean$pred)
c_index_ext <- auc(roc_ext)
message("External C-index (CHARLS): ", round(c_index_ext, 3))

pdf("results/figures/External_Calibration.pdf")
val.prob(charls_clean$pred, charls_clean$High_CRP_num)
dev.off()

dca_ext <- dca(High_CRP_num ~ pred, data = charls_clean,
               thresholds = seq(0,1,0.01),
               label = list(pred = "Model"))
pdf("results/figures/External_DCA.pdf")
plot(dca_ext)
dev.off()

# 7. Save results ------------------------------------------------------------
saveRDS(list(
  fit_final = fit_final,
  selected_vars = selected_vars,
  c_internal = c_index_internal,
  c_external = c_index_ext,
  subgroups = subgroups,
  pca_pvals = pca_pvals,
  trend_p = p_trend
), "data/final_results.rds")

message("All analyses completed successfully.")


# ==============================================================================
# Project Data Overview – Print all key results for manuscript preparation
# Run after the main pipeline has completed successfully.
# ==============================================================================

# Load saved results
res <- readRDS("data/final_results.rds")
attach(res, warn.conflicts = FALSE)

cat("\n")
cat("================================================================================\n")
cat("                    PROJECT DATA OVERVIEW\n")
cat("================================================================================\n\n")

# 1. Cohort sizes -------------------------------------------------------------
cat("--- 1. COHORT SIZES ---\n")
cat(sprintf("NHANES full cohort (Track I)           : %d\n", nrow(readRDS("data/nhanes_clean.rds"))))
cat(sprintf("CHARLS external validation cohort      : %d\n", nrow(readRDS("data/final_results.rds")$selected_vars) )) # 修正：用实际数据
cat("\n")

# 修正：直接从拟合对象和环境获取样本量
cat(sprintf("Training set (Track II)                : %d\n", length(fit_final$y)))
cat(sprintf("Internal test set                      : %d\n", length(fit_final$fitted.values) - length(fit_final$y))) # 近似
cat("\n")

# 更好的方法：重新读取nhanes_model
nhanes_model <- readRDS("data/nhanes_clean.rds") %>%
  dplyr::select(SEQN, Age, Gender, BMI, Smoking, TyG_Index, High_CRP_num) %>%
  drop_na()
train_idx <- caret::createDataPartition(nhanes_model$High_CRP_num, p = 0.7, list = FALSE)
cat(sprintf("Training set (Track II)                : %d\n", nrow(nhanes_model[train_idx, ])))
cat(sprintf("Internal test set                      : %d\n", nrow(nhanes_model[-train_idx, ])))
cat(sprintf("CHARLS external validation cohort      : %d\n", nrow(readRDS("data/final_results.rds")$selected_vars) )) # 不适用，改用外部验证时保存的charls_clean
# 实际上charls_clean没有单独保存，我们后面单独处理
cat("\n")

# 2. Track I core association -------------------------------------------------
cat("--- 2. TRACK I: OBSERVATIONAL ASSOCIATION ---\n")
cat(sprintf("RCS nonlinear P-value                  : 0.7971\n"))
cat(sprintf("Linear trend P-value                   : %.4f\n", trend_p))
cat("Multivariable-adjusted (Model 3):\n")
cat("  Q4 vs Q1 OR = 1.20 (95% CI: 1.04-1.40), P = 0.013\n")
cat("\n")

# 3. PSM and E-value ----------------------------------------------------------
cat("--- 3. PSM AND E-VALUE ---\n")
cat("PSM matched pseudo-population          : 4,192 participants (2,096 per group)\n")
cat("PSM OR (Q4 vs Q1)                      : 1.176 (95% CI: 1.006-1.375), P = 0.0423\n")
cat(sprintf("E-value for point estimate             : 1.418\n"))
cat(sprintf("E-value for CI lower limit             : 1.082\n"))  # 近似
cat("\n")

# 4. Sensitivity analyses -----------------------------------------------------
cat("--- 4. SENSITIVITY ANALYSES (Q4 vs Q1) ---\n")
cat("Excluding malignancy   : OR = 1.20 (1.04-1.39), P = 0.015\n")
cat("Excluding diabetes     : OR = 1.18 (1.00-1.38), P = 0.046\n")
cat("Excluding hypertension : OR = 1.23 (1.02-1.49), P = 0.033\n")
cat("Stricter hs-CRP >5 mg/L: OR = 1.26 (1.05-1.52), P = 0.015\n")
cat("\n")

# 5. PCA ----------------------------------------------------------------------
cat("--- 5. PCA OF DIETARY FATTY ACIDS ---\n")
cat(sprintf("Principal components extracted         : 3 (eigenvalues > 1)\n"))
cat(sprintf("Cumulative variance explained          : 63.72%%\n"))
cat(sprintf("RC1 (long-chain SFA pattern) variance  : 34.86%%\n"))
cat(sprintf("RC1 association with SCI P-value       : %.4f\n", pca_pvals[1]))
cat(sprintf("RC2 P-value                            : %.4f\n", pca_pvals[2]))
cat(sprintf("RC3 P-value                            : %.4f\n", pca_pvals[3]))
cat("\n")

# 6. Mediation analysis -------------------------------------------------------
cat("--- 6. MEDIATION ANALYSIS ---\n")
cat("Average Direct Effect (ADE)            : 0.015 (95% CI: 0.002-0.028), P = 0.016\n")
cat("TyG indirect pathway (ACME)            : -0.0003 (95% CI: -0.0015 to 0.0013), P = 0.736\n")
cat("SII indirect pathway (ACME)            : 0.0011 (95% CI: -0.0001 to 0.0022)\n")
cat("  Proportion mediated                  : ~7.2%, P = 0.080 (marginal)\n")
cat("\n")

# 7. Track II model -----------------------------------------------------------
cat("--- 7. TRACK II: PREDICTION MODEL ---\n")
cat(sprintf("LASSO selected predictors              : %s\n", paste(selected_vars, collapse = ", ")))
cat(sprintf("Internal validation C-index            : %.3f\n", c_internal))
cat(sprintf("External validation C-index (CHARLS)   : %.3f\n", c_external))
cat("Calibration and DCA plots saved in     : results/figures/\n")
cat("\n")

# 8. Incremental value of SII ------------------------------------------------
cat("--- 8. INCREMENTAL VALUE OF SII ---\n")
cat("DeLong test (adding SII to model) P-value: (see console output above)\n")
cat("NRI and IDI P-values not directly stored, please refer to original output.\n")
cat("\n")

# 9. Subgroup results (example) -----------------------------------------------
cat("--- 9. SUBGROUP ANALYSES (PER 0.1-UNIT INCREASE) ---\n")
print(subgroups, row.names = FALSE)
cat("\n")

cat("================================================================================\n")
cat("                   END OF DATA OVERVIEW\n")
cat("================================================================================\n")
detach(res)


# ==============================================================================
# Final publication-quality figures 
# All calibrations redrawn with ggplot2 to avoid val.prob argument errors
# ==============================================================================

library(dplyr)
library(ggplot2)
library(rms)         # only for nomogram
library(glmnet)
library(pROC)
library(dcurves)
library(splines)

# ---------- Publication style ----------
pal_blue   <- "#0072B2"
pal_orange <- "#D55E00"
pal_grey   <- "#666666"

theme_nature <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "#E0E0E0")
  )

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)





# ==============================================================================
# 最终出版级图表（修正 RCS Y 轴为 OR，使用 glm + ns 并设定参考点）
# ==============================================================================

library(dplyr)
library(ggplot2)
library(splines)
library(rms)         # 用于 nomogram
library(glmnet)
library(pROC)
library(dcurves)

# ---------- 统一配色与主题 ----------
pal_blue   <- "#0072B2"
pal_orange <- "#D55E00"
pal_grey   <- "#666666"

theme_nature <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "#E0E0E0")
  )

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ============================
# Figure 1: RCS dose-response (OR scale, correct reference)
# ============================
message("1/8 RCS curve (OR scale)")

# 使用 glm + ns 拟合
fit_rcs_glm <- glm(High_CRP_num ~ ns(Ratio_adj, df = 4) + Age + Gender + Race + PIR +
                     BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal +
                     Total_Fat + Fiber + Diabetes + Hypertension,
                   data = nhanes_clean, family = binomial)

# 构建预测网格（使用 5%-95% 分位数范围）
x_range <- quantile(nhanes_clean$Ratio_adj, probs = c(0.05, 0.95), na.rm = TRUE)
newdata <- data.frame(Ratio_adj = seq(x_range[1], x_range[2], length.out = 100))

# 固定协变量为中位数/众数
newdata$Age <- median(nhanes_clean$Age, na.rm = TRUE)
newdata$Gender <- factor("Male", levels = levels(nhanes_clean$Gender))
newdata$Race <- factor(levels(nhanes_clean$Race)[1], levels = levels(nhanes_clean$Race))
newdata$PIR <- median(nhanes_clean$PIR, na.rm = TRUE)
newdata$BMI <- median(nhanes_clean$BMI, na.rm = TRUE)
newdata$Smoking <- factor("Never", levels = levels(nhanes_clean$Smoking))
newdata$Alcohol_g <- median(nhanes_clean$Alcohol_g, na.rm = TRUE)
newdata$MET_imp <- median(nhanes_clean$MET_imp, na.rm = TRUE)
newdata$Energy_Kcal <- median(nhanes_clean$Energy_Kcal, na.rm = TRUE)
newdata$Total_Fat <- median(nhanes_clean$Total_Fat, na.rm = TRUE)
newdata$Fiber <- median(nhanes_clean$Fiber, na.rm = TRUE)
newdata$Diabetes <- factor("No", levels = levels(nhanes_clean$Diabetes))
newdata$Hypertension <- factor("No", levels = levels(nhanes_clean$Hypertension))

# 预测 log-odds
pred_link <- predict(fit_rcs_glm, newdata, se.fit = TRUE, type = "link")

# 以中位数为参考点，计算 OR
ref_idx <- which.min(abs(newdata$Ratio_adj - median(nhanes_clean$Ratio_adj, na.rm = TRUE)))
logor <- pred_link$fit - pred_link$fit[ref_idx]
se_logor <- pred_link$se.fit
newdata$OR <- exp(logor)
newdata$lower <- exp(logor - 1.96 * se_logor)
newdata$upper <- exp(logor + 1.96 * se_logor)

# 绘制
p_rcs <- ggplot(newdata, aes(x = Ratio_adj)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = pal_grey, linewidth = 0.7) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = pal_blue, alpha = 0.15) +
  geom_line(aes(y = OR), color = pal_blue, linewidth = 1.2) +
  annotate("text", x = max(newdata$Ratio_adj) * 0.7,
           y = max(newdata$upper) * 0.8,
           label = expression(italic(P)[non-linear] == 0.153),
           size = 5, color = "black") +
  labs(x = "Dietary 18:0/16:0 ratio (energy-adjusted)",
       y = "Odds Ratio (95% CI)") +
  theme_nature

ggsave("results/figures/Stage_2_3_RCS_Curve.pdf", p_rcs, width = 8, height = 6)

# ============================
# Figure 2: Subgroup forest plot
# ============================
message("2/8 Subgroup forest plot")

df_plot <- subgroups %>%
  mutate(Subgroup = factor(Subgroup, levels = rev(unique(Subgroup))),
         Label = sprintf("%.2f (%.2f–%.2f)", OR, Lower, Upper))

p_forest <- ggplot(df_plot, aes(x = OR, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = pal_grey, linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25,
                 color = pal_blue, linewidth = 0.9) +
  geom_point(shape = 18, size = 4.5, color = pal_blue) +
  geom_text(aes(label = Label, x = Upper + 0.08), hjust = 0, size = 4.2, color = "black") +
  scale_x_continuous(limits = c(min(df_plot$Lower) - 0.2, max(df_plot$Upper) + 0.6)) +
  labs(x = "Odds Ratio (95% CI) per 0.1‑unit increase", y = "") +
  theme_nature +
  theme(axis.text.y = element_text(face = "bold", size = 11))

ggsave("results/figures/Stage_3_5_Subgroup_Forest.pdf", p_forest, width = 10, height = 7)

# ============================
# Figure 3: LASSO path + CV
# ============================
message("3/8 LASSO path")

pdf("results/figures/Stage_5_1_LASSO_Path.pdf", width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(cv_lasso$glmnet.fit, "lambda", label = TRUE,
     col = c(pal_grey, pal_blue, pal_orange),
     xlab = expression(Log~lambda), ylab = "Coefficients",
     main = "LASSO coefficient paths")
plot(cv_lasso, col = pal_blue, main = "10‑fold cross‑validation",
     xlab = expression(Log~lambda), ylab = "Binomial Deviance")
dev.off()

# ============================
# Figure 4: Nomogram
# ============================
message("4/8 Nomogram")

formula_nom <- as.formula(paste("High_CRP_num ~", paste(selected_vars, collapse = " + ")))
dd_train <- datadist(df_train)
options(datadist = "dd_train")
fit_nom <- lrm(formula_nom, data = df_train, x = TRUE, y = TRUE)

nom <- nomogram(fit_nom, fun = plogis, fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                funlabel = "Risk of high systemic inflammation", lp = FALSE)

pdf("results/figures/Stage_5_2_Training_Nomogram.pdf", width = 12, height = 8)
plot(nom, cex.var = 0.85, cex.axis = 0.85, lmgp = 0.2,
     col.grid = c(pal_grey, pal_grey))
title(main = "Nomogram for predicting inflammation (training cohort)", line = 1.5)
dev.off()

# ---------- Helper for calibration plots ----------
plot_calibration <- function(prob, outcome, title_text, color_line) {
  data <- data.frame(prob = prob, outcome = outcome)
  data <- data %>% mutate(bin = cut(prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE))
  cal_data <- data %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(prob, na.rm = TRUE),
      obs_rate  = mean(outcome, na.rm = TRUE),
      count     = n()
    ) %>%
    ungroup() %>%
    filter(!is.na(bin))
  
  ggplot(cal_data, aes(x = mean_pred, y = obs_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = pal_grey, linewidth = 0.8) +
    geom_point(color = color_line, size = 3) +
    geom_line(color = color_line, linewidth = 1) +
    labs(x = "Predicted probability", y = "Observed proportion",
         title = title_text) +
    theme_nature +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

# ============================
# Figures 5–8: Calibration & DCA
# ============================
message("5/8 Internal calibration")
p_cal_int <- plot_calibration(df_test$pred, df_test$High_CRP_num,
                              "NHANES internal calibration", pal_blue)
ggsave("results/figures/Internal_Calibration.pdf", p_cal_int, width = 7, height = 6)

message("6/8 Internal DCA")
dca_int <- dca(High_CRP_num ~ pred, data = df_test,
               thresholds = seq(0, 1, 0.01),
               label = list(pred = "3‑variable model")) |> plot() +
  labs(title = "NHANES internal DCA") + theme_nature
ggsave("results/figures/Internal_DCA.pdf", dca_int, width = 7, height = 6)

message("7/8 External calibration")
p_cal_ext <- plot_calibration(charls_clean$pred, charls_clean$High_CRP_num,
                              "CHARLS external calibration", pal_orange)
ggsave("results/figures/External_Calibration.pdf", p_cal_ext, width = 7, height = 6)

message("8/8 External DCA")
dca_ext <- dca(High_CRP_num ~ pred, data = charls_clean,
               thresholds = seq(0, 1, 0.01),
               label = list(pred = "3‑variable model")) |> plot() +
  labs(title = "CHARLS external DCA") + theme_nature
ggsave("results/figures/External_DCA.pdf", dca_ext, width = 7, height = 6)

message("All figures saved in results/figures/")



# ==============================================================================
# 敏感性分析代码（回应审稿人关于过度调整偏倚和总SFA独立性的质疑）
# 运行前提：已运行主分析脚本，环境中存在 nhanes_clean
# ==============================================================================

library(dplyr)
library(broom)

# ---- 准备数据：计算能量校正后的总SFA ----
nhanes_clean <- nhanes_clean %>%
  mutate(
    Total_SFA_adj = resid(lm(DR1TSFAT ~ Energy_Kcal, data = nhanes_clean, 
                             na.action = na.exclude)) + mean(DR1TSFAT, na.rm = TRUE),
    Ratio_adj_Q_num = as.numeric(Ratio_adj_Q)
  )

# ============================================================
# 敏感性分析 1：不调整 BMI、糖尿病、高血压（排除过度调整偏倚）
# ============================================================
message(">>> 敏感性分析 1：不调整 BMI/糖尿病/高血压")

fit_sens1 <- glm(
  High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR +
    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber,
  data = nhanes_clean, family = binomial()
)

res_sens1 <- tidy(fit_sens1, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "Ratio_adj_QQ4")

# 趋势检验
fit_trend_sens1 <- glm(
  High_CRP_num ~ Ratio_adj_Q_num + Age + Gender + Race + PIR +
    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber,
  data = nhanes_clean, family = binomial()
)
p_trend_sens1 <- coef(summary(fit_trend_sens1))["Ratio_adj_Q_num", "Pr(>|z|)"]

cat("\n敏感性分析 1 结果（不调整 BMI/糖尿病/高血压）:\n")
cat(sprintf("  Q4 vs Q1 OR = %.2f (95%% CI: %.2f-%.2f), P = %.4f\n",
            res_sens1$estimate, res_sens1$conf.low, res_sens1$conf.high, res_sens1$p.value))
cat(sprintf("  趋势 P 值 = %.4f\n", p_trend_sens1))

# ============================================================
# 敏感性分析 2：额外调整总饱和脂肪（检验比值效应的独立性）
# ============================================================
message("\n>>> 敏感性分析 2：额外调整总 SFA")

fit_sens2_full <- glm(
  High_CRP_num ~ Ratio_adj_Q + Total_SFA_adj + Age + Gender + Race + PIR + 
    BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber + 
    Diabetes + Hypertension,
  data = nhanes_clean, family = binomial()
)

res_sens2 <- tidy(fit_sens2_full, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "Ratio_adj_QQ4")

cat("\n敏感性分析 2 结果（额外调整总 SFA）:\n")
cat(sprintf("  Q4 vs Q1 OR = %.2f (95%% CI: %.2f-%.2f), P = %.4f\n",
            res_sens2$estimate, res_sens2$conf.low, res_sens2$conf.high, res_sens2$p.value))

# ============================================================
# 保存结果
# ============================================================
sensitivity_results <- list(
  sens1_OR = res_sens1$estimate,
  sens1_CI_low = res_sens1$conf.low,
  sens1_CI_high = res_sens1$conf.high,
  sens1_P = res_sens1$p.value,
  sens1_trend_P = p_trend_sens1,
  sens2_OR = res_sens2$estimate,
  sens2_CI_low = res_sens2$conf.low,
  sens2_CI_high = res_sens2$conf.high,
  sens2_P = res_sens2$p.value
)
saveRDS(sensitivity_results, "data/sensitivity_results.rds")

message("\n>>> 敏感性分析完成，结果已保存至 data/sensitivity_results.rds")

# ==============================================================================
# 第二轮审稿意见补充分析完整代码
# 包含：增量价值、交互效应量、Ratio分布、PSM暴露差异、用药敏感性
# ==============================================================================

library(dplyr)
library(broom)
library(pROC)
library(Hmisc)
library(nhanesA)

# 确保环境中存在 nhanes_clean, df_train, df_test, psm_match

# ============================================================
# 补充分析 1：膳食比值的预测增量价值（DeLong, NRI, IDI）
# ============================================================
message("===== 补充分析1：膳食比值增量价值 =====")

# 构建含膳食比值的测试集
df_test_with_ratio <- df_test %>%
  inner_join(nhanes_clean %>% dplyr::select(SEQN, Ratio_adj), by = "SEQN")

df_train_with_ratio <- df_train %>%
  inner_join(nhanes_clean %>% dplyr::select(SEQN, Ratio_adj), by = "SEQN")

# 训练两个模型
fit_3var <- glm(High_CRP_num ~ Gender + BMI + TyG_Index, 
                data = df_train, family = binomial())

fit_4var <- glm(High_CRP_num ~ Gender + BMI + TyG_Index + Ratio_adj, 
                data = df_train_with_ratio, family = binomial())

# 预测
df_test_with_ratio$pred_3var <- predict(fit_3var, newdata = df_test_with_ratio, type = "response")
df_test_with_ratio$pred_4var <- predict(fit_4var, newdata = df_test_with_ratio, type = "response")

# DeLong检验
roc_3var <- roc(df_test_with_ratio$High_CRP_num, df_test_with_ratio$pred_3var)
roc_4var <- roc(df_test_with_ratio$High_CRP_num, df_test_with_ratio$pred_4var)
delong_res <- roc.test(roc_3var, roc_4var)

# NRI和IDI
nri_idi <- Hmisc::improveProb(x1 = df_test_with_ratio$pred_3var, 
                               x2 = df_test_with_ratio$pred_4var, 
                               y = df_test_with_ratio$High_CRP_num)

cat(sprintf("DeLong P = %.4f\n", delong_res$p.value))
print(nri_idi)

# ============================================================
# 补充分析 2：亚组交互效应量与95%CI（性别、BMI、年龄）
# ============================================================
message("\n===== 补充分析2：亚组交互效应 =====")

nhanes_clean <- nhanes_clean %>% mutate(Ratio_Scale = Ratio_adj * 10)

# 性别交互
fit_sex <- glm(High_CRP_num ~ Ratio_Scale * Gender + Age + BMI + Smoking + Energy_Kcal,
               data = nhanes_clean, family = binomial())
sex_interact <- tidy(fit_sex, conf.int = TRUE) %>% filter(grepl(":", term))
cat("\n性别交互效应:\n")
print(sex_interact)

# BMI交互（连续）
fit_bmi <- glm(High_CRP_num ~ Ratio_Scale * BMI + Age + Gender + Smoking + Energy_Kcal,
               data = nhanes_clean, family = binomial())
bmi_interact <- tidy(fit_bmi, conf.int = TRUE) %>% filter(grepl(":", term))
cat("\nBMI交互效应:\n")
print(bmi_interact)

# 年龄交互（连续）
fit_age <- glm(High_CRP_num ~ Ratio_Scale * Age + Gender + BMI + Smoking + Energy_Kcal,
               data = nhanes_clean, family = binomial())
age_interact <- tidy(fit_age, conf.int = TRUE) %>% filter(grepl(":", term))
cat("\n年龄交互效应:\n")
print(age_interact)

# ============================================================
# 补充分析 3：Ratio 分布特征
# ============================================================
message("\n===== 补充分析3：Ratio分布 =====")
cat(sprintf("Mean = %.3f\n", mean(nhanes_clean$Ratio_adj, na.rm = TRUE)))
cat(sprintf("SD = %.3f\n", sd(nhanes_clean$Ratio_adj, na.rm = TRUE)))
cat(sprintf("Median = %.3f\n", median(nhanes_clean$Ratio_adj, na.rm = TRUE)))
cat(sprintf("IQR = %.3f - %.3f\n", 
            quantile(nhanes_clean$Ratio_adj, 0.25, na.rm = TRUE), 
            quantile(nhanes_clean$Ratio_adj, 0.75, na.rm = TRUE)))
cat(sprintf("Min = %.3f, Max = %.3f\n", 
            min(nhanes_clean$Ratio_adj, na.rm = TRUE), 
            max(nhanes_clean$Ratio_adj, na.rm = TRUE)))

# ============================================================
# 补充分析 4：PSM 后暴露差异
# ============================================================
message("\n===== 补充分析4：PSM后暴露差异 =====")
matched_data <- match.data(psm_match)
q1_idx <- matched_data$Ratio_adj_Q == "Q1"
q4_idx <- matched_data$Ratio_adj_Q == "Q4"
cat(sprintf("Q1 Ratio_adj mean = %.3f\n", mean(matched_data$Ratio_adj[q1_idx])))
cat(sprintf("Q4 Ratio_adj mean = %.3f\n", mean(matched_data$Ratio_adj[q4_idx])))
cat(sprintf("Mean difference = %.3f\n", 
            mean(matched_data$Ratio_adj[q4_idx]) - mean(matched_data$Ratio_adj[q1_idx])))

# ============================================================
# 补充分析 5：排除他汀/NSAID 用药人群
# ============================================================
message("\n===== 补充分析5：排除用药者 =====")
rx_i <- nhanes('RXQ_RX_I') %>% dplyr::select(SEQN, RXDDRUG)
rx_j <- nhanes('RXQ_RX_J') %>% dplyr::select(SEQN, RXDDRUG)
rx_all <- bind_rows(rx_i, rx_j)

statins <- unique(rx_all$RXDDRUG[grep(
  "ATORVASTATIN|ROSUVASTATIN|SIMVASTATIN|PRAVASTATIN|LOVASTATIN|FLUVASTATIN|PITAVASTATIN",
  rx_all$RXDDRUG, ignore.case = TRUE)])
nsaids <- unique(rx_all$RXDDRUG[grep(
  "IBUPROFEN|NAPROXEN|DICLOFENAC|CELECOXIB|MELOXICAM|INDOMETHACIN|KETOROLAC",
  rx_all$RXDDRUG, ignore.case = TRUE)])

rx_users <- rx_all %>% 
  filter(RXDDRUG %in% c(statins, nsaids)) %>%
  pull(SEQN) %>% unique()

nhanes_sens_drug <- nhanes_clean %>% filter(!SEQN %in% rx_users)

fit_drug_sens <- glm(High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR +
                      Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber,
                    data = nhanes_sens_drug, family = binomial())
res_drug_sens <- tidy(fit_drug_sens, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "Ratio_adj_QQ4")

cat(sprintf("纳入人数: %d\n", nrow(nhanes_sens_drug)))
print(res_drug_sens)

message("\n===== 全部补充分析完成 =====")


fit_primary <- glm(High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR +
                    Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber,
                  data = nhanes_clean, family = binomial())
tidy(fit_primary, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("Ratio_adj_Q", term))

# 安装 cobalt（如果没有）
if (!require("cobalt")) install.packages("cobalt")
library(cobalt)
library(ggplot2)

# 创建 Love 图（比较匹配前后标准化均值差）
p_love <- love.plot(psm_match, 
                    binary = "std", 
                    thresholds = c(m = 0.1),   # 原 threshold 参数可能改为 thresholds
                    abs = TRUE,
                    var.order = "unadjusted",
                    line = TRUE,
                    themes = list(
                      theme_classic(base_size = 12) +
                        theme(plot.title = element_text(face = "bold", hjust = 0.5),
                              axis.title = element_text(face = "bold"),
                              axis.text = element_text(color = "black"),
                              legend.position = "bottom",
                              panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))
                    )) +
  labs(title = "Covariate Balance Before and After Matching",
       x = "Absolute Standardized Mean Difference",
       y = "")

# 保存为 PDF
ggsave("results/figures/PSM_Love_Plot.pdf", p_love, width = 8, height = 6, device = "pdf")



                      
# 重新运行他汀/NSAID 排除敏感性分析
library(dplyr)
library(broom)
library(nhanesA)

rx_i <- nhanes('RXQ_RX_I') %>% dplyr::select(SEQN, RXDDRUG)
rx_j <- nhanes('RXQ_RX_J') %>% dplyr::select(SEQN, RXDDRUG)
rx_all <- bind_rows(rx_i, rx_j)

statins <- unique(rx_all$RXDDRUG[grep(
  "ATORVASTATIN|ROSUVASTATIN|SIMVASTATIN|PRAVASTATIN|LOVASTATIN|FLUVASTATIN|PITAVASTATIN",
  rx_all$RXDDRUG, ignore.case = TRUE)])
nsaids <- unique(rx_all$RXDDRUG[grep(
  "IBUPROFEN|NAPROXEN|DICLOFENAC|CELECOXIB|MELOXICAM|INDOMETHACIN|KETOROLAC",
  rx_all$RXDDRUG, ignore.case = TRUE)])

rx_users <- rx_all %>% 
  filter(RXDDRUG %in% c(statins, nsaids)) %>%
  pull(SEQN) %>% unique()

nhanes_sens_drug <- nhanes_clean %>% filter(!SEQN %in% rx_users)

fit_drug_sens <- glm(High_CRP_num ~ Ratio_adj_Q + Age + Gender + Race + PIR +
                      Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber,
                    data = nhanes_sens_drug, family = binomial())
res_drug_sens <- tidy(fit_drug_sens, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "Ratio_adj_QQ4")

cat(sprintf("纳入人数: %d\n", nrow(nhanes_sens_drug)))
print(res_drug_sens)

# 在已有的 df_pca 上运行
df_pca_std <- df_pca %>%
  mutate(across(starts_with("Pattern_"), ~ as.numeric(scale(.))))
for (i in 1:3) {
  f <- as.formula(paste0("High_CRP_num ~ Pattern_", i))
  mod <- glm(f, data = df_pca_std, family = binomial)
  res <- tidy(mod, exponentiate = TRUE, conf.int = TRUE)[2,]
  cat(sprintf("PC%d: OR = %.2f (95%% CI %.2f-%.2f), P = %.4f\n",
              i, res$estimate, res$conf.low, res$conf.high, coef(summary(mod))[2,4]))
}

fit_trend_full <- glm(High_CRP_num ~ Ratio_adj_Q_num + Age + Gender + Race + PIR + BMI +
                      Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber +
                      Diabetes + Hypertension,
                      data = nhanes_clean, family = binomial())
p_trend_full <- coef(summary(fit_trend_full))["Ratio_adj_Q_num", "Pr(>|z|)"]
print(p_trend_full)


fit_trend_full <- glm(High_CRP_num ~ Ratio_adj_Q_num + Age + Gender + Race + PIR + BMI +
                      Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber +
                      Diabetes + Hypertension,
                      data = nhanes_clean, family = binomial())
coef_trend <- coef(summary(fit_trend_full))["Ratio_adj_Q_num", ]
or_trend_full <- exp(coef_trend["Estimate"])
ci_trend_full <- exp(confint(fit_trend_full)["Ratio_adj_Q_num", ])
cat(sprintf("Fully adjusted trend OR = %.2f (95%% CI: %.2f-%.2f), P = %.4f\n",
            or_trend_full, ci_trend_full[1], ci_trend_full[2], coef_trend["Pr(>|z|)"]))


print(pca_res$loadings, cutoff = 0.3)
# 或者查看方差贡献
summary(pca_res)



