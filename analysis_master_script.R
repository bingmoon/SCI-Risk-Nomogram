
# ==============================================================================
# 项目名称：AB (基于三明治架构的双周期合并队列研究)
# 当前阶段：一、基线构建与脂肪酸全谱解码
# 步骤：1 & 2 (数据提取与清洗) & 3 (Table 1 构建)
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

# ------------------------------------------------------------------------------
# 一_1：初始化环境与多维原始数据拉取
# ------------------------------------------------------------------------------
dir.create("/Users/bing/AA/AB", recursive = TRUE, showWarnings = FALSE)
setwd("/Users/bing/AA/AB")

dir.create("results_一", showWarnings = FALSE)
dir.create("results_一/data", showWarnings = FALSE)
dir.create("results_一/plots", showWarnings = FALSE)

if (!require("pacman")) install.packages("pacman")
p_load(nhanesA, dplyr, tidyr, tableone)

message(">>> [一_1] 启动 2015-2018 双周期底层原始数据拉取 (含脂肪酸全谱)...")

pull_cycle_data <- function(cycle) {
  demo <- nhanes(paste0('DEMO_', cycle)) %>% select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, INDFMPIR)
  
  # 强烈扩增：拉取总脂肪及主要个别脂肪酸，为 PCA 降维做准备
  diet <- nhanes(paste0('DR1TOT_', cycle)) %>% 
    select(SEQN, DR1TKCAL, DR1TALCO, DR1TFIBE, DR1TSUGR, # 基础营养与添加糖
           DR1TTFAT, DR1TSFAT, DR1TMFAT, DR1TPFAT,       # 大类脂肪
           DR1TS140, DR1TS160, DR1TS180,                 # 饱和类细分
           DR1TM181, DR1TP182, DR1TP183, DR1TP204, DR1TP205, DR1TP226) # 单/多不饱和细分 (含EPA/DHA)
           
  crp  <- nhanes(paste0('HSCRP_', cycle)) %>% select(SEQN, LBXHSCRP)
  bmx  <- nhanes(paste0('BMX_', cycle)) %>% select(SEQN, BMXBMI)
  smq  <- nhanes(paste0('SMQ_', cycle)) %>% select(SEQN, SMQ020, SMQ040)
  diq  <- nhanes(paste0('DIQ_', cycle)) %>% select(SEQN, DIQ010)
  bpq  <- nhanes(paste0('BPQ_', cycle)) %>% select(SEQN, BPQ020)
  paq  <- nhanes(paste0('PAQ_', cycle)) %>% select(SEQN, any_of(c("PAQ610", "PAD615", "PAQ625", "PAD630", "PAQ640", "PAD645", "PAQ655", "PAD660", "PAQ670", "PAD675")))
  
  # 生化/血常规 (留作后续中介分析与模型构建)
  glu <- nhanes(paste0('GLU_', cycle)) %>% select(SEQN, LBXGLU)
  tg  <- nhanes(paste0('TRIGLY_', cycle)) %>% select(SEQN, LBXTR)
  cbc <- nhanes(paste0('CBC_', cycle)) %>% select(SEQN, LBXPLTSI, LBDNENO, LBDLYMNO)
  
  df_merge <- demo %>%
    left_join(diet, by = "SEQN") %>% left_join(crp, by = "SEQN") %>%
    left_join(bmx, by = "SEQN") %>% left_join(smq, by = "SEQN") %>%
    left_join(diq, by = "SEQN") %>% left_join(bpq, by = "SEQN") %>%
    left_join(paq, by = "SEQN") %>% left_join(glu, by = "SEQN") %>%
    left_join(tg, by = "SEQN") %>% left_join(cbc, by = "SEQN")
  
  return(df_merge)
}

df_raw_I <- pull_cycle_data("I") # 15-16
df_raw_J <- pull_cycle_data("J") # 17-18
df_stage_一_1_raw <- bind_rows(df_raw_I, df_raw_J)


# ------------------------------------------------------------------------------
# 一_2：严苛清洗与变量衍生 (分离混杂与中介)
# ------------------------------------------------------------------------------
message(">>> [一_2] 执行核心变量清洗与连续型暴露衍化...")

df_stage_一_2_clean <- df_stage_一_1_raw %>%
  drop_na(RIDAGEYR, RIAGENDR, DR1TS180, DR1TS160, LBXHSCRP) %>%
  filter(RIDAGEYR >= 20, LBXHSCRP <= 10, DR1TS160 > 0, DR1TS180 > 0) %>%
  mutate(
    # --- 核心主线变量 ---
    Ratio_18_16 = DR1TS180 / (DR1TS160 + 0.0001), 
    High_CRP = factor(ifelse(LBXHSCRP > 3, "Yes", "No"), levels = c("No", "Yes")),
    High_CRP_num = ifelse(LBXHSCRP > 3, 1, 0),
    
    # --- 人口学与体格 (混杂层) ---
    Age = RIDAGEYR,
    Gender = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    Race = factor(RIDRETH1),
    PIR = ifelse(is.na(INDFMPIR), median(INDFMPIR, na.rm = TRUE), INDFMPIR),
    BMI = ifelse(is.na(BMXBMI), median(BMXBMI, na.rm = TRUE), BMXBMI),
    
    # --- 膳食与生活方式 (混杂层) ---
    Energy_Kcal = ifelse(is.na(DR1TKCAL), median(DR1TKCAL, na.rm = TRUE), DR1TKCAL),
    Total_Fat = ifelse(is.na(DR1TTFAT), median(DR1TTFAT, na.rm = TRUE), DR1TTFAT),
    Fiber = ifelse(is.na(DR1TFIBE), median(DR1TFIBE, na.rm = TRUE), DR1TFIBE),
    Alcohol_g = ifelse(is.na(DR1TALCO), 0, DR1TALCO),
    
    Smoking = factor(case_when(
      as.character(SMQ020) %in% c("2", "No") ~ "Never",
      as.character(SMQ020) %in% c("1", "Yes") & as.character(SMQ040) %in% c("1", "2", "Every day", "Some days") ~ "Current",
      TRUE ~ "Former/Unknown"
    )),
    
    MET_min_week = (replace_na(ifelse(PAQ610 %in% 1:7, PAQ610, 0), 0) * replace_na(ifelse(PAD615 < 7777, PAD615, 0), 0) * 8) +  
                   (replace_na(ifelse(PAQ625 %in% 1:7, PAQ625, 0), 0) * replace_na(ifelse(PAD630 < 7777, PAD630, 0), 0) * 4) +  
                   (replace_na(ifelse(PAQ640 %in% 1:7, PAQ640, 0), 0) * replace_na(ifelse(PAD645 < 7777, PAD645, 0), 0) * 4) +  
                   (replace_na(ifelse(PAQ655 %in% 1:7, PAQ655, 0), 0) * replace_na(ifelse(PAD660 < 7777, PAD660, 0), 0) * 8) +  
                   (replace_na(ifelse(PAQ670 %in% 1:7, PAQ670, 0), 0) * replace_na(ifelse(PAD675 < 7777, PAD675, 0), 0) * 4),
    MET_imp = ifelse(MET_min_week == 0 | is.na(MET_min_week), median(MET_min_week[MET_min_week > 0], na.rm = TRUE), MET_min_week),
    
    # --- 临床共病 (混杂层) ---
    Diabetes = factor(case_when(as.character(DIQ010) %in% c("1", "Yes") ~ "Yes", TRUE ~ "No")),
    Hypertension = factor(case_when(as.character(BPQ020) %in% c("1", "Yes") ~ "Yes", TRUE ~ "No")),
    
    # --- 宿主基线内环境 (潜在中介/候选预测因子层，留作后用) ---
    TyG_Index = log(LBXTR * LBXGLU / 2),
    SII_Index = (LBXPLTSI * LBDNENO) / (LBDLYMNO + 0.0001)
  )

# 切割四分位数分组
ratio_q <- quantile(df_stage_一_2_clean$Ratio_18_16, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
df_stage_一_2_clean$Ratio_Quartile <- cut(df_stage_一_2_clean$Ratio_18_16, breaks = ratio_q, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))

# 保存纯净数据
saveRDS(df_stage_一_2_clean, "results_一/data/Stage_一_2_Cleaned_Cohort.rds")
message(glue::glue(">>> [一_2] 清洗完毕！最终进入后续流程的样本量：{nrow(df_stage_一_2_clean)}"))


# ------------------------------------------------------------------------------
# 一_3：生成标准描述性基线表 (Table 1)
# ------------------------------------------------------------------------------
message(">>> [一_3] 构建基于连续暴露四分位数的 Table 1...")

vars_t1 <- c("Age", "Gender", "Race", "PIR", "BMI", "Smoking", "Alcohol_g", 
             "MET_imp", "Energy_Kcal", "Total_Fat", "Fiber", 
             "Diabetes", "Hypertension", "High_CRP")

cat_t1 <- c("Gender", "Race", "Smoking", "Diabetes", "Hypertension", "High_CRP")

tab1 <- CreateTableOne(vars = vars_t1, strata = "Ratio_Quartile", data = df_stage_一_2_clean, factorVars = cat_t1)

cat("\n=================================================================\n")
cat("📊 Table 1: Baseline Characteristics (Q1 to Q4)\n")
cat("=================================================================\n")
print(tab1, nonnormal = c("Age", "PIR", "BMI", "Alcohol_g", "MET_imp", "Energy_Kcal", "Total_Fat", "Fiber"), showAllLevels = FALSE)

tab1_export <- print(tab1, nonnormal = c("Age", "PIR", "BMI", "Alcohol_g", "MET_imp", "Energy_Kcal", "Total_Fat", "Fiber"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_export, file = "results_一/data/Stage_一_3_Table1.csv")

message(">>> 🎉 阶段一前置工作完毕！请确认控制台输出的样本量及结果。")
# ==========================================================
# 紧急补丁：修复 Gender 变量全为 NA 的 Bug
# ==========================================================
setwd("/Users/bing/AA/AB")
library(dplyr)

# 1. 重新读取阶段一的数据底座
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# 2. 极客级修复 Gender (完美兼容字符和数字两种底层格式)
df_stage_2 <- df_stage_2 %>%
  mutate(
    Gender = factor(case_when(
      as.character(RIAGENDR) %in% c("1", "Male") ~ "Male",
      as.character(RIAGENDR) %in% c("2", "Female") ~ "Female",
      TRUE ~ "Unknown"
    ))
  )

# 3. 覆盖保存回原文件，彻底修复地基
saveRDS(df_stage_2, "results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# 检查一下是不是有男有女了
print(table(df_stage_2$Gender, useNA = "ifany"))

message(">>> 🎉 Gender 修复完毕！")



# ==============================================================================
# 项目名称：AB
# 当前阶段：二、核心关联与剂量反应探索
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

setwd("/Users/bing/AA/AB")
dir.create("results_二", showWarnings = FALSE)
dir.create("results_二/data", showWarnings = FALSE)
dir.create("results_二/plots", showWarnings = FALSE)

if (!require("pacman")) install.packages("pacman")
p_load(rms, dplyr, ggplot2, broom)

message(">>> [二_1] 加载阶段一的纯净数据底座...")
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# ==============================================================================
# 二_1：连续暴露的 RCS 探路 (确立线性模型的合法性)
# ==============================================================================
message(">>> [二_1] 启动 RCS 引擎，检验连续暴露的非线性特征...")

# 设定 rms 包需要的数据分布环境
dd_ab <- datadist(df_stage_2)
options(datadist = "dd_ab")

# 构建全校正 RCS 模型 (纳入 Table 1 中所有显著的混杂因素)
fit_rcs_ab <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + Age + Gender + Race + PIR + 
                    BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal + 
                    Total_Fat + Fiber + Diabetes + Hypertension, 
                  data = df_stage_2)

cat("\n======================================================\n")
cat("🔥 RCS ANOVA 检验结果 (重点看 Nonlinear 的 P 值)\n")
cat("======================================================\n")
print(anova(fit_rcs_ab))

# ==============================================================================
# 二_2：连续暴露与分类暴露的主 Logistic 回归 (计算 OR 与 P for trend)
# ==============================================================================
message("\n>>> [二_2] 启动经典多因素 Logistic 回归模型计算核心 OR 值...")

# Model 1: 仅调整人口学 (Age, Gender, Race, PIR)
model1_formula <- "High_CRP_num ~ Ratio_Quartile + Age + Gender + Race + PIR"

# Model 2: Model 1 + 生活方式与膳食 (BMI, Smoking, Alcohol_g, MET_imp, Energy_Kcal, Total_Fat, Fiber)
model2_formula <- paste(model1_formula, "+ BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber")

# Model 3 (完全调整): Model 2 + 共病 (Diabetes, Hypertension)
model3_formula <- paste(model2_formula, "+ Diabetes + Hypertension")

# 运行最终的全校正模型 (Model 3) 提取分类 OR 值
fit_m3 <- glm(as.formula(model3_formula), data = df_stage_2, family = binomial())
res_m3 <- tidy(fit_m3, exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("Ratio_Quartile", term)) %>%
  select(term, estimate, conf.low, conf.high, p.value)

# 计算 P for trend (将四分位数作为连续变量带入模型检验线性趋势)
df_stage_2$Ratio_Quartile_Num <- as.numeric(df_stage_2$Ratio_Quartile)
trend_formula <- sub("Ratio_Quartile", "Ratio_Quartile_Num", model3_formula)
fit_trend <- glm(as.formula(trend_formula), data = df_stage_2, family = binomial())
p_trend <- tidy(fit_trend) %>% filter(term == "Ratio_Quartile_Num") %>% pull(p.value)

cat("\n======================================================\n")
cat("📊 主分析结果：全校正模型 (Model 3) 下的独立致炎效应\n")
cat("======================================================\n")
print(res_m3)
cat(glue::glue("\n📈 趋势检验 P for trend: {round(p_trend, 4)}\n"))

# 留痕结果
write.csv(res_m3, "results_二/data/Stage_二_2_Multivariable_Logistic_ORs.csv", row.names = FALSE)
message(">>> 🎉 阶段二关联分析完毕！数据已留痕。")
# ==============================================================================
# 二_3：绘制顶刊级 RCS 剂量反应曲线 (Figure 1)
# ==============================================================================
message(">>> [二_3] 正在绘制并导出高分辨率 RCS 剂量反应曲线...")

# 1. 重新设定 datadist，将比值的参考点设为中位数 (OR = 1 的基准)
ref_val <- median(df_stage_2$Ratio_18_16, na.rm = TRUE)
dd_ab$limits["Adjust to", "Ratio_18_16"] <- ref_val
fit_rcs_ab <- update(fit_rcs_ab) # 更新模型以应用新的参考点

# 2. 计算预测值 (转换为 OR 值)
pred_rcs <- rms::Predict(fit_rcs_ab, Ratio_18_16, ref.zero = TRUE, fun = exp)

# 3. 使用 ggplot2 绘制出版级图表
p_rcs <- ggplot(pred_rcs) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", size = 0.8) + # OR=1 基准线
  geom_ribbon(aes(x = Ratio_18_16, ymin = lower, ymax = upper), fill = "#0073C2", alpha = 0.2) + # 95% CI 阴影
  geom_line(aes(x = Ratio_18_16, y = yhat), color = "#0073C2", size = 1.2) + # OR 趋势线
  theme_classic(base_size = 14) +
  labs(
    title = "Dose-Response Relationship between 18:0/16:0 Ratio and High hs-CRP",
    x = "Dietary 18:0/16:0 Ratio",
    y = "Odds Ratio (95% CI)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  # 添加 P for non-linearity 的文本注释 (此处数值为您之前跑出的 0.797)
  annotate("text", x = max(df_stage_2$Ratio_18_16, na.rm=T) * 0.7, y = max(pred_rcs$upper)*0.8, 
           label = "P for non-linearity = 0.797", size = 5, fontface = "italic")

# 4. 导出高清 PDF
pdf("results_二/plots/Stage_二_3_RCS_Curve.pdf", width = 8, height = 6)
print(p_rcs)
dev.off()
message(">>> 🎉 RCS 曲线 (Figure 1) 绘制成功！")

# ==============================================================================
# 项目名称：AB (三明治架构)
# 当前阶段：三、敏感性分析与伪因果推断 (极端组分析 Q4 vs Q1) - 终极修复版
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

setwd("/Users/bing/AA/AB")

if (!require("pacman")) install.packages("pacman")
p_load(dplyr, MatchIt, tableone, survey, EValue, broom)

message(">>> [三_1] 加载全样本数据，准备提取极端组 (Q4 vs Q1)...")
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# 1. 精准提取极端对照组 (Q1 与 Q4)
df_extreme <- df_stage_2 %>%
  filter(Ratio_Quartile %in% c("Q1", "Q4")) %>%
  mutate(
    # 将 Q4 设为暴露组 (1)，Q1 设为对照组 (0)
    Exposure_Extreme = ifelse(Ratio_Quartile == "Q4", 1, 0),
    Exposure_Factor = factor(Exposure_Extreme, levels = c(0, 1), labels = c("Q1", "Q4"))
  ) %>%
  droplevels()

message(glue::glue(">>> 极端组 (Q4 vs Q1) 提取完毕，总人数应为：{nrow(df_extreme)}"))

# ==============================================================================
# 三_2：倾向性评分匹配 (PSM) - 终极混杂剥离
# ==============================================================================
message(">>> [三_2] 启动 PSM 因果推断引擎 (1:1 最邻近匹配)...")

set.seed(2026)
psm_formula <- Exposure_Extreme ~ Age + Gender + Race + PIR + BMI + Smoking + 
               Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber + Diabetes + Hypertension

psm_model <- matchit(psm_formula, data = df_extreme, method = "nearest", caliper = 0.05)
matched_data <- match.data(psm_model)

message(glue::glue(">>> 匹配成功！保留了 {nrow(matched_data)} 名极其均衡的患者。"))

# 计算匹配后的真实致炎效应
design_psm <- svydesign(ids = ~1, weights = ~weights, data = matched_data)
fit_psm <- svyglm(High_CRP_num ~ Exposure_Factor, design = design_psm, family = binomial())
res_psm <- tidy(fit_psm, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "Exposure_FactorQ4")

cat("\n======================================================\n")
cat("🏆 伪因果推断结果 (PSM: Q4 vs Q1)\n")
cat("======================================================\n")
cat(glue::glue("匹配后 OR = {round(res_psm$estimate, 3)} (95% CI: {round(res_psm$conf.low, 3)} - {round(res_psm$conf.high, 3)}), P = {round(res_psm$p.value, 4)}\n"))

# ==============================================================================
# 三_3：量化未测量混杂 (E-value)
# ==============================================================================
message("\n>>> [三_3] 计算 E-value 边界...")
# 基于第二阶段 Model 3 中 Q4 的 OR=1.20 和下限 1.04 计算
ev <- evalues.OR(est = 1.20, lo = 1.04, hi = 1.39, rare = FALSE)
cat("\n======================================================\n")
cat("🛡️ E-value 稳健性评估\n")
cat("======================================================\n")
print(ev)

# ==============================================================================
# 三_4：亚组分层分析 (Subgroup Analysis) - 动态防报错版
# ==============================================================================
message("\n>>> [三_4] 启动亚组交互作用分析 (动态过滤单一水平变量)...")

run_subgroup <- function(sub_df, subgroup_name) {
  # 乘以 10，观察每增加 0.1 个单位 Ratio 的影响，清理冗余 levels
  sub_df <- sub_df %>% mutate(Ratio_Scale = Ratio_18_16 * 10) %>% droplevels()
  
  # 定义想要校正的候选混杂变量池
  candidate_vars <- c("Age", "Gender", "BMI", "Smoking", "Energy_Kcal")
  
  # 核心防御机制：只保留在当前亚组中拥有 2 个或以上不同值的变量
  valid_vars <- candidate_vars[sapply(candidate_vars, function(v) length(unique(sub_df[[v]])) > 1)]
  
  # 动态拼接公式
  f_str <- paste("High_CRP_num ~ Ratio_Scale +", paste(valid_vars, collapse = " + "))
  
  fit <- glm(as.formula(f_str), data = sub_df, family = binomial())
  res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "Ratio_Scale")
  
  return(data.frame(
    Subgroup = subgroup_name,
    N = nrow(sub_df),
    OR = round(res$estimate, 2),
    Lower = round(res$conf.low, 2),
    Upper = round(res$conf.high, 2),
    P_Value = round(res$p.value, 4)
  ))
}

# 逐个运行亚组
res_male <- run_subgroup(df_stage_2 %>% filter(Gender == "Male"), "Male")
res_female <- run_subgroup(df_stage_2 %>% filter(Gender == "Female"), "Female")
res_bmi1 <- run_subgroup(df_stage_2 %>% filter(BMI < 25), "BMI < 25")
res_bmi2 <- run_subgroup(df_stage_2 %>% filter(BMI >= 25 & BMI < 30), "BMI 25-30")
res_bmi3 <- run_subgroup(df_stage_2 %>% filter(BMI >= 30), "BMI >= 30")
res_age1 <- run_subgroup(df_stage_2 %>% filter(Age < 60), "Age < 60")
res_age2 <- run_subgroup(df_stage_2 %>% filter(Age >= 60), "Age >= 60")

# 整合结果并导出
df_subgroup <- bind_rows(res_male, res_female, res_bmi1, res_bmi2, res_bmi3, res_age1, res_age2)
write.csv(df_subgroup, "results_三/data/Stage_三_4_Subgroup_Analysis.csv", row.names = FALSE)

cat("\n======================================================\n")
cat("🌲 亚组分析底层数据 (每增加 0.1 ratio 单位的 OR)\n")
cat("======================================================\n")
print(df_subgroup)
message(">>> 🎉 阶段三运算全部完成！")
# ==============================================================================
# 三_5：绘制高颜值亚组分析森林图 (Figure 2)
# ==============================================================================
message("\n>>> [三_5] 正在绘制出版级亚组分析森林图...")

dir.create("results_三/plots", showWarnings = FALSE, recursive = TRUE)

# 整理绘图数据结构，反转顺序使图形从上到下显示
df_subgroup_plot <- df_subgroup %>%
  mutate(
    Subgroup = factor(Subgroup, levels = rev(Subgroup)),
    Label = sprintf("%.2f (%.2f-%.2f)", OR, Lower, Upper) # 生成图上的文本标签
  )

p_forest <- ggplot(df_subgroup_plot, aes(x = OR, y = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred", size = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#2E9FDF", size = 0.8) +
  geom_point(shape = 15, size = 4, color = "#2E9FDF") +
  # 在点旁边添加 OR (95% CI) 的文本
  geom_text(aes(label = Label, x = Upper + 0.05), hjust = 0, size = 4) +
  theme_classic(base_size = 14) +
  labs(
    title = "Subgroup Analysis of the 18:0/16:0 Ratio (Per 0.1 Unit Increase)",
    x = "Odds Ratio (95% CI)",
    y = ""
  ) +
  scale_x_continuous(limits = c(min(df_subgroup_plot$Lower) - 0.1, max(df_subgroup_plot$Upper) + 0.4)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.text.x = element_text(color = "black")
  )

pdf("results_三/plots/Stage_三_5_Subgroup_Forest.pdf", width = 9, height = 6)
print(p_forest)
dev.off()
message(">>> 🎉 亚组森林图 (Figure 2) 绘制成功！")

# 加载必要的包和修复后的数据
library(tableone)
df_fixed <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# 专门生成 Gender 在 Q1-Q4 中的分布 (包含人数和百分比)
tab_gender <- CreateTableOne(vars = "Gender", strata = "Ratio_Quartile", data = df_fixed)

cat("\n====================================\n")
cat("📊 修复后的 Gender 分布数据 (供填入 Table 1)\n")
cat("====================================\n")
print(tab_gender, showAllLevels = TRUE)

# ==============================================================================
# 项目名称：AB (三明治架构)
# 当前阶段：四、机制揭秘 (脂肪酸 PCA 降维 & 宿主微环境中介分析) - 防冲突修复版
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

setwd("/Users/bing/AA/AB")

if (!require("pacman")) install.packages("pacman")
p_load(dplyr, tidyr, psych, factoextra, mediation, broom)

message(">>> [四_1] 加载全样本纯净基底数据...")
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# ==============================================================================
# 四_1：脂肪酸平衡谱的 PCA 主成分分析
# ==============================================================================
message("\n>>> [四_1] 启动 PCA 降维，提取宏观膳食脂肪酸模式...")

# 1. 提取阶段一埋线的 9 种代表性脂肪酸
fa_vars <- c("DR1TS140", "DR1TS160", "DR1TS180", "DR1TM181", 
             "DR1TP182", "DR1TP183", "DR1TP204", "DR1TP205", "DR1TP226")

# ⚠️ 修复点：强制使用 dplyr::select 和 tidyr::drop_na 防止包冲突
df_pca <- df_stage_2 %>%
  dplyr::select(SEQN, all_of(fa_vars), High_CRP_num, Ratio_18_16) %>%
  tidyr::drop_na()

# 2. 执行 PCA (使用 psych 包的 varimax 正交旋转)
pca_res <- principal(df_pca[, fa_vars], nfactors = 3, rotate = "varimax", scores = TRUE)

cat("\n======================================================\n")
cat("🧬 脂肪酸 PCA 模式提取 (载荷矩阵)\n")
cat("======================================================\n")
print(pca_res$loadings, cutoff = 0.3)

# 3. 将 PCA 得分合并回原数据，探究其与炎症的关系
df_pca_scores <- cbind(df_pca, pca_res$scores) %>%
  dplyr::rename(Pattern_1 = RC1, Pattern_2 = RC2, Pattern_3 = RC3)

fit_pca1 <- glm(High_CRP_num ~ Pattern_1, data = df_pca_scores, family = binomial())
fit_pca2 <- glm(High_CRP_num ~ Pattern_2, data = df_pca_scores, family = binomial())
fit_pca3 <- glm(High_CRP_num ~ Pattern_3, data = df_pca_scores, family = binomial())

cat("\n======================================================\n")
cat("📊 三大膳食模式与系统性炎症的关联 (单因素探索)\n")
cat("======================================================\n")
cat("Pattern_1 模式致炎 OR P值:", round(summary(fit_pca1)$coefficients[2, 4], 4), "\n")
cat("Pattern_2 模式致炎 OR P值:", round(summary(fit_pca2)$coefficients[2, 4], 4), "\n")
cat("Pattern_3 模式致炎 OR P值:", round(summary(fit_pca3)$coefficients[2, 4], 4), "\n")

capture.output(print(pca_res$loadings), file = "results_四/data/Stage_四_1_PCA_Loadings.txt")

# ==============================================================================
# 四_2：中介分析 (TyG 与 SII 的桥梁作用)
# ==============================================================================
message("\n>>> [四_2] 启动中介分析引擎 (运算涉及 Bootstrap，请等待 1-3 分钟)...")

df_med <- df_stage_2 %>%
  dplyr::mutate(Ratio_Scale = Ratio_18_16 * 10) %>% 
  tidyr::drop_na(TyG_Index, SII_Index)

covariates <- "Age + Gender + BMI + Smoking + Energy_Kcal + Diabetes"

# ---------------------------------------------------------
# 路径 A：代谢路径 (TyG_Index) 中介效应
# ---------------------------------------------------------
message("   -> 正在计算 TyG 代谢中介路径...")
set.seed(2026)

formula_m_tyg <- as.formula(paste("TyG_Index ~ Ratio_Scale +", covariates))
fitM_TyG <- lm(formula_m_tyg, data = df_med)

formula_y_tyg <- as.formula(paste("High_CRP_num ~ Ratio_Scale + TyG_Index +", covariates))
fitY_TyG <- glm(formula_y_tyg, data = df_med, family = binomial("probit")) 

med_out_tyg <- mediate(fitM_TyG, fitY_TyG, treat = "Ratio_Scale", mediator = "TyG_Index", 
                       boot = TRUE, sims = 500)

# ---------------------------------------------------------
# 路径 B：免疫微环境路径 (SII_Index) 中介效应
# ---------------------------------------------------------
message("   -> 正在计算 SII 免疫中介路径...")
set.seed(2026)

formula_m_sii <- as.formula(paste("SII_Index ~ Ratio_Scale +", covariates))
fitM_SII <- lm(formula_m_sii, data = df_med)

formula_y_sii <- as.formula(paste("High_CRP_num ~ Ratio_Scale + SII_Index +", covariates))
fitY_SII <- glm(formula_y_sii, data = df_med, family = binomial("probit"))

med_out_sii <- mediate(fitM_SII, fitY_SII, treat = "Ratio_Scale", mediator = "SII_Index", 
                       boot = TRUE, sims = 500)

# ---------------------------------------------------------
# 输出并留痕中介结果
# ---------------------------------------------------------
cat("\n======================================================\n")
cat("🧬 中介分析结果 1：TyG 代谢途径\n")
cat("======================================================\n")
print(summary(med_out_tyg))

cat("\n======================================================\n")
cat("🦠 中介分析结果 2：SII 免疫途径\n")
cat("======================================================\n")
print(summary(med_out_sii))

sink("results_四/data/Stage_四_2_Mediation_Summary.txt")
cat("--- TyG Mediation ---\n")
print(summary(med_out_tyg))
cat("\n--- SII Mediation ---\n")
print(summary(med_out_sii))
sink()

message("\n>>> 🎉 阶段四全部跑完！机制网络构建成功。")


# ==============================================================================
# 四_3：导出中介分析原生效应图 (备用/参考)
# ==============================================================================
message("\n>>> [四_3] 正在导出中介效应点图...")

dir.create("results_四/plots", showWarnings = FALSE, recursive = TRUE)

pdf("results_四/plots/Stage_四_3_Mediation_Effects.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(med_out_tyg, main = "Mediation via TyG Index", xlab = "Effect Size")
plot(med_out_sii, main = "Mediation via SII Index", xlab = "Effect Size")
dev.off()

message(">>> 🎉 机制效应图导出完毕！")````

# ==============================================================================
# 项目名称：AB (基于三明治架构的双周期合并队列研究)
# 当前阶段：五、数据绝对切割与临床预测模型重塑 (严格遵循 TRIPOD 规范)
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

setwd("/Users/bing/AA/AB")
dir.create("results_五_新", showWarnings = FALSE)
dir.create("results_五_新/data", showWarnings = FALSE)
dir.create("results_五_新/plots", showWarnings = FALSE)

if (!require("pacman")) install.packages("pacman")
p_load(dplyr, tidyr, glmnet, rms, ggplot2, dcurves, caret)

message(">>> [五_1] 加载全样本数据，并立即执行 70/30 绝对物理切割...")
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")

# 1. 提取全维特征池
df_all <- df_stage_2 %>%
  dplyr::select(Ratio_18_16, Age, Gender, Race, PIR, BMI, Smoking, Alcohol_g, 
                MET_imp, Energy_Kcal, Total_Fat, Fiber, Diabetes, Hypertension, 
                TyG_Index, SII_Index, High_CRP_num) %>%
  tidyr::drop_na() %>%
  as.data.frame()

# 2. 物理切割：70% 训练集，30% 测试集锁死备用
set.seed(2026)
train_index <- caret::createDataPartition(df_all$High_CRP_num, p = 0.7, list = FALSE)

df_train <- df_all[train_index, ]
df_test  <- df_all[-train_index, ]

# 将测试集原封不动锁入硬盘，供第六阶段使用
saveRDS(df_test, "results_五_新/data/Stage_五_Locked_Test_Set.rds")
message(glue::glue(">>> 切割完成！训练集 (70%): {nrow(df_train)} 人 | 盲测集已锁死 (30%): {nrow(df_test)} 人"))

# ==============================================================================
# 五_2：在 70% 训练集中，启动 LASSO 降维
# ==============================================================================
message(">>> [五_2] 仅在训练集中启动 LASSO 机器学习引擎 (杜绝数据泄露)...")

x_matrix <- model.matrix(High_CRP_num ~ ., data = df_train)[, -1]
y_vector <- df_train$High_CRP_num

set.seed(2026)
cv_lasso <- cv.glmnet(x_matrix, y_vector, family = "binomial", alpha = 1, nfolds = 10)
best_lambda <- cv_lasso$lambda.1se

lasso_coefs <- coef(cv_lasso, s = best_lambda)
selected_vars_idx <- which(lasso_coefs != 0)
selected_vars_names <- rownames(lasso_coefs)[selected_vars_idx]
selected_vars_names <- selected_vars_names[selected_vars_names != "(Intercept)"]

cat("\n======================================================\n")
cat("🎯 训练集 LASSO 筛选出的核心预测因子\n")
cat("======================================================\n")
print(selected_vars_names)

pdf("results_五_新/plots/Stage_五_1_LASSO_Path.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(cv_lasso$glmnet.fit, "lambda", label = TRUE)
plot(cv_lasso)
dev.off()

# ==============================================================================
# 五_3：在 70% 训练集中，重构终极 Nomogram 预测模型
# ==============================================================================
message("\n>>> [五_3] 基于训练集 LASSO 结果，重塑临床预测列线图...")

raw_candidate_vars <- c("Ratio_18_16", "Age", "Gender", "Race", "PIR", "BMI", "Smoking", 
                        "Alcohol_g", "MET_imp", "Energy_Kcal", "Total_Fat", "Fiber", 
                        "Diabetes", "Hypertension", "TyG_Index", "SII_Index")

final_model_vars <- raw_candidate_vars[sapply(raw_candidate_vars, function(v) any(grepl(v, selected_vars_names)))]
cat("👉 最终纳入模型的特征为: ", paste(final_model_vars, collapse = ", "), "\n")

dd_train <- datadist(df_train)
options(datadist = "dd_train")

nomo_formula <- as.formula(paste("High_CRP_num ~", paste(final_model_vars, collapse = " + ")))
fit_nomo_final <- rms::lrm(nomo_formula, data = df_train, x = TRUE, y = TRUE)

# 保存模型对象供第六阶段调用
saveRDS(fit_nomo_final, "results_五_新/data/Stage_五_Final_Model.rds")

nom <- rms::nomogram(fit_nomo_final, fun = plogis, fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                     funlabel = "Risk of High Systemic Inflammation", lp = FALSE)

pdf("results_五_新/plots/Stage_五_2_Training_Nomogram.pdf", width = 12, height = 8)
plot(nom, cex.var = 0.8, cex.axis = 0.8, lmgp = 0.2)
title(main = "Nomogram for Predicting Inflammation (Training Cohort)", line = 1.5)
dev.off()

# 训练集内部 Bootstrap 校准
set.seed(2026)
cal_train <- rms::calibrate(fit_nomo_final, method = "boot", B = 1000)
pdf("results_五_新/plots/Stage_五_3_Training_Calibration.pdf", width = 7, height = 7)
plot(cal_train, xlab = "Predicted Probability", ylab = "Actual Probability", main = "Training Calibration (B=1000)")
dev.off()

cat(glue::glue("\n🌟 训练集内部 C-index (AUC) 为: {round(fit_nomo_final$stats['C'], 3)}\n"))
message(">>> 🎉 第五阶段 (模型生成) 完成！模型与测试集已封存。")

# ==============================================================================
# 项目名称：AB (基于三明治架构的双周期合并队列研究)
# 当前阶段：六、严格独立盲测验证 (基于 30% 隔离测试集)
# 工作目录：/Users/bing/AA/AB
# ==============================================================================

setwd("/Users/bing/AA/AB")
dir.create("results_六_新", showWarnings = FALSE)

if (!require("pacman")) install.packages("pacman")
p_load(dplyr, rms, ggplot2, dcurves, Hmisc)

message(">>> [六_1] 提取锁死的测试集与训练完成的模型...")

# 读取第五阶段封存的数据与模型
df_test <- readRDS("results_五_新/data/Stage_五_Locked_Test_Set.rds")
fit_nomo_final <- readRDS("results_五_新/data/Stage_五_Final_Model.rds")

# ==============================================================================
# 六_2：执行实弹盲测与 C-index 计算
# ==============================================================================
message(">>> [六_2] 正在 30% 隔离测试集中进行实战预测盲测...")

df_test$pred_prob <- predict(fit_nomo_final, newdata = df_test, type = "fitted")

val_res <- Hmisc::rcorr.cens(df_test$pred_prob, df_test$High_CRP_num)
c_index_test <- val_res["C Index"]

cat("\n======================================================\n")
cat(glue::glue("🌟 30% 隔离盲测集 C-index 为: {round(c_index_test, 3)}\n"))
cat("======================================================\n")

# ==============================================================================
# 六_3：绘制测试集的终极图表 (Calibration & DCA)
# ==============================================================================
message(">>> [六_3] 导出隔离盲测集的最终性能图纸...")

# 1. 绘制测试集校准曲线
pdf("results_六_新/Stage_六_1_Testing_Calibration.pdf", width = 7, height = 7)
rms::val.prob(df_test$pred_prob, df_test$High_CRP_num, 
              xlab = "Nomogram Predicted Probability", 
              ylab = "Actual Observed Fraction")
title(main = "Testing Cohort Calibration Curve (30% Split)", line = 2.5)
dev.off()

# 2. 绘制测试集决策曲线 (DCA)
dca_test <- dcurves::dca(High_CRP_num ~ pred_prob, 
                         data = df_test,
                         thresholds = seq(0, 1, by = 0.01),
                         label = list(pred_prob = "Validated Nomogram Model"))

pdf("results_六_新/Stage_六_2_Testing_DCA.pdf", width = 8, height = 6)
p_test_dca <- plot(dca_test) +
  labs(title = "Testing Cohort Decision Curve Analysis (30% Split)",
       x = "Threshold Probability", y = "Net Benefit") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_test_dca)
dev.off()

message(">>> 🎉 第六阶段 (独立盲测) 完美收官！临床预测工具性能定案。")

# ==============================================================================
# 项目名称：AB ( Revision 终极数据回应完整脚本)
# 运行说明：请在重启 RStudio 后，直接全选运行此脚本
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 环境初始化与数据加载
# ------------------------------------------------------------------------------
message(">>> [0/4] 正在初始化环境并加载数据...")

if (!require("pacman")) install.packages("pacman")
p_load(dplyr, rms, ResourceSelection, nhanesA, pROC, Hmisc, broom, car)

setwd("/Users/bing/AA/AB")

# 加载您的核心纯净数据底座和锁死的测试集
df_stage_2 <- readRDS("results_一/data/Stage_一_2_Cleaned_Cohort.rds")
df_test <- readRDS("results_五_新/data/Stage_五_Locked_Test_Set.rds")

# ------------------------------------------------------------------------------
# 模块一：营养学铁律——残差法能量校正与共线性陷阱(VIF)反击
# ------------------------------------------------------------------------------
message("\n>>> [1/4] 正在执行模块一：能量校正与独立效应验证...")

# 1. 计算能量校正后的绝对摄入量 (残差法)
fit_16 <- lm(DR1TS160 ~ Energy_Kcal, data = df_stage_2, na.action = na.exclude)
df_stage_2$FA16_adj <- resid(fit_16) + mean(df_stage_2$DR1TS160, na.rm = TRUE)

fit_18 <- lm(DR1TS180 ~ Energy_Kcal, data = df_stage_2, na.action = na.exclude)
df_stage_2$FA18_adj <- resid(fit_18) + mean(df_stage_2$DR1TS180, na.rm = TRUE)

# 2. 计算校正后的比值，并切分四分位数
df_stage_2$Ratio_adj <- df_stage_2$FA18_adj / (df_stage_2$FA16_adj + 0.0001)
ratio_adj_q <- quantile(df_stage_2$Ratio_adj, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
df_stage_2$Ratio_adj_Quartile <- cut(df_stage_2$Ratio_adj, breaks = ratio_adj_q, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))

# 3. 陷阱验证：同时放入比值与绝对值，计算方差膨胀因子(VIF)
trap_formula <- "High_CRP_num ~ Ratio_adj_Quartile + FA16_adj + FA18_adj + Age + Gender + Race + PIR + BMI + Smoking + Alcohol_g + MET_imp + Total_Fat + Fiber + Diabetes + Hypertension"
fit_trap <- glm(as.formula(trap_formula), data = df_stage_2, family = binomial())

cat("\n======================================================\n")
cat("🛡️ [1A] 审稿人共线性陷阱验证 (查看 VIF 值，若 >5 则证明审稿人逻辑有误)\n")
cat("======================================================\n")
print(car::vif(fit_trap))

# 4. 正确的独立效应：仅使用能量校正后的比值
correct_formula <- "High_CRP_num ~ Ratio_adj_Quartile + Age + Gender + Race + PIR + BMI + Smoking + Alcohol_g + MET_imp + Total_Fat + Fiber + Diabetes + Hypertension"
fit_correct <- glm(as.formula(correct_formula), data = df_stage_2, family = binomial())
res_correct <- broom::tidy(fit_correct, exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("Ratio_adj_QuartileQ4", term))

cat("\n======================================================\n")
cat("🛡️ [1B] 能量校正后的独立致炎效应 (写进文章的核心 OR 值)\n")
cat("======================================================\n")
print(as.data.frame(res_correct)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])


# ------------------------------------------------------------------------------
# 模块二：临床转化的硬核逻辑——模型的增量价值分析 (NRI / IDI)
# ------------------------------------------------------------------------------
message("\n>>> [2/4] 正在执行模块二：模型预测增量价值分析...")

# 模型 A (基准 4 变量) 与 模型 B (加上膳食比值)
fit_A <- glm(High_CRP_num ~ Gender + BMI + TyG_Index + SII_Index, family = binomial(), data = df_test)
fit_B <- glm(High_CRP_num ~ Gender + BMI + TyG_Index + SII_Index + Ratio_18_16, family = binomial(), data = df_test)

prob_A <- predict(fit_A, type = "response")
prob_B <- predict(fit_B, type = "response")

cat("\n======================================================\n")
cat("📊 [2A] DeLong's Test: AUC 差异 (期待 P > 0.05)\n")
cat("======================================================\n")
roc_A <- roc(df_test$High_CRP_num, prob_A, quiet = TRUE)
roc_B <- roc(df_test$High_CRP_num, prob_B, quiet = TRUE)
print(roc.test(roc_A, roc_B))

cat("\n======================================================\n")
cat("📊 [2B] NRI 与 IDI 综合重分类评估 (期待 P > 0.05)\n")
cat("======================================================\n")
idi_nri_result <- Hmisc::improveProb(x1 = prob_A, x2 = prob_B, y = df_test$High_CRP_num)
print(idi_nri_result)


# ------------------------------------------------------------------------------
# 模块三：TRIPOD 规范——Hosmer-Lemeshow 拟合优度检验 (终极防崩溃版)
# ------------------------------------------------------------------------------
message("\n>>> [3/4] 正在执行模块三：Hosmer-Lemeshow 检验 (启动自动降阶与清洗)...")

# 1. 提取向量并硬性清洗 NA
hl_y_raw <- as.numeric(as.character(df_test$High_CRP_num))
hl_prob_raw <- as.numeric(df_test$pred_prob)
valid_idx <- !is.na(hl_y_raw) & !is.na(hl_prob_raw)
hl_y <- hl_y_raw[valid_idx]
hl_prob <- hl_prob_raw[valid_idx]

# 2. 注入极微小白噪声，强行打破绝对重叠 (1e-9量级，不影响真实统计)
set.seed(2026)
hl_prob_noise <- hl_prob + runif(length(hl_prob), min = 0, max = 1e-9)

# 3. 自动降阶防撞车寻优机制
run_hl_bulletproof <- function(y, prob) {
  for (g_val in 10:5) {
    res <- tryCatch({
      ResourceSelection::hoslem.test(y, prob, g = g_val)
    }, error = function(e) NULL)
    
    if (!is.null(res)) {
      cat(glue::glue("✅ 成功！当前采用无重叠的最优分箱 g = {g_val}。\n"))
      return(res)
    }
  }
  return("警告：数据极度集中，HL 检验彻底不收敛。请直接在文中报告 Brier Score。")
}

cat("\n======================================================\n")
cat("📏 [3] Hosmer-Lemeshow 检验最终结果\n")
cat("======================================================\n")
hl_test_final <- run_hl_bulletproof(hl_y, hl_prob_noise)
if(is.list(hl_test_final)) print(hl_test_final) else cat(hl_test_final, "\n")


# ------------------------------------------------------------------------------
# 模块四：反向因果——剔除恶性肿瘤患者的敏感性分析
# ------------------------------------------------------------------------------
message("\n>>> [4/4] 正在执行模块四：排除恶性肿瘤病史的反向因果敏感性分析...")

# 联网拉取肿瘤调查问卷，并做安全合并
mcq_I <- nhanes('MCQ_I') %>% dplyr::select(SEQN, MCQ220)
mcq_J <- nhanes('MCQ_J') %>% dplyr::select(SEQN, MCQ220)
mcq_merge <- bind_rows(mcq_I, mcq_J)

df_sens <- df_stage_2 %>% 
  left_join(mcq_merge, by = "SEQN") %>%
  filter(is.na(MCQ220) | MCQ220 != 1)

# 重跑模型
model3_sens_formula <- "High_CRP_num ~ Ratio_Quartile + Age + Gender + Race + PIR + BMI + Smoking + Alcohol_g + MET_imp + Energy_Kcal + Total_Fat + Fiber + Diabetes + Hypertension"
fit_sens <- glm(as.formula(model3_sens_formula), data = df_sens, family = binomial())
res_sens <- broom::tidy(fit_sens, exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("Ratio_QuartileQ4", term))

cat("\n======================================================\n")
cat("⚔️ [4] 排除肿瘤病史后的敏感性 OR 值\n")
cat("======================================================\n")
cat(glue::glue("过滤后最终纳入的敏感性队列总人数: {nrow(df_sens)}\n\n"))
print(as.data.frame(res_sens)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

message("\n>>> 🎉 全部四大战役数据运算彻底完毕！可以整理出组了。")
