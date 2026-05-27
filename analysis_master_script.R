# =============================================================================
# Master Script: 硬脂酸/棕榈酸比值与系统性炎症网络 (小而美路线)
# 包含：数据拉取 -> 特征工程 -> ML 降维 -> SHAP 解释 -> 自动留痕
# 工作目录：/Users/bing/AA
# =============================================================================


# 初始化工作目录与结果文件夹
setwd("/Users/bing/AA")
dir.create("results", showWarnings = FALSE)
dir.create("results/data", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

# =============================================================================
# =============================================================================
# 第一阶段：提取底层核心指标 (涵盖饮食、生化、血常规) + 新增混杂因素
# =============================================================================
if (!require("pacman")) install.packages("pacman")
p_load(nhanesA, dplyr, tidyr)

message(">>> 🚀 [Stage 1] 开始从 CDC 在线拉取底层原始数据与协变量...")

# 1. 拉取核心表
demo <- nhanes("DEMO_J")
diet <- nhanes("DR1TOT_J")
cbc  <- nhanes("CBC_J")
glu  <- nhanes("GLU_J")
tg   <- nhanes("TRIGLY_J")
crp  <- nhanes("HSCRP_J")  

# 2. 拉取第三方要求的新增敏感性分析协变量表
bmx  <- nhanes("BMX_J")    # 体格检查 (用于 BMI)
smq  <- nhanes("SMQ_J")    # 吸烟问卷
diq  <- nhanes("DIQ_J")    # 糖尿病问卷
paq  <- nhanes("PAQ_J")    # 体力活动问卷

# 3. 提取变量子集
demo_sub <- demo %>% select(SEQN, RIDAGEYR, RIAGENDR)
diet_sub <- diet %>% select(SEQN, DR1TS180, DR1TS160, DR1TKCAL, DR1TALCO)
cbc_sub  <- cbc  %>% select(SEQN, LBXPLTSI, LBDNENO, LBDLYMNO)
glu_sub  <- glu  %>% select(SEQN, LBXGLU)
tg_sub   <- tg   %>% select(SEQN, LBXTR)
crp_sub  <- crp  %>% select(SEQN, LBXHSCRP)
bmx_sub  <- bmx  %>% select(SEQN, BMXBMI)
smq_sub  <- smq  %>% select(SEQN, SMQ020, SMQ040)
diq_sub  <- diq  %>% select(SEQN, DIQ010)
paq_sub  <- paq  %>% select(SEQN, any_of(c("PAQ610", "PAD615", "PAQ625", "PAD630", "PAQ640", "PAD645", "PAQ655", "PAD660", "PAQ670", "PAD675")))

# 4. 合并数据
df_raw <- demo_sub %>%
  left_join(diet_sub, by = "SEQN") %>%
  left_join(cbc_sub, by = "SEQN") %>%
  left_join(glu_sub, by = "SEQN") %>%
  left_join(tg_sub, by = "SEQN") %>%
  left_join(crp_sub, by = "SEQN") %>%
  left_join(bmx_sub, by = "SEQN") %>%
  left_join(smq_sub, by = "SEQN") %>%
  left_join(diq_sub, by = "SEQN") %>%
  left_join(paq_sub, by = "SEQN")      

# 5. 核心清洗与衍生变量计算 (已融合 Smoking 和 Diabetes 的类型修复补丁！)
df_clean <- df_raw %>%
  # 🚨 极其关键的防御：只剔除“核心变量”存在 NA 的行，确保 1973 这个数字绝对不变！
  drop_na(RIDAGEYR, RIAGENDR, DR1TS180, DR1TS160, LBXPLTSI, LBDNENO, LBDLYMNO, LBXGLU, LBXTR, LBXHSCRP) %>%
  filter(RIDAGEYR >= 20) %>%
  filter(LBXHSCRP <= 10) %>% # 剔除急性感染人群
  mutate(
    # --- 新增的 6 大混杂因素 (已完美修复因子编码 Bug) ---
    BMI = BMXBMI,
    Energy_Kcal = DR1TKCAL,
    Alcohol_g = ifelse(is.na(DR1TALCO), 0, DR1TALCO),
    Diabetes = case_when(
      as.character(DIQ010) %in% c("1", "Yes") ~ "Yes",
      as.character(DIQ010) %in% c("2", "3", "No", "Borderline") ~ "No/Borderline",
      TRUE ~ "Unknown"
    ),
    Smoking = case_when(
      as.character(SMQ020) %in% c("2", "No") ~ "Never",
      as.character(SMQ020) %in% c("1", "Yes") & as.character(SMQ040) %in% c("3", "Not at all") ~ "Former",
      as.character(SMQ020) %in% c("1", "Yes") & as.character(SMQ040) %in% c("1", "2", "Every day", "Some days") ~ "Current",
      TRUE ~ "Unknown"
    ),
    
    # 极其严谨的 NHANES 2017-2018 MET-min/week 计算 (完美处理 Skip Pattern 陷阱，将合法 NA 转换为 0)
    MET_min_week = 
      (replace_na(ifelse(PAQ610 %in% 1:7, PAQ610, 0), 0) * replace_na(ifelse(PAD615 < 7777, PAD615, 0), 0) * 8) +  
      (replace_na(ifelse(PAQ625 %in% 1:7, PAQ625, 0), 0) * replace_na(ifelse(PAD630 < 7777, PAD630, 0), 0) * 4) +  
      (replace_na(ifelse(PAQ640 %in% 1:7, PAQ640, 0), 0) * replace_na(ifelse(PAD645 < 7777, PAD645, 0), 0) * 4) +  
      (replace_na(ifelse(PAQ655 %in% 1:7, PAQ655, 0), 0) * replace_na(ifelse(PAD660 < 7777, PAD660, 0), 0) * 8) +  
      (replace_na(ifelse(PAQ670 %in% 1:7, PAQ670, 0), 0) * replace_na(ifelse(PAD675 < 7777, PAD675, 0), 0) * 4),

    # --- 临床核心预测指数 ---
    Ratio_18_16 = DR1TS180 / (DR1TS160 + 0.0001), 
    TyG_Index = log(LBXTR * LBXGLU / 2),
    SII_Index = (LBXPLTSI * LBDNENO) / (LBDLYMNO + 0.0001),
    High_CRP = factor(ifelse(LBXHSCRP > 3, "High", "Normal"), levels = c("Normal", "High")),
    High_CRP_num = ifelse(LBXHSCRP > 3, 1, 0), # 供 Logistic 使用的 0/1 变量
    Age = RIDAGEYR,
    Gender = as.factor(RIAGENDR)
  )

message(glue::glue(">>> 第一步基础清洗完成！获得纯净样本 {nrow(df_clean)} 人。"))

# 💾 留痕：保存第一阶段清洗后的全变量数据
write.csv(df_clean, "results/data/Stage1_NHANES_Cleaned_Raw_with_Covariates.csv", row.names = FALSE)

# =============================================================================
# 补丁：生成终极版 Table 1 (包含所有底层代谢指标 + 新增混杂因素)
# =============================================================================
if (!require("tableone")) install.packages("tableone")
library(tableone)

message("\n>>> 📊 正在生成满足顶级期刊要求的终极版 Table 1...")

myVars <- c(
  "Age", "Gender",
  "BMI", "Smoking", "Alcohol_g", "MET_min_week", "Diabetes", "Energy_Kcal",
  "DR1TS180", "DR1TS160", "Ratio_18_16",
  "LBXGLU", "LBXTR", "TyG_Index",
  "LBXPLTSI", "LBDNENO", "LBDLYMNO", "SII_Index"
)

catVars <- c("Gender", "Smoking", "Diabetes")

nonNormalVars <- c(
  "Age", "BMI", "Alcohol_g", "MET_min_week", "Energy_Kcal",
  "DR1TS180", "DR1TS160", "Ratio_18_16",
  "LBXGLU", "LBXTR", "TyG_Index",
  "LBXPLTSI", "LBDNENO", "LBDLYMNO", "SII_Index"
)

tab1_ult <- CreateTableOne(vars = myVars, strata = "High_CRP", data = df_clean, factorVars = catVars)

cat("\n==================================================================\n")
cat("📋 终极版 Table 1: Baseline Characteristics (包含所有新增协变量)\n")
cat("==================================================================\n")
print(tab1_ult, nonnormal = nonNormalVars, showAllLevels = TRUE)

tab1_mat_ult <- print(tab1_ult, nonnormal = nonNormalVars, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_mat_ult, file = "results/data/Table1_Ultimate_Matched.csv")

message(">>> 🎉 终极版 Table 1 已经成功导出至 results/data 目录！")

df_features <- df_clean

# =============================================================================
# 第二阶段：特征工程与 ML 探路
# =============================================================================
if (!require("randomForest")) install.packages("randomForest")
library(randomForest)

message("\n>>> 🚀 [Stage 2] 开始计算多因素网络的核心临床指数 (Feature Engineering)...")

df_features <- df_clean %>%
  mutate(
    Ratio_18_16 = DR1TS180 / (DR1TS160 + 0.0001), 
    TyG_Index = log(LBXTR * LBXGLU / 2),
    SII_Index = (LBXPLTSI * LBDNENO) / (LBDLYMNO + 0.0001),
    High_CRP = factor(ifelse(LBXHSCRP > 3, "High", "Normal"), levels = c("Normal", "High")),
    Age = RIDAGEYR,
    Gender = as.factor(RIAGENDR)
  )

# 💾 留痕：保存第二阶段特征工程数据
write.csv(df_features, "results/data/Stage2_Features_Engineered.csv", row.names = FALSE)

message(">>> 🤖 启动 Random Forest 引擎，评估各指标对炎性微环境的非线性贡献度...")
set.seed(2026) 
rf_model <- randomForest(
  High_CRP ~ Ratio_18_16 + TyG_Index + SII_Index + Age + Gender, 
  data = df_features,
  importance = TRUE, 
  ntree = 500
)

# 💾 留痕：保存第二阶段训练好的 RF 模型
saveRDS(rf_model, "results/data/Stage2_RandomForest_Model.rds")

# 画并保存 Gini 重要性图
pdf("results/plots/Stage2_RF_Gini_Importance.pdf", width = 8, height = 6)
varImpPlot(rf_model, main = "Random Forest Variable Importance for High hs-CRP")
dev.off()


# =============================================================================
# 第三阶段：SHAP 黑盒可解释性分析
# =============================================================================
p_load(fastshap, shapviz, ggplot2)

message("\n>>> 🚀 [Stage 3] 启动 SHAP 解释引擎 (计算需要几十秒，请稍候)...")

X_features <- df_features %>% select(Ratio_18_16, TyG_Index, SII_Index, Age, Gender)

pfun <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, "High"]
}

set.seed(2026)
shap_values <- explain(
  rf_model, 
  X = X_features, 
  pred_wrapper = pfun, 
  nsim = 50
)

sv <- shapviz(shap_values, X = X_features)

# 绘制 SHAP 图表
p1 <- sv_importance(sv, kind = "bar", show_numbers = TRUE) +
  theme_bw() +
  ggtitle("SHAP Feature Importance")

p2 <- sv_importance(sv, kind = "beeswarm") +
  theme_bw() +
  ggtitle("SHAP Beeswarm Plot: Impact on High hs-CRP Risk")

# 💾 留痕：导出高分辨率 SHAP 图表
ggsave("results/plots/Stage3_SHAP_Bar_Importance.pdf", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("results/plots/Stage3_SHAP_Beeswarm_Direction.pdf", plot = p2, width = 8, height = 6, dpi = 300)

message("\n>>> 🎉 全流程运行完毕！所有代码、数据 (CSV/RDS) 与高清图表已自动留痕至 /Users/bing/AA/results 目录！")

# =============================================================================
# 第四阶段：传统临床统计学落地 - 限制性立方样条 (RCS) 剂量反应分析
# =============================================================================

# 1. 加载高阶统计建模包 rms (Regression Modeling Strategies)
if (!require("rms")) install.packages("rms")
p_load(rms, ggplot2, dplyr)

message("\n>>> 🚀 [Stage 4] 启动传统临床统计学引擎 (RCS 建模)...")

# 2. 数据预处理：将结局转为 rms 包需要的 0/1 数值型
df_rcs <- df_features %>%
  mutate(
    High_CRP_num = ifelse(High_CRP == "High", 1, 0)
  )

# 3. 设定数据分布环境 (rms 包建模的必需前置步骤)
dd <- datadist(df_rcs)
options(datadist = "dd")

# 4. 构建包含 RCS 的多因素 Logistic 回归模型
# 校正 Age 和 Gender，探究 Ratio_18_16 的独立非线性总效应
# rcs(..., 4) 表示设置 4 个节点 (Knots)，足够捕捉复杂的非线性趋势
fit_rcs <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + Age + Gender, data = df_rcs)

# 打印模型基本信息
print(fit_rcs)

# 5. 生成用于绘图的预测数据
# 这里我们将结局转化为 Odds Ratio (OR)，并以人群的中位数作为参考点 (OR=1)
pred_rcs <- Predict(fit_rcs, Ratio_18_16, fun = exp, ref.zero = TRUE)

message(">>> 正在绘制 RCS 剂量反应曲线...")

# 6. 使用 ggplot2 绘制高颜值 RCS 曲线
p_rcs <- ggplot(pred_rcs) +
  geom_line(color = "#e74c3c", size = 1.2) + # 主曲线 (红色)
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#e74c3c") + # 95% 置信区间
  geom_hline(yintercept = 1, linetype = "dashed", color = "#2c3e50", size = 0.8) + # OR = 1 的参考线
  labs(
    title = "Dose-Response Relationship between Stearic/Palmitic Acid Ratio and High hs-CRP Risk",
    x = "Dietary Stearic / Palmitic Acid Ratio (18:0 / 16:0)",
    y = "Odds Ratio (OR) for High hs-CRP"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# 打印图表
print(p_rcs)

# 7. 💾 留痕：导出高分辨率 RCS 曲线图
ggsave("results/plots/Stage4_RCS_Curve.pdf", plot = p_rcs, width = 8, height = 6, dpi = 300)
message(">>> 🎉 RCS 曲线已生成并留痕至 results/plots 目录！")
# =============================================================================
# 补丁：导出 RCS 的具体 OR 值与 P 值统计量
# =============================================================================

message(">>> 正在导出 RCS 具体数值与非线性检验 P 值...")

# 1. 导出曲线上的所有 (X, Y) 坐标及其置信区间 (直接作为表格保存)
rcs_data_export <- as.data.frame(pred_rcs) %>%
  rename(
    Ratio_Value = Ratio_18_16,
    OR_Value = yhat,     # 预测的 OR 值
    CI_Lower = lower,    # 95% 置信区间下限
    CI_Upper = upper     # 95% 置信区间上限
  )
write.csv(rcs_data_export, "results/data/Stage4_RCS_Detailed_Values.csv", row.names = FALSE)

# 2. 计算并打印 P 值 (非常重要：用于判断是不是真的"非线性")
rcs_anova <- anova(fit_rcs)
cat("\n====================================================================\n")
cat("📊 RCS 模型的 P 值检验结果 (请写入论文 Results 部分)\n")
cat("====================================================================\n")
print(rcs_anova)

# 将 ANOVA 结果也保存下来
capture.output(print(rcs_anova), file = "results/data/Stage4_RCS_ANOVA_Pvalues.txt")
message(">>> 🎉 RCS 数值表和 P 值文本已保存至 results/data 目录！")

# =============================================================================
# 第五阶段：多维风险评估 - 构建并验证临床列线图 (Nomogram)
# =============================================================================
library(rms)
library(ggplot2)
library(dplyr)

message("\n>>> 🚀 [Stage 5] 开始整合多维变量，构建临床风险列线图...")

# 1. 动态刷新数据分布环境 (rms 建模必需的前置步骤，固化此步以确保绝对可复现)
# 显式将新纳入模型的 TyG_Index, SII_Index 以及修复的 Gender 固化进分布矩阵
dd <- datadist(df_rcs)
options(datadist = "dd")

# 2. # 构建包含多变量限制性立方样条 (RCS) 的 Logistic 回归预测模型
# 严谨性检查：
#   - 必须使用 rcs(Ratio_18_16, 4) 以继承第四阶段的非线性发现，否则模型前后矛盾
#   - 必须带回 Gender (性别) 这一具有强独立效应的协变量
#   - 必须指定 x = TRUE, y = TRUE，否则后续的 Bootstrap 重抽样验证 (calibrate) 会直接报错
fit_nomo <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + TyG_Index + SII_Index + Age + Gender,
                data = df_rcs, x = TRUE, y = TRUE)
# 🛡️ 增加：逻辑验证检查 (无需改动原有逻辑，仅用于确保你对结果的绝对掌控)
message(">>> 正在核对模型逻辑...")
if(any(is.na(coef(fit_nomo)))) {
  stop("🚨 警告：模型存在共线性导致的 NA 系数，请检查 SII 与其他变量是否共线！")
}
# 打印模型基本统计量 (用于审稿人核对各变量的 Beta 系数与 Wald 检验)
print(fit_nomo)

# 3. 构建具有非线性尺度的 Nomogram 对象
# 将预测的目标结局概率精确切分为临床常用的阶梯区间
nom <- nomogram(fit_nomo,
                fun = plogis, # 采用标准Logistic逻辑斯蒂反函数将 Log-odds 转换为真实发病概率
                fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                funlabel = "Risk of High hs-CRP (> 3 mg/L)",
                lp = FALSE) # 隐藏线性预测轴 (Linear Predictor)，大幅提升临床图表的可读性与美观度

# 4. 🎨 导出高清晰度多维临床列线图
pdf("results/plots/Stage5_Multivariable_Nomogram.pdf", width = 11, height = 8)
# 稍微调整画布边缘与字体缩放系数，防止变量标签由于非线性弯曲而发生重叠
plot(nom, cex.var = 0.85, cex.axis = 0.8, lmgp = 0.2)
title(main = "Nomogram for Predicting High Systemic Inflammation Risk", line = 1.5, cex.main = 1.2)
dev.off()

message(">>> 🎉 包含非线性项与全协变量的列线图已成功生成！")
message(">>> 🔄 正在启动内部模型验证：基于 1,000 次 Bootstrap 重抽样的校准曲线 (Calibration)...")

# 5. 模型内部验证：高阶校准曲线 (Calibration Curve)
# 设置确定的随机数种子，确保任何学者在任何机器上重抽样 1000 次的结果完全一致
set.seed(2026)
cal <- calibrate(fit_nomo, method = "boot", B = 1000)

# 6. 🎨 导出严谨的临床校准度验证图
pdf("results/plots/Stage5_Calibration_Curve.pdf", width = 7, height = 7)
plot(cal, 
     xlab = "Nomogram Predicted Probability", 
     ylab = "Actual Observed Probability",
     main = "Calibration Curve for the Inflammation Risk Model",
     subtitles = FALSE, # 关闭默认的冗余子标题，保持工业级排版洁净度
     cex.axis = 0.9,
     cex.lab = 1.0)
dev.off()

# 7. 💾 区分度留痕：提取并量化模型的 C-index (等价于多因素模型下的 AUC)
c_index <- fit_nomo$stats["C"]

cat("\n==================================================================\n")
cat(glue::glue("🌟 经非线性校正后的预测模型 C-index (区分度 AUC) 为: {round(c_index, 3)}\n"))
cat("   (评估基准: C-index > 0.70 提示模型具备优良的临床预测与鉴别能力)\n")
cat("==================================================================\n")

# 将模型整体的各项评估指标 (含 R2, Brier 分数, C-index) 完整导出为文本
sink("results/data/Stage5_Model_Performance_Stats.txt")
print(fit_nomo$stats)
sink()

message(">>> 🎉 [Stage 5] 全部核心临床落地图表及可复现统计量已成功保存至 results/ 目录！")
# =============================================================================

# =============================================================================
# 第六阶段：临床净收益评估 - 决策曲线分析 (DCA)
# =============================================================================
# 首次运行请先安装 dcurves 包 (如果已安装请忽略)
# install.packages("dcurves")

library(dcurves)
library(ggplot2)
library(dplyr)

message("\n>>> 🚀 [Stage 6] 启动终极临床价值评估 (DCA 建模)...")

# 1. 提取第五阶段 Nomogram 模型对每个患者预测的“真实发病概率”
# 这里直接从 fit_nomo 提取拟合值 (fitted probabilities)
df_rcs$pred_prob <- predict(fit_nomo, type = "fitted")

# 2. 构建 DCA 模型
# 评估我们的 Nomogram 模型在 0 到 1 的各种风险阈值下的临床净收益
dca_model <- dca(High_CRP_num ~ pred_prob, 
                 data = df_rcs,
                 thresholds = seq(0, 1, by = 0.01),
                 label = list(pred_prob = "Nomogram Prediction Model"))

# 3. 绘制高颜值的 DCA 曲线
p_dca <- plot(dca_model) +
  labs(
    title = "Decision Curve Analysis for the Inflammation Risk Nomogram",
    x = "Threshold Probability",
    y = "Net Benefit"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# 4. 💾 导出高分辨率 DCA 曲线
pdf("results/plots/Stage6_DCA_Curve.pdf", width = 8, height = 6)
print(p_dca)
dev.off()

message(">>> 🎉 [Stage 6] DCA 决策曲线已生成并留痕至 results/plots 目录！")
# =======================================================================

# =======================================================================
# 补丁：导出多因素 Logistic 回归的 OR 和 95% CI (完美避开 RCS 陷阱版)
# =============================================================================
message("\n>>> 📊 正在强制提取并美化多因素回归的 OR 值...")

# 1. 提取底层数据
betas <- coef(fit_nomo)
ses <- sqrt(diag(vcov(fit_nomo)))

# 2. 计算 OR 值与 95% CI
or_values <- exp(betas)
ci_lower  <- exp(betas - 1.96 * ses)
ci_upper  <- exp(betas + 1.96 * ses)

table_s1 <- data.frame(
  Variable = names(betas),
  Odds_Ratio = sprintf("%.3f", or_values),
  CI_Lower   = sprintf("%.3f", ci_lower),
  CI_Upper   = sprintf("%.3f", ci_upper)
)

table_s1$OR_with_95CI <- paste0(table_s1$Odds_Ratio, " (", table_s1$CI_Lower, "-", table_s1$CI_Upper, ")")
table_s1_clean <- table_s1[, c("Variable", "OR_with_95CI")]

# ---------------------------------------------------------
# 3. 终极洗表操作：剔除无意义的数学节点与截距项！
# ---------------------------------------------------------
# 删除包含 Intercept 或 Ratio_18_16 的所有行
table_s1_clean <- table_s1_clean[!grepl("Intercept|Ratio_18_16", table_s1_clean$Variable), ]

# 手动添加一行极其专业的非线性变量说明
ratio_row <- data.frame(
  Variable = "Dietary 18:0/16:0 Ratio", 
  OR_with_95CI = "Non-linear association (Refer to Figure 2 for dose-response curve)"
)

# 合并并重命名规范的变量名
table_s1_final <- rbind(ratio_row, table_s1_clean)

# 让变量名字更好看
table_s1_final$Variable <- c(
  "Dietary 18:0/16:0 Ratio",
  "TyG Index",
  "SII Index",
  "Age (years)",
  "Gender (Female vs. Male)"
)

# 4. 导出为 CSV
write.csv(table_s1_final, "results/data/Supplementary_Table_S1.csv", row.names = FALSE)
message(">>> 🎉 Supplementary Table S1 已经完美剔除数学节点并重新导出！")

# =============================================================================
# 终极除错版 Stage 7：彻底解决 rms 的 "length zero" 报错
# =============================================================================
library(rms)
library(dplyr)
library(tidyr)

message("\n>>> 🚀 开始执行修正数据类型的敏感性分析...")

# 1. 彻底修复数据类型（解决 anova 报错的绝对元凶！）
# rms 包的 anova 函数极度挑剔，分类变量必须是 factor，且必须清除未使用的 levels
df_sens_clean <- df_clean %>%
  mutate(
    # 将第三方要求的协变量重新提取，并强制转换为 factor
    Diabetes_Clean = factor(case_when(
      as.character(DIQ010) %in% c("1", "Yes") ~ "Yes",
      as.character(DIQ010) %in% c("2", "3", "No", "Borderline") ~ "No",
      TRUE ~ "Unknown"
    )),
    Smoking_Clean = factor(case_when(
      as.character(SMQ020) %in% c("2", "No") ~ "Never",
      as.character(SMQ020) %in% c("1", "Yes") & as.character(SMQ040) %in% c("3", "Not at all") ~ "Former",
      as.character(SMQ020) %in% c("1", "Yes") & as.character(SMQ040) %in% c("1", "2", "Every day", "Some days") ~ "Current",
      TRUE ~ "Unknown"
    )),
    Gender_Clean = factor(RIAGENDR),
    
    # 对缺失的连续变量进行中位数插补，保住核心样本量！
    BMI_imp = ifelse(is.na(BMI), median(BMI, na.rm = TRUE), BMI),
    Alcohol_imp = ifelse(is.na(Alcohol_g), 0, Alcohol_g),
    MET_imp = ifelse(is.na(MET_min_week), median(MET_min_week, na.rm = TRUE), MET_min_week)
  )

# =============================================================================
# 2. 真实敏感性分析 2：严苛混杂因素调整
# =============================================================================
# 我们剔除这两个新协变量不明的人群，并清除空因子（droplevels 是防崩溃核心）
df_sens2 <- df_sens_clean %>%
  filter(Smoking_Clean != "Unknown" & Diabetes_Clean != "Unknown") %>%
  droplevels()

# 重新打包数据环境
dd_sens2 <- datadist(df_sens2)
options(datadist = "dd_sens2")

# 运行模型
fit_confounders <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + Age + Gender_Clean + 
                         BMI_imp + Smoking_Clean + Alcohol_imp + MET_imp + Diabetes_Clean, 
                       data = df_sens2)

cat("\n======================================================\n")
cat("✅ 真实敏感性分析 2: 调整所有生活方式/代谢混杂因素\n")
cat("======================================================\n")
print(anova(fit_confounders))

# =============================================================================
# 3. 真实敏感性分析 1：能量校正（协变量法）
# =============================================================================
df_energy <- df_sens_clean %>% 
  filter(!is.na(Energy_Kcal) & Energy_Kcal > 500 & Energy_Kcal < 5000) %>%
  droplevels()
  
dd_energy <- datadist(df_energy)
options(datadist = "dd_energy")

fit_energy <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + Age + Gender_Clean + Energy_Kcal, data = df_energy)

cat("\n======================================================\n")
cat("✅ 真实敏感性分析 1: 能量摄入校正后的模型\n")
cat("======================================================\n")
print(anova(fit_energy))

# =============================================================================
# 修正版：真实敏感性分析 3 (使用 75% 分位数作为极高危截点)
# =============================================================================
message("\n>>> 🚀 开始执行敏感性分析 3: 阈值替换法...")

# 1. 计算 hs-CRP 的 75% 分位数
cutoff_75 <- quantile(df_clean$LBXHSCRP, 0.75, na.rm = TRUE)
cat("当前队列的 hs-CRP 75% 分位数为:", round(cutoff_75, 2), "mg/L\n")

# 2. 生成新的更为严苛的结局变量
df_sens3 <- df_clean %>% 
  mutate(High_CRP_75 = ifelse(LBXHSCRP > cutoff_75, 1, 0))

dd_sens3 <- datadist(df_sens3)
options(datadist = "dd_sens3")

# 3. 运行模型 (注意：这里将变量名改回了原始的 Gender)
fit_75 <- lrm(High_CRP_75 ~ rcs(Ratio_18_16, 4) + Age + Gender, data = df_sens3)

cat("\n======================================================\n")
cat("✅ 真实敏感性分析 3: 使用 75% 分位数作为结局\n")
cat("======================================================\n")
print(rms:::anova.rms(fit_75))
