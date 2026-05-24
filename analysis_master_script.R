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
# 第一阶段：提取底层核心指标 (涵盖饮食、生化、血常规)
# =============================================================================
if (!require("pacman")) install.packages("pacman")
p_load(nhanesA, dplyr, tidyr)

message(">>> 🚀 [Stage 1] 开始从 CDC 在线拉取底层原始数据...")

demo <- nhanes("DEMO_J")   
diet <- nhanes("DR1TOT_J") 
cbc  <- nhanes("CBC_J")    
glu  <- nhanes("GLU_J")    
tg   <- nhanes("TRIGLY_J") 
crp  <- nhanes("HSCRP_J")  

demo_sub <- demo %>% select(SEQN, RIDAGEYR, RIAGENDR)
diet_sub <- diet %>% select(SEQN, DR1TS180, DR1TS160)
cbc_sub  <- cbc  %>% select(SEQN, LBXPLTSI, LBDNENO, LBDLYMNO) 
glu_sub  <- glu  %>% select(SEQN, LBXGLU)
tg_sub   <- tg   %>% select(SEQN, LBXTR)
crp_sub  <- crp  %>% select(SEQN, LBXHSCRP)

df_raw <- demo_sub %>%
  left_join(diet_sub, by = "SEQN") %>%
  left_join(cbc_sub, by = "SEQN") %>%
  left_join(glu_sub, by = "SEQN") %>%
  left_join(tg_sub, by = "SEQN") %>%
  left_join(crp_sub, by = "SEQN")

df_clean <- df_raw %>%
  drop_na() %>%
  filter(RIDAGEYR >= 20) 

message(glue::glue(">>> 第一步基础清洗完成！获得纯净样本 {nrow(df_clean)} 人。"))

# 💾 留痕：保存第一阶段原始清洗数据
write.csv(df_clean, "results/data/Stage1_NHANES_Cleaned_Raw.csv", row.names = FALSE)

# =============================================================================
# 补丁：生成并导出标准 Table 1 (基线特征表)
# =============================================================================
if (!require("tableone")) install.packages("tableone")
library(tableone)

message("\n>>> 📊 正在生成标准的 Table 1 基线特征统计表...")

# 1. 明确我们需要放入 Table 1 的变量
# (注意：必须使用第二阶段生成的 df_features，因为那里才有算好的 TyG, SII 等指数)
myVars <- c("Age", "Gender", "Ratio_18_16", "TyG_Index", "SII_Index")

# 2. 明确哪些是分类变量
catVars <- c("Gender")

# 3. 按结局变量 (High_CRP) 进行分组统计
tab1 <- CreateTableOne(vars = myVars, strata = "High_CRP", data = df_features, factorVars = catVars)

# 4. 打印结果到控制台
# (因为比值、TyG、SII 在临床上通常呈非正态分布，我们指定它们输出中位数 [四分位距])
cat("\n====================================================================\n")
cat("📋 Table 1: Baseline Characteristics (请将此表内容放入论文)\n")
cat("====================================================================\n")
print(tab1, nonnormal = c("Ratio_18_16", "TyG_Index", "SII_Index"), showAllLevels = TRUE)

# 5. 自动导出为 CSV 格式，方便你直接复制进 Word 或 LaTeX 模板
tab1_mat <- print(tab1, nonnormal = c("Ratio_18_16", "TyG_Index", "SII_Index"), 
                  quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_mat, file = "results/data/Table1_Baseline_Characteristics.csv")

message(">>> 🎉 Table 1 已经成功导出到 results/data 目录！")

# =============================================================================
# 补丁：生成终极版 Table 1 (包含所有底层基础代谢/免疫/膳食指标)
# =============================================================================
if (!require("tableone")) install.packages("tableone")
library(tableone)

message("\n>>> 📊 正在生成终极镜像版 Table 1...")

# 1. 自动合并原始底层数据与分组/计算指标数据
# (确保既有血常规/饮食的绝对值，又有 TyG, SII 和 High_CRP)
df_table1_ultimate <- merge(df_clean, df_features[, c("SEQN", "High_CRP", "TyG_Index", "SII_Index")], by = "SEQN")

# 2. 罗列需要展示的所有变量 (严格对应 Methods 描述)
myVars <- c(
  "RIDAGEYR", "RIAGENDR",              # 人口学
  "DR1TS180", "DR1TS160",              # 膳食：硬脂酸、棕榈酸绝对摄入量 (克)
  "LBXGLU", "LBXTR", "TyG_Index",      # 代谢：空腹血糖、甘油三酯、TyG
  "LBXPLTSI", "LBDNENO", "LBDLYMNO", "SII_Index" # 免疫：血小板、中性、淋巴、SII
)

catVars <- c("RIAGENDR")

# 因为真实世界的血清和饮食数据往往极度偏态，我们统一使用中位数 [IQR] 和非参数检验，这在统计学上是最严谨的
nonNormalVars <- c(
  "RIDAGEYR", "DR1TS180", "DR1TS160", 
  "LBXGLU", "LBXTR", "TyG_Index", 
  "LBXPLTSI", "LBDNENO", "LBDLYMNO", "SII_Index"
)

# 3. 生成统计表
tab1_ult <- CreateTableOne(vars = myVars, strata = "High_CRP", data = df_table1_ultimate, factorVars = catVars)

cat("\n====================================================================\n")
cat("📋 终极版 Table 1: Baseline Characteristics (与 Methods 完美镜像)\n")
cat("====================================================================\n")
print(tab1_ult, nonnormal = nonNormalVars, showAllLevels = TRUE)

# 4. 导出 CSV
tab1_mat_ult <- print(tab1_ult, nonnormal = nonNormalVars, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_mat_ult, file = "results/data/Table1_Ultimate_Matched.csv")
message(">>> 🎉 终极版 Table 1 已经成功导出至 results/data 目录！")


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

# 2. 构建包含非线性项与完整协变量的多因素 Logistic 回归模型
# 严谨性检查：
#   - 必须使用 rcs(Ratio_18_16, 4) 以继承第四阶段的非线性发现，否则模型前后矛盾
#   - 必须带回 Gender (性别) 这一具有强独立效应的协变量
#   - 必须指定 x = TRUE, y = TRUE，否则后续的 Bootstrap 重抽样验证 (calibrate) 会直接报错
fit_nomo <- lrm(High_CRP_num ~ rcs(Ratio_18_16, 4) + TyG_Index + SII_Index + Age + Gender,
                data = df_rcs, x = TRUE, y = TRUE)

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
````

````bash
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
# =============================================================================
````
````bash
# =============================================================================
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
