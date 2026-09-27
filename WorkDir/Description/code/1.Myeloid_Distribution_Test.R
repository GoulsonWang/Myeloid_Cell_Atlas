#2025/11/27
# 分析Response和PreOrPost对髓系免疫细胞比例的影响
# 本脚本使用混合效应模型来分析纵向数据，其中Patient.ID作为随机效应



#input
#   1.Description1(clean).RDS     注意改正P039患者的Response.Comprehensive值，Na改为No

#output
#   0.Description1(clean)_model_InputData.RDS   输入数据(PreOrPost列为Pre或Post、Response.Comprehensive为Yes或No的样本)
#   1.macrophage_proportion_boxplot.pdf   显示不同组合下的髓系细胞比例
#   2.model_residuals.pdf   残差图————考虑了交互项的model1
#   3.model2_residuals.pdf   残差图————删除了交互项的model2
#   4.patient_trajectories.pdf   显示每个患者的轨迹图
#   5.model2.RDS  model2
#   6.model2_summary.RDS    model2的summary



# 加载必要的库
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(tidyr))
suppressMessages(library(gridExtra))  # 添加gridExtra用于组合多个图

# 设置工作目录
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code/")

# 定义路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")

# 载入数据
SData <- readRDS(file.path(data_path, "Description1(clean).RDS"))
dim(SData)



# 改正P039患者的Response.Comprehensive值错误
idx_to_modify <- which(SData@meta.data$Patient.ID == "P039" & is.na(SData@meta.data$Response.Comprehensive))
SData@meta.data[idx_to_modify, "Response.Comprehensive"] <- "No"

# 仅保留PreOrPost列为Pre或Post、Response.Comprehensive为Yes或No的样本
SData <- subset(SData, PreOrPost %in% c("Pre", "Post") & Response.Comprehensive %in% c("Yes", "No"))
#saveRDS(SData, file.path(data_path, "Description1(clean)_model_InputData.RDS"))
#SData <- readRDS(file.path(data_path, "Description1(clean)_model_InputData.RDS"))    #该行用做调试用
# 计算每个orig.ident中Macrophage细胞的比例
sample_stats <- SData@meta.data %>%     #各样本下髓系细胞比例
  group_by(orig.ident) %>%
  summarize(
    total_cells = n(),
    macrophage_cells = sum(Lineage == "Macrophage", na.rm = TRUE),
    proportion = macrophage_cells / total_cells
  )

# 将样本统计信息与元数据合并
meta_data <- merge(sample_stats, unique(SData@meta.data[, c("orig.ident", "Response.Comprehensive", "PreOrPost", "Patient.ID")]), by = "orig.ident")

# 检查重复测量情况
cat("检查重复测量情况:\n")
duplicate_check <- meta_data %>%
  count(Patient.ID, PreOrPost) %>%
  filter(n > 1)

if (nrow(duplicate_check) > 0) {
  cat("发现以下患者在同个时间点有重复测量:\n")
  print(duplicate_check)
  
  # 对于重复测量，计算均值
  cat("对重复测量数据取均值...\n")
  meta_data <- meta_data %>%
    group_by(Patient.ID, PreOrPost, Response.Comprehensive) %>%
    summarize(
      proportion = mean(proportion),
      .groups = "drop"
    )
} else {
  cat("未发现重复测量情况\n")
}

# 提供两种数据分析策略选项：
# 策略：使用所有可用数据（包括只有Pre或只有Post的患者）

cat("\n数据策略选项:\n")
cat("1. 完整数据策略: 仅分析同时具有Pre和Post测量值的患者\n")
cat("2. 所有数据策略: 使用所有可用数据，包括仅有Pre或Post测量值的患者\n")

# 查看每个Patient.ID在不同时间点的样本数
patient_timepoints <- meta_data %>%
  count(Patient.ID, PreOrPost) %>%
  tidyr::pivot_wider(names_from = PreOrPost, values_from = n, values_fill = 0)

# 显示数据分布情况
cat("\n数据分布情况:\n")
print(patient_timepoints, n = 50)

# 统计各类患者数量
only_pre <- sum(patient_timepoints$Pre > 0 & patient_timepoints$Post == 0)
only_post <- sum(patient_timepoints$Pre == 0 & patient_timepoints$Post > 0)
both <- sum(patient_timepoints$Pre > 0 & patient_timepoints$Post > 0)

cat("\n患者分组统计:\n")
cat("仅Pre测量患者:", only_pre, "人\n")
cat("仅Post测量患者:", only_post, "人\n")
cat("同时具有Pre和Post测量患者:", both, "人\n")
cat("总患者数:", nrow(patient_timepoints), "人\n")

# 默认使用所有数据，因为混合效应模型可以处理不匹配数据
meta_data_all <- meta_data

# 显示使用的数据概况
cat("\n使用所有可用数据进行分析:\n")
cat("数据包括", nrow(meta_data_all), "个时间点观测值\n")
cat("来自", length(unique(meta_data_all$Patient.ID)), "个患者\n")



#==============================================================================
# 第一部分：考虑交互效应的混合效应模型
#==============================================================================

cat("\n=== 考虑交互效应的混合效应模型 ===\n")
# 拟合考虑交互效应的混合效应模型
model1 <- lmer(proportion ~ PreOrPost * Response.Comprehensive + (1|Patient.ID), 
               data = meta_data_all)
summary_model1 <- summary(model1)
cat("\n混合效应模型结果 (考虑交互效应):\n")
print(summary_model1)

# 优雅地输出model1中三个回归系数的p值
cat("\nModel 1 回归系数 p 值:\n")
coef_table <- summary_model1$coefficients
coef_names <- rownames(coef_table)
# 只输出前三个系数（截距项除外的主要效应和交互效应）
for (i in 2:min(4, nrow(coef_table))) {
  coef_name <- coef_names[i]
  p_value <- coef_table[i, "Pr(>|t|)"]
  formatted_p <- ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3f", p_value))
  cat(sprintf("  %s: p %s\n", coef_name, formatted_p))
}

# 模型诊断：检查残差
# 绘制残差图
pdf(file.path(plot_path, "model1_residuals.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# 残差vs拟合值图
plot(fitted(model1), residuals(model1), 
     xlab = "fitted value", ylab = "residual error",
     main = "residual error vs fitted value")
abline(h = 0, lty = 2, col = "red")

# Q-Q图检查正态性
qqnorm(residuals(model1), main = "residual error normal distribution test")
qqline(residuals(model1), col = "red")

dev.off()

# 检查随机效应的方差
cat("\n模型1的随机效应方差成分:\n")
print(VarCorr(model1))

# 提取最佳模型的估计边际均值进行事后分析
if (!requireNamespace("emmeans", quietly = TRUE)) {
  cat("emmeans包未安装，跳过后事分析\n")
} else {
  emm <- emmeans::emmeans(model1, ~ PreOrPost * Response.Comprehensive)
  contrasts <- emmeans::contrast(emm, interaction = "pairwise")
  cat("\n估计边际均值:\n")
  print(contrasts)
}



#==============================================================================
# 第二部分：不考虑交互效应的混合效应模型
#==============================================================================

cat("\n=== 不考虑交互效应的混合效应模型 ===\n")
# 拟合不考虑交互效应的混合效应模型
model2 <- lmer(proportion ~ PreOrPost + Response.Comprehensive + (1|Patient.ID), 
               data = meta_data_all)
summary_model2 <- summary(model2)
cat("\n混合效应模型结果 (不考虑交互效应):\n")
print(summary_model2)

# 优雅地输出model2中回归系数的p值
cat("\nModel 2 回归系数 p 值:\n")
coef_table2 <- summary_model2$coefficients
coef_names2 <- rownames(coef_table2)
# 只输出主要效应（截距项除外）
for (i in 2:min(3, nrow(coef_table2))) {
  coef_name <- coef_names2[i]
  p_value <- coef_table2[i, "Pr(>|t|)"]
  formatted_p <- ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3f", p_value))
  cat(sprintf("  %s: p %s\n", coef_name, formatted_p))
}

#模型诊断：检查残差
# 绘制残差图
pdf(file.path(plot_path, "model2_residuals.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# 残差vs拟合值图
plot(fitted(model2), residuals(model2), 
     xlab = "fitted value", ylab = "residual error",
     main = "residual error vs fitted value")
abline(h = 0, lty = 2, col = "red")

# Q-Q图检查正态性
qqnorm(residuals(model2), main = "residual error normal distribution test")
qqline(residuals(model2), col = "red")

dev.off()

# 检查随机效应的方差
cat("\n模型2的随机效应方差成分:\n")
print(VarCorr(model2))



#==============================================================================
# 模型比较与结果输出
#==============================================================================

# 比较两个模型
cat("\n=== 模型比较 ===\n")
anova_result <- anova(model2, model1)
print(anova_result)
cat("如果p值大于0.05，则支持不考虑交互效应的简化模型\n")

# 输出两个模型的AIC值以便比较
cat("\n模型AIC比较:\n")
cat("考虑交互效应的模型AIC:", AIC(model1), "\n")
cat("不考虑交互效应的模型AIC:", AIC(model2), "\n")
cat("AIC值较小的模型拟合更好\n")



#==============================================================================
# 可视化结果
#==============================================================================

# 绘制箱线图显示不同组合间的差异
pdf(file.path(plot_path, "macrophage_proportion_boxplot.pdf"), width = 8, height = 6)

# 首先检查 model2 的系数名称
cat("检查 Model 2 的系数名称:\n")
print(rownames(summary_model2 $ coefficients))
cat("\n")

# 尝试查找包含 "PreOrPost" 和 "Response.Comprehensive" 的行名
pre_or_post_rows <- grep("PreOrPost", rownames(summary_model2 $ coefficients), value = TRUE)
response_rows <- grep("Response.Comprehensive", rownames(summary_model2 $ coefficients), value = TRUE)

# 检查是否找到了预期的行
pre_or_post_p <- NA
response_p <- NA
if(length(pre_or_post_rows) > 0) {
  # 如果找到多个 PreOrPost 相关行（理论上不应有，因为是二分类），取第一个
  target_pre_or_post_row <- pre_or_post_rows[1]
  cat("找到 PreOrPost 相关系数行名:", target_pre_or_post_row, "\n")
  pre_or_post_p <- summary_model2 $ coefficients[target_pre_or_post_row, "Pr(>|t|)"]
} else {
  cat("警告: 未找到 PreOrPost 相关系数行!\n")
}

if(length(response_rows) > 0) {
  # 如果找到多个 Response 相关行（理论上不应有），取第一个
  target_response_row <- response_rows[1]
  cat("找到 Response.Comprehensive 相关系数行名:", target_response_row, "\n")
  response_p <- summary_model2 $ coefficients[target_response_row, "Pr(>|t|)"]
} else {
  cat("警告: 未找到 Response.Comprehensive 相关系数行!\n")
}

# --- 检查是否成功提取p值，如果没有则给出提示 ---
if(is.na(pre_or_post_p) || is.null(pre_or_post_p)) {
  cat("错误: 无法提取 PreOrPost 的 p 值。\n")
  # 可以选择停止或使用默认值
  # stop("无法提取 PreOrPost 的 p 值，检查模型系数名称。")
  pre_or_post_p <- NA # 或者赋予一个默认值，但要注意后续处理
}

if(is.na(response_p) || is.null(response_p)) {
  cat("错误: 无法提取 Response.Comprehensive 的 p 值。\n")
  # stop("无法提取 Response.Comprehensive 的 p 值，检查模型系数名称。")
  response_p <- NA # 或者赋予一个默认值，但要注意后续处理
}



# 创建森林图来可视化模型系数
cat("\n=== 绘制模型系数森林图 ===\n")

# 提取系数表
coef1 <- summary_model1$coefficients
coef2 <- summary_model2$coefficients

# 获取系数名称（去除截距）
coef1_names <- rownames(coef1)[-1]  # 去掉截距
coef2_names <- rownames(coef2)[-1]  # 去掉截距

# 获取估计值和置信区间
library(broom)
confint1 <- confint(model1)
confint2 <- confint(model2)

# 为模型1创建森林图数据
forest_data1 <- data.frame(
  term = coef1_names,
  estimate = coef1[-1, "Estimate"],  # 去掉截距的估计值
  se = coef1[-1, "Std. Error"],
  p_value = coef1[-1, "Pr(>|t|)"],
  lower_ci = NA,
  upper_ci = NA
)

# 计算置信区间（使用更精确的自由度估算）
for(i in 1:nrow(forest_data1)) {
  est <- forest_data1$estimate[i]
  se <- forest_data1$se[i]
  # 使用更精确的自由度估算方法
  # 在混合效应模型中，自由度的计算较为复杂
  # 一种常用近似是使用残差自由度
  df_val <- length(residuals(model1)) - length(fixef(model1))  # 样本数 - 固定效应参数数
  t_val <- qt(0.025, df = df_val, lower.tail = FALSE)
  forest_data1$lower_ci[i] <- est - t_val * se
  forest_data1$upper_ci[i] <- est + t_val * se
}

# 为模型2创建森林图数据
forest_data2 <- data.frame(
  term = coef2_names,
  estimate = coef2[-1, "Estimate"],  # 去掉截距的估计值
  se = coef2[-1, "Std. Error"],
  p_value = coef2[-1, "Pr(>|t|)"],
  lower_ci = NA,
  upper_ci = NA
)

# 计算置信区间
for(i in 1:nrow(forest_data2)) {
  est <- forest_data2$estimate[i]
  se <- forest_data2$se[i]
  # 使用更精确的自由度估算方法
  df_val <- length(residuals(model2)) - length(fixef(model2))  # 样本数 - 固定效应参数数
  t_val <- qt(0.025, df = df_val, lower.tail = FALSE)
  forest_data2$lower_ci[i] <- est - t_val * se
  forest_data2$upper_ci[i] <- est + t_val * se
}

# 创建带p值的标签数据框
forest_data1$label <- paste0(rownames(forest_data1), "\n(p=", 
                             ifelse(forest_data1$p_value < 0.001, "<0.001", 
                                    sprintf("%.3f", forest_data1$p_value)), ")")
forest_data2$label <- paste0(rownames(forest_data2), "\n(p=", 
                             ifelse(forest_data2$p_value < 0.001, "<0.001", 
                                    sprintf("%.3f", forest_data2$p_value)), ")")

# 添加显著性指示列
forest_data1$significance <- ifelse(forest_data1$lower_ci <= 0 & forest_data1$upper_ci >= 0, "Not Significant", "Significant")
forest_data2$significance <- ifelse(forest_data2$lower_ci <= 0 & forest_data2$upper_ci >= 0, "Not Significant", "Significant")

# 为每个数据框添加模型标识
forest_data1$model <- "Model 1 (with Interaction)"
forest_data2$model <- "Model 2 (without Interaction)"

# 合并数据
combined_data <- rbind(
  data.frame(
    estimate = forest_data1$estimate,
    lower_ci = forest_data1$lower_ci,
    upper_ci = forest_data1$upper_ci,
    label = forest_data1$label,
    significance = forest_data1$significance,
    model = forest_data1$model,
    term = rownames(forest_data1),
    stringsAsFactors = FALSE
  ),
  data.frame(
    estimate = forest_data2$estimate,
    lower_ci = forest_data2$lower_ci,
    upper_ci = forest_data2$upper_ci,
    label = forest_data2$label,
    significance = forest_data2$significance,
    model = forest_data2$model,
    term = rownames(forest_data2),
    stringsAsFactors = FALSE
  )
)

# 为y轴标签重新编码term，确保它们在绘图中按正确顺序显示
combined_data$term_factor <- factor(combined_data$term, levels = unique(c(rownames(forest_data1), rownames(forest_data2))))

# 创建模型1的森林图
p1 <- ggplot(subset(combined_data, model == "Model 1 (with Interaction)"), 
             aes(x = estimate, y = term_factor, color = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # 零值参考线
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), color = "black", height = 0.2) +  # 置信区间线段
  geom_point(size = 3, shape = 15) +  # 正方形点
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "red")) +
  labs(title = "Model 1 Coefficients (with Interaction)",
       x = "Coefficient Estimate",
       y = "Variables") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# 创建模型2的森林图
p2 <- ggplot(subset(combined_data, model == "Model 2 (without Interaction)"), 
             aes(x = estimate, y = term_factor, color = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # 零值参考线
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), color = "black", height = 0.2) +  # 置信区间线段
  geom_point(size = 3, shape = 15) +  # 正方形点
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "red")) +
  labs(title = "Model 2 Coefficients (without Interaction)",
       x = "Coefficient Estimate",
       y = "Variables") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# 获取所有系数估计值的范围，以确保两个图使用相同的x轴范围
x_range <- c(min(c(forest_data1$lower_ci, forest_data2$lower_ci), na.rm = TRUE),
             max(c(forest_data1$upper_ci, forest_data2$upper_ci), na.rm = TRUE))

# 应用统一的x轴范围
p1 <- p1 + xlim(x_range)
p2 <- p2 + xlim(x_range)

# 组合并保存图片
pdf(file.path(plot_path, "model_coefficients_forest_plot.pdf"), width = 8, height = 12)
gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

cat("森林图已保存至:", file.path(plot_path, "model_coefficients_forest_plot.pdf"), "\n")


#==============================================================================
# 保存模型结果
#==============================================================================

# 保存考虑交互效应模型的结果到data_path
saveRDS(summary_model1, file = file.path(data_path, "model1_summary.rds"))
cat("\n考虑交互效应的模型结果已保存至:", file.path(data_path, "model1_summary.rds"), "\n")

# 保存不考虑交互效应模型的结果到data_path
saveRDS(summary_model2, file = file.path(data_path, "model2_summary.rds"))
cat("\n不考虑交互效应的模型结果已保存至:", file.path(data_path, "model2_summary.rds"), "\n")

# 同时保存模型本身，便于后续进一步分析
saveRDS(model1, file = file.path(data_path, "model1.rds"))
cat("考虑交互效应的模型已保存至:", file.path(data_path, "model1.rds"), "\n")

saveRDS(model2, file = file.path(data_path, "model2.rds"))
cat("不考虑交互效应的模型已保存至:", file.path(data_path, "model2.rds"), "\n")



#==============================================================================
# 输出结果总结
#==============================================================================

cat("\n=== 分析总结 ===\n")
cat("1. 使用混合效应模型分析了Response和PreOrPost对Macrophage比例的影响\n")
cat("2. 考虑了患者内部相关性(Patient.ID作为随机效应)\n")
cat("3. 检测了两个因子的交互作用\n")
cat("4. 通过AIC、BIC、似然比检验等方法比较了不同模型\n")
cat("5. 生成了可视化图表辅助结果解释\n")
cat("6. 利用了所有可用数据（包括不完整数据）以提高统计效力\n")