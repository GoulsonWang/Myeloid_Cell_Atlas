#2025/11/27     
#本脚本用于运用多元logistic回归分析不同变量(CancerType、Response、PreOrPost)对细胞亚类分布的影响



#input
#   1.SData_Cluster030.RDS
#   2.

#output
#   1.


# 加载必要的库
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(nnet))  # 用于多项logistic回归
suppressMessages(library(tidyr))  # 用于数据整理



# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code/")



# 载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/SData_Cluster030.RDS")    #6GB
colnames(SData@meta.data)



# 查看数据基本信息
cat("Seurat对象基本信息:\n")
cat("细胞总数:", ncol(SData), "\n")
cat("基因总数:", nrow(SData), "\n")

# 查看关键变量的信息
cluster_ids <- as.character(SData@meta.data$SCT_snn_res.0.3)
cancer_types <- as.character(SData@meta.data$CancerType)
response <- as.character(SData@meta.data$Response.Comprehensive)
pre_or_post <- as.character(SData@meta.data$PreOrPost)

cat("细胞亚类(cluster)数量:", length(unique(cluster_ids)), "\n")
cat("不同癌症类型数量:", length(unique(cancer_types)), "\n")
cat("不同Response数量:", length(unique(response)), "\n")
cat("不同PreOrPost数量:", length(unique(pre_or_post)), "\n")
cat("各亚类名称:", sort(unique(cluster_ids)), "\n")
cat("各癌症类型:", sort(unique(cancer_types)), "\n")
cat("Response类型:", sort(unique(response)), "\n")
cat("PreOrPost类型:", sort(unique(pre_or_post)), "\n")

# 准备多元logistic回归分析所需的数据
# 创建一个包含所有相关变量的数据框
meta_data <- SData@meta.data %>%
  dplyr::select(SCT_snn_res.0.3, CancerType, Response.Comprehensive, PreOrPost) %>%
  dplyr::filter(!is.na(SCT_snn_res.0.3) & !is.na(CancerType) & 
                !is.na(Response.Comprehensive) & !is.na(PreOrPost))

cat("\n用于多元logistic回归分析的细胞数量:", nrow(meta_data), "\n")

# 查看数据分布
cat("\n各变量组合的细胞数量:\n")
data_summary <- meta_data %>%
  dplyr::group_by(CancerType, Response.Comprehensive, PreOrPost, SCT_snn_res.0.3) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = SCT_snn_res.0.3, values_from = count, values_fill = 0)

print(data_summary)

# 执行多元logistic回归分析
# 以细胞亚类为因变量，CancerType、Response、PreOrPost为自变量
cat("\n开始执行多元logistic回归分析...\n")

# 设置参考组（baseline）
meta_data$CancerType <- relevel(as.factor(meta_data$CancerType), ref = "ESCC")
meta_data$Response.Comprehensive <- relevel(as.factor(meta_data$Response.Comprehensive), ref = "No")
meta_data$PreOrPost <- relevel(as.factor(meta_data$PreOrPost), ref = "Pre")

# 拟合多项logistic回归模型
multinom_model <- multinom(SCT_snn_res.0.3 ~ CancerType + Response.Comprehensive + PreOrPost, 
                          data = meta_data, trace = FALSE)

# 输出模型摘要
cat("\n多元logistic回归模型结果:\n")
model_summary <- summary(multinom_model)
print(model_summary)

# 计算 Odds Ratios 和置信区间
cat("\nOdds Ratios (相对于参考组):\n")
ORs <- exp(coef(multinom_model))
print(ORs)

# 获取标准误并计算置信区间
ses <- model_summary$standard.errors
cat("\n标准误:\n")
print(ses)

# 计算95%置信区间
ci_lower <- exp(coef(multinom_model) - 1.96 * ses)
ci_upper <- exp(coef(multinom_model) + 1.96 * ses)

cat("\nOdds Ratios的95%置信区间:\n")
ci_table <- cbind(exp(coef(multinom_model)), ci_lower, ci_upper)
colnames(ci_table) <- c("OR", "CI Lower", "CI Upper")
print(ci_table)

# 执行模型诊断
cat("\n模型拟合统计信息:\n")
cat("AIC:", AIC(multinom_model), "\n")
cat("残差自由度:", df.residual(multinom_model), "\n")

# 保存模型结果
saveRDS(multinom_model, file.path(data_path, "multinom_model_result.rds"))
saveRDS(model_summary, file.path(data_path, "multinom_model_summary.rds"))
saveRDS(ORs, file.path(data_path, "odds_ratios.rds"))
saveRDS(ci_table, file.path(data_path, "confidence_intervals.rds"))

cat("\n多元logistic回归分析结果已保存到:", data_path, "\n")

# 输出重要变量的显著性结果
cat("\n=== 分析结论 ===\n")
p_values <- model_summary$coefficients / model_summary$standard.errors
significant_vars <- which(abs(p_values) > 1.96, arr.ind = TRUE)

if(nrow(significant_vars) > 0) {
  cat("发现以下变量对细胞亚类分布有显著影响 (p < 0.05):\n")
  for(i in 1:nrow(significant_vars)) {
    row_idx <- significant_vars[i, 1]
    col_idx <- significant_vars[i, 2]
    cat(sprintf("细胞亚类 %s 相对于参考类: %s (p = %f)\n", 
                rownames(p_values)[row_idx], 
                colnames(p_values)[col_idx],
                2 * pnorm(-abs(p_values[row_idx, col_idx]))))
  }
} else {
  cat("未发现对细胞亚类分布有显著影响的变量 (p >= 0.05)\n")
}

# 可视化结果
# 1. 绘制Odds Ratios森林图
# 将系数整理成适合绘图的格式
coef_df <- data.frame(
  Variable = colnames(coef(multinom_model))[-1],  # 排除截距
  Cluster = rep(rownames(coef(multinom_model)), each = ncol(coef(multinom_model))-1),
  OR = as.numeric(ORs[, -1]),  # 排除截距列
  CI_lower = as.numeric(ci_lower[, -1]),
  CI_upper = as.numeric(ci_upper[, -1])
)

# 添加显著性标记
coef_df$Significant <- ifelse(coef_df$CI_lower > 1 | coef_df$CI_upper < 1, "Yes", "No")

# 绘制森林图
png(file.path(plot_path, "odds_ratios_forest_plot.png"), width = 1000, height = 800)
ggplot(coef_df, aes(x = OR, y = reorder(paste(Cluster, Variable), OR))) +
  geom_point(aes(color = Significant), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Significant), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  labs(title = "多元logistic回归分析 - Odds Ratios森林图",
       x = "Odds Ratio (log scale)",
       y = "细胞亚类和变量组合") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()

# 2. 绘制各变量对细胞亚类分布的影响热图
# 计算每个变量在各类别中的平均Odds Ratio
or_matrix <- ORs[, -1]  # 排除截距
or_matrix_df <- as.data.frame(as.table(or_matrix))
colnames(or_matrix_df) <- c("Cluster", "Variable", "OR")

png(file.path(plot_path, "odds_ratios_heatmap.png"), width = 1000, height = 600)
ggplot(or_matrix_df, aes(x = Variable, y = Cluster, fill = log(OR))) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, name = "log(Odds Ratio)") +
  labs(title = "各变量对细胞亚类分布的影响 (log Odds Ratio)",
       x = "变量",
       y = "细胞亚类") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 3. 显示各类别中细胞数量的分布
cluster_counts <- meta_data %>%
  group_by(SCT_snn_res.0.3) %>%
  summarise(count = n(), .groups = "drop")

png(file.path(plot_path, "cluster_distribution.png"), width = 800, height = 600)
ggplot(cluster_counts, aes(x = SCT_snn_res.0.3, y = count, fill = SCT_snn_res.0.3)) +
  geom_bar(stat = "identity") +
  labs(title = "各细胞亚类的分布",
       x = "细胞亚类",
       y = "细胞数量") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

cat("\n可视化图表已保存到:", plot_path, "\n")

# 保存分析结果总结
result_summary <- list(
  Model = multinom_model,
  Summary = model_summary,
  OddsRatios = ORs,
  ConfidenceIntervals = ci_table,
  CoefficientData = coef_df,
  ClusterCounts = cluster_counts
)

saveRDS(result_summary, file.path(data_path, "multinomial_logistic_regression_results.rds"))
cat("完整分析结果已保存到:", file.path(data_path, "multinomial_logistic_regression_results.rds"), "\n")

cat("\n=== 分析完成 ===\n")
cat("多元logistic回归分析已完成，结果包括:\n")
cat("1. 模型参数估计\n")
cat("2. Odds Ratios及其置信区间\n")
cat("3. 各变量的显著性检验\n")
cat("4. 可视化图表\n")
cat("所有结果已保存到数据目录中\n")

# 构建列联表
contingency_table <- table(cancer_types, cluster_ids)
cat("\n列联表:\n")
print(contingency_table)

# 检查期望频数
chi_expected <- chisq.test(contingency_table)$expected
low_expected <- sum(chi_expected < 5)
cat("\n期望频数小于5的单元格数量:", low_expected, "\n")
cat("总单元格数量:", nrow(contingency_table) * ncol(contingency_table), "\n")

# 执行卡方检验
chi_result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
cat("\n卡方检验结果:\n")
print(chi_result)

# 计算标准化残差以识别主要差异来源
std_residuals <- chi_result$stdres
cat("\n标准化残差 (>1.96 或 <-1.96 表示显著偏离):\n")
print(std_residuals)

# 计算各癌症类型中各类细胞的比例
proportion_table <- prop.table(contingency_table, margin = 1)
cat("\n各癌症类型中各类细胞的比例:\n")
print(round(proportion_table, 4))

# 可视化结果
# 1. 堆叠柱状图显示比例
proportion_df <- as.data.frame(proportion_table)
colnames(proportion_df) <- c("CancerType", "Cluster", "Proportion")

png(file.path(plot_path, "subtype_distribution_stacked_barplot_by_cancer.png"), width = 1000, height = 600)
ggplot(proportion_df, aes(fill=Cluster, y=Proportion, x=CancerType)) + 
  geom_bar(position="fill", stat="identity") +
  labs(title="各癌症类型中细胞亚类分布比例", 
       x="CancerType", 
       y="Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 2. 热图显示标准化残差
# 将标准化残差转换为数据框
residuals_df <- as.data.frame(as.table(std_residuals))
colnames(residuals_df) <- c("CancerType", "Cluster", "StdResidual")

# 创建热图
png(file.path(plot_path, "subtype_distribution_residuals_heatmap_by_cancer.png"), width = 1000, height = 600)
ggplot(residuals_df, aes(x=CancerType, y=Cluster, fill=StdResidual)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", 
                       name="Standardized Residual") +
  labs(title="标准化残差热图 (绝对值>1.96表示显著偏离期望)",
       x="CancerType", 
       y="Cell Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 保存结果
saveRDS(contingency_table, file.path(data_path, "contingency_table_by_cancer.rds"))
saveRDS(chi_result, file.path(data_path, "chi_square_test_result_by_cancer.rds"))
saveRDS(proportion_table, file.path(data_path, "proportion_table_by_cancer.rds"))
saveRDS(std_residuals, file.path(data_path, "standardized_residuals_by_cancer.rds"))

cat("\n结果已保存到:", data_path, "\n")
cat("保存的文件包括:\n")
cat("1. contingency_table_by_cancer.rds - 列联表\n")
cat("2. chi_square_test_result_by_cancer.rds - 卡方检验结果\n")
cat("3. proportion_table_by_cancer.rds - 比例表\n")
cat("4. standardized_residuals_by_cancer.rds - 标准化残差\n")

# 输出结论
cat("\n=== 分析结论 ===\n")
cat("卡方统计量:", round(chi_result$statistic, 2), "\n")
cat("自由度:", chi_result$parameter, "\n")
cat("p值:", format.pval(chi_result$p.value, digits = 3), "\n")

if(chi_result$p.value < 0.05) {
  cat("结论: 不同癌症类型下的细胞亚类分布存在显著差异 (p < 0.05)\n")
} else {
  cat("结论: 不同癌症类型下的细胞亚类分布无显著差异 (p >= 0.05)\n")
}

# 检查是否需要Fisher精确检验的替代方案
if(low_expected > 0.2 * (nrow(contingency_table) * ncol(contingency_table))) {
  cat("警告: 超过20%的单元格期望频数小于5，卡方检验结果可能不可靠。\n")
  cat("建议考虑使用Fisher精确检验或其他替代方法。\n")
}


