#2025/12/5    
#本脚本用于构建三维列联表，每一列代表一个亚型，行属性分别为PreOrPost,ROrNR。最终目的是找出R与NR组之间主要差异的亚型.
#注意，analyze_cancer_type函数太长，导致在VScode中无法正常闭合，在命令行中通过Rscript 3.High_dimensional_contingency_table.R 来得到该脚本的输出结果.



#input
#   1.SData_Cluster030.RDS

#output
#   1.high_dimensional_contingency_table.csv    #R、NR的分组结果是Yes/No
#   2.contingency_table_ESCC.csv      以下四个数据的治疗效果评估是按照Response.RECIST Or Response.Pathologic的严格结果，属于有序变量
#   3.contingency_table_PC.csv
#   4.contingency_table_NSCLC.csv
#   5.contingency_table_TNBC1.csv     
#   6.各癌症类型下的列联表分析结果



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#2026/1/3修改：将Minor_Cell_Type的亚类结果作为输入，重新跑该脚本代码。结果均保存在Minor_Cell_Type_Input文件夹之下。


#input:/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS

#output
#   1.high_dimensional_contingency_table.csv    #R、NR的分组结果是Yes/No
#   2.contingency_table_ESCC.csv      以下四个数据的治疗效果评估是按照Response.RECIST Or Response.Pathologic的严格结果，属于有序变量
#   3.contingency_table_PC.csv
#   4.contingency_table_NSCLC.csv
#   5.contingency_table_TNBC1.csv     
#   6.各癌症类型下的列联表分析结果



# 加载必要的库
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))  # 用于数据整理



# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
data_path <- paste(file_path, "/data/Minor_Cell_Type_Input", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code/")

# 载入数据(改名后的数据，为亚型命名后的！)
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_Rename.RDS")    #6GB
colnames(SData@meta.data)
table(SData$Minor_Cell_Type)


#构建四维列联表

# 提取元数据并清理缺失值，仅保留PreOrPost为"Pre"或"Post"的内容
meta_data <- SData@meta.data %>%
  dplyr::select(PreOrPost, Response.Comprehensive, CancerType, Minor_Cell_Type) %>%
  dplyr::filter(!is.na(PreOrPost) & !is.na(Response.Comprehensive) & 
                !is.na(CancerType) & !is.na(Minor_Cell_Type)) %>%
  dplyr::filter(PreOrPost %in% c("Pre", "Post"))

# 重命名列以提高可读性
colnames(meta_data) <- c("PreOrPost", "ROrNR", "CancerType", "Subtype")

# 构建四维列联表
contingency_table <- table(meta_data$PreOrPost, meta_data$ROrNR, meta_data$CancerType, meta_data$Subtype)
dimnames(contingency_table) <- list(PreOrPost = levels(as.factor(meta_data$PreOrPost)),
                                    ROrNR = levels(as.factor(meta_data$ROrNR)),
                                    CancerType = levels(as.factor(meta_data$CancerType)),
                                    Subtype = levels(as.factor(meta_data$Subtype)))

# 显示列联表维度
cat("列联表维度:", dim(contingency_table), "\n")
# 将高维表转换为数据框格式便于输出到CSV
contingency_df <- as.data.frame.table(contingency_table)
# 写入CSV文件
write.csv(contingency_df, file = file.path(data_path, "high_dimensional_contingency_table.csv"), row.names = FALSE)
cat("四维列联表已保存到:", file.path(data_path, "high_dimensional_contingency_table.csv"), "\n")



#以CancerType分层，构建三维列联表

# 提取包含更多响应变量的元数据
meta_data_all <- SData@meta.data %>%
  dplyr::select(PreOrPost, CancerType, Minor_Cell_Type, Response.RECIST, Response.Pathologic) %>%
  dplyr::filter(!is.na(PreOrPost) & !is.na(CancerType) & !is.na(Minor_Cell_Type)) %>%
  dplyr::filter(PreOrPost %in% c("Pre", "Post"))

# 重命名列以提高可读性
colnames(meta_data_all)[colnames(meta_data_all) == "Minor_Cell_Type"] <- "Subtype"

# 获取所有唯一的癌症类型
cancer_types <- unique(meta_data_all$CancerType)

# 为每种癌症类型构建三维列联表
for (cancer_type in cancer_types) {
  # 筛选特定癌症类型的数据
  cancer_data <- meta_data_all[meta_data_all$CancerType == cancer_type, ]
  
  # 根据癌症类型选择适当的响应变量
  if (cancer_type %in% c("PC", "TNBC1")) {
    # 对于PC类型，使用Response.RECIST作为响应变量
    ror_nr_var <- cancer_data$Response.RECIST
    response_name <- "Response.RECIST"
  } else {
    # 对于其他类型，使用Response.Pathologic作为响应变量
    ror_nr_var <- cancer_data$Response.Pathologic
    response_name <- "Response.Pathologic"
  }
  
  # 构建三维列联表 (PreOrPost x Response x Subtype)
  three_d_table <- table(cancer_data$PreOrPost, ror_nr_var, cancer_data$Subtype)
  dimnames(three_d_table) <- list(PreOrPost = levels(as.factor(cancer_data$PreOrPost)),
                                  Response = levels(as.factor(ror_nr_var)),
                                  Subtype = levels(as.factor(cancer_data$Subtype)))
  
  # 将三维表转换为数据框
  three_d_df <- as.data.frame.table(three_d_table)
  
  # 生成文件名并写入CSV
  filename <- paste0("contingency_table_", cancer_type, ".csv")
  write.csv(three_d_df, file = file.path(data_path, filename), row.names = FALSE)
  
  cat(sprintf("%s癌症类型的三维列联表已保存到: %s，使用的响应变量为: %s\n", cancer_type, file.path(data_path, filename), response_name))
}



#三维列联表建模分析（线性对数模型）
library(MASS)
# 对三维列联表进行建模分析的通用函数
analyze_contingency_table <- function(filename, cancer_type) {

  # 读入数据
  data <- read.csv(file.path(data_path, filename))
  # 处理Freq中结构0
  data$Freq[data$Freq == 0] <- 1
  # 设定列属性
  data$PreOrPost <- as.factor(data$PreOrPost)
  data$Response <- as.factor(data$Response)
  data$Subtype <- as.factor(data$Subtype)
  data$Freq <- as.numeric(data$Freq)
  
  # 构建列联表模型
  data_model <- xtabs(Freq ~ ., data = data)
  
  # 全独立模型
  model1 <- loglm(~ 1 + 2 + 3, data = data_model)
  # 单交互模型
  model2 <- update(model1, .~. + 1*2)
  model3 <- update(model1, .~. + 1*3)
  model4 <- update(model1, .~. + 2*3)
  # 双交互模型
  model5 <- update(model1, .~. + 1*2 + 1*3)
  model6 <- update(model1, .~. + 1*2 + 2*3)
  model7 <- update(model1, .~. + 1*3 + 2*3)
  # 三交互模型
  model8 <- update(model1, .~. + 1*2 + 1*3 + 2*3)
  
  # 输出方差分析表
  anova_result <- anova(model1, model2, model3, model4, model5, model6, model7, model8)
  print(sprintf("%s: 模型比较ANOVA", cancer_type))
  print(anova_result)
  
  # 整理结果
  model_names <- paste0(cancer_type, "_model", 1:8)
  models <- list(model1, model2, model3, model4, model5, model6, model7, model8)
  results <- data.frame(
    CancerType = rep(cancer_type, 8),
    Model = model_names,
    Likelihood_Ratio_Stat = numeric(8),
    Chi_Square_Stat = numeric(8),
    Degrees_of_Freedom = numeric(8),
    P_Value = numeric(8)
  )
  
  for (i in 1:8) {
    smry <- summary(models[[i]])
    results$Likelihood_Ratio_Stat[i] <- smry$tests[1,1]
    results$Chi_Square_Stat[i] <- smry$tests[2,1]
    results$Degrees_of_Freedom[i] <- smry$tests[1,2]
    results$P_Value[i] <- smry$tests[1,3]
  }
  
  return(results)
}

# 分别对每种癌症类型进行分析
PC_results <- analyze_contingency_table("contingency_table_PC.csv", "PC")
ESCC_results <- analyze_contingency_table("contingency_table_ESCC.csv", "ESCC")
NSCLC_results <- analyze_contingency_table("contingency_table_NSCLC.csv", "NSCLC")
TNBC1_results <- analyze_contingency_table("contingency_table_TNBC1.csv", "TNBC1")

# 合并所有模型结果并保存
all_model_results <- rbind(PC_results, ESCC_results, NSCLC_results, TNBC1_results)
write.csv(all_model_results, file = file.path(data_path, "model_comparison_results.csv"), row.names = FALSE)
cat("模型比较结果已保存到:", file.path(data_path, "model_comparison_results.csv"), "\n")

# 输出总结果
print("所有癌症类型的模型比较结果:")
print(all_model_results)



#结果解读：
#NSCLC中，model8双交互模型（同时考虑1*2，1*3，2*3）可以拟合数据。即：两两相关模型（若按Subtype分层，则每层的相合方向与相合程度都相同）
#TNBC1中，model7双交互模型（考虑1*3 + 2*3）可以较好的拟合数据。即：给定变量3，则变量1和变量2相互独立
#其余癌症均无特殊的条件独立性成立
#TNBC model7的拟合值结果：已保存为txt文档，同该脚本同目录下
#NSCLC model8的拟合值结果：已保存为txt文档，同该脚本同目录下





# 定义通用函数来处理各种癌症类型的压缩和分层分析
analyze_cancer_type <- function(cancer_type, file_name) {
  cat("\n====================", cancer_type, "分析====================\n")
  
  # 读取数据，指定编码以避免MBCS字符问题
  cancer_data <- read.csv(file.path(data_path, file_name), encoding = "UTF-8")
  # 不再预先处理Freq中的0，改为在压缩或分层后的二维列联表中处理
  # 设定列属性
  cancer_data$PreOrPost <- as.factor(cancer_data$PreOrPost)
  cancer_data$Response <- as.factor(cancer_data$Response)
  cancer_data$Subtype <- as.factor(cancer_data$Subtype)
  cancer_data$Freq <- as.numeric(cancer_data$Freq)
  
  cat("数据维度:", dim(cancer_data), "\n")
  
  # 讨论PreOrPost * Subtype关系
  cat("\n---------- 讨论PreOrPost * Subtype关系 ----------\n")
  
  # 压缩Response维度
  compressed_response <- aggregate(Freq ~ PreOrPost + Subtype, data = cancer_data, FUN = sum)
  twod_table_pre_post_subtype <- xtabs(Freq ~ PreOrPost + Subtype, data = compressed_response)
  
  # 处理压缩后的二维列联表中的结构0
  twod_table_pre_post_subtype[twod_table_pre_post_subtype == 0] <- 1
  
  # 显示二维列联表和卡方检验结果
  cat("\n", cancer_type, "压缩后的二维列联表（行：PreOrPost，列：Subtype）:\n")
  print(twod_table_pre_post_subtype)
  chi_test_pre_post_subtype <- chisq.test(twod_table_pre_post_subtype)
  cat("\n", cancer_type, "二维列联表的卡方检验结果 (PreOrPost vs Subtype):\n")
  print(chi_test_pre_post_subtype)
  
  # 分层分析：按Response分层，分析PreOrPost和Subtype的关系
  cat("\n---------- 按Response分层分析PreOrPost * Subtype关系 ----------\n")
  response_levels <- levels(cancer_data$Response)
  cat(cancer_type, "数据中的Response水平:", response_levels, "\n")
  
  # 创建一个列表来存储每个Response层的结果
  stratified_results <- list()
  
  for (response_level in response_levels) {
    # 筛选数据并创建二维列联表
    subset_data <- cancer_data[cancer_data$Response == response_level, ]
    twod_table <- xtabs(Freq ~ PreOrPost + Subtype, data = subset_data)
    twod_table[twod_table == 0] <- 1  # 处理结构0
    
    # 显示该层的二维列联表
    cat("\n", cancer_type, "- Response level:", response_level, "\n")
    print(twod_table)
    
    # 进行卡方检验并保存结果
    if (any(twod_table > 0)) {
      chi_test_result <- chisq.test(twod_table)
      cat("卡方检验结果:\n")
      print(chi_test_result)
      stratified_results[[response_level]] <- list(table = twod_table, chi_square_test = chi_test_result)
    } else {
      cat("该层数据全为0，跳过卡方检验。\n")
    }
  }
  
  # 总结分层分析结果
  cat("\n=== ", cancer_type, "按Response分层的二维列联表分析总结 ===\n")
  for (response_level in response_levels) {
    if (!is.null(stratified_results[[response_level]])) {
      cat("\nResponse level:", response_level, "\n")
      cat("  Chi-square statistic:", stratified_results[[response_level]]$chi_square_test$statistic, "\n")
      cat("  Degrees of freedom:", stratified_results[[response_level]]$chi_square_test$parameter, "\n")
      cat("  P-value:", stratified_results[[response_level]]$chi_square_test$p.value, "\n")
    }
  }
  
  # 讨论Response和Subtype之间的关系
  cat("\n---------- 讨论Response * Subtype关系 ----------\n")
  
  # 压缩PreOrPost维度
  compressed_pre_or_post <- aggregate(Freq ~ Response + Subtype, data = cancer_data, FUN = sum)
  twod_table_response_subtype <- xtabs(Freq ~ Response + Subtype, data = compressed_pre_or_post)
  
  # 处理压缩后的二维列联表中的结构0
  twod_table_response_subtype[twod_table_response_subtype == 0] <- 1
  
  # 显示二维列联表和卡方检验结果
  cat("\n", cancer_type, "压缩后的二维列联表（行：Response，列：Subtype）:\n")
  print(twod_table_response_subtype)
  chi_test_response_subtype <- chisq.test(twod_table_response_subtype)
  cat("\n", cancer_type, "二维列联表的卡方检验结果 (Response vs Subtype):\n")
  print(chi_test_response_subtype)
  
  # 分层分析：按PreOrPost分层，分析Response和Subtype的关系
  cat("\n---------- 按PreOrPost分层分析Response * Subtype关系 ----------\n")
  pre_or_post_levels <- levels(cancer_data$PreOrPost)
  cat(cancer_type, "数据中的PreOrPost水平:", pre_or_post_levels, "\n")
  
  # 创建一个列表来存储每个PreOrPost层的结果
  pre_or_post_results <- list()
  
  for (pre_or_post_level in pre_or_post_levels) {
    # 筛选数据并创建二维列联表
    subset_data <- cancer_data[cancer_data$PreOrPost == pre_or_post_level, ]
    twod_table <- xtabs(Freq ~ Response + Subtype, data = subset_data)
    twod_table[twod_table == 0] <- 1  # 处理结构0
    
    # 显示该层的二维列联表
    cat("\n", cancer_type, "- PreOrPost level:", pre_or_post_level, "\n")
    print(twod_table)
    
    # 进行卡方检验并保存结果
    if (any(twod_table > 0)) {
      chi_test_result <- chisq.test(twod_table)
      cat("卡方检验结果:\n")
      print(chi_test_result)
      pre_or_post_results[[pre_or_post_level]] <- list(table = twod_table, chi_square_test = chi_test_result)
    } else {
      cat("该层数据全为0，跳过卡方检验。\n")
    }
  }
  
  # 总结按PreOrPost分层的分析结果
  cat("\n=== ", cancer_type, "按PreOrPost分层的二维列联表分析总结 ===\n")
  for (pre_or_post_level in pre_or_post_levels) {
    if (!is.null(pre_or_post_results[[pre_or_post_level]])) {
      cat("\nPreOrPost level:", pre_or_post_level, "\n")
      cat("  Chi-square statistic:", pre_or_post_results[[pre_or_post_level]]$chi_square_test$statistic, "\n")
      cat("  Degrees of freedom:", pre_or_post_results[[pre_or_post_level]]$chi_square_test$parameter, "\n")
      cat("  P-value:", pre_or_post_results[[pre_or_post_level]]$chi_square_test$p.value, "\n")
    }
  }
  
  # 返回结果列表
  return(list(
    data = cancer_data,
    pre_post_subtype_table = twod_table_pre_post_subtype,
    pre_post_subtype_test = chi_test_pre_post_subtype,
    response_subtype_table = twod_table_response_subtype,
    response_subtype_test = chi_test_response_subtype,
    stratified_by_response = stratified_results,
    stratified_by_pre_or_post = pre_or_post_results
  ))
}

# 对PC、ESCC、NSCLC、TNBC1分别进行分析
PC_results <- analyze_cancer_type("PC", "contingency_table_PC.csv")
ESCC_results <- analyze_cancer_type("ESCC", "contingency_table_ESCC.csv")
NSCLC_results <- analyze_cancer_type("NSCLC", "contingency_table_NSCLC.csv")
TNBC1_results <- analyze_cancer_type("TNBC1", "contingency_table_TNBC1.csv")

# 总结所有癌症类型的分析结果
cat("\n\n==================== 所有癌症类型分析总结 ====================\n")

cancer_types <- c("PC", "ESCC", "NSCLC", "TNBC1")
cancer_results <- list(PC_results, ESCC_results, NSCLC_results, TNBC1_results)

for (i in 1:length(cancer_types)) {
  cat("\n", cancer_types[i], ":\n")
  cat("  PreOrPost vs Subtype - p-value:", cancer_results[[i]]$pre_post_subtype_test$p.value, "\n")
  cat("  Response vs Subtype - p-value:", cancer_results[[i]]$response_subtype_test$p.value, "\n")
}



#重要结果汇总：

#1.PC
#PreOrPost与Subtype不独立，分层后各层均不独立
#Response与Subtype不独立，分层后各层均不独立

#2.ESCC
#PreOrPost与Subtype不独立，分层后各层均不独立
#Response与Subtype不独立，分层后各层均不独立

#3.NSCLC
#PreOrPost与Subtype不独立，分层后各层均不独立
#Response与Subtype不独立，分层后各层均不独立

#4.TNBC1
#PreOrPost与Subtype不独立，分层后Partial层独立，Stable层不独立
#Response与Subtype不独立，分层后各层均不独立
#另：证明得到给定Subtype后，PreOrPost与Response相互独立。


#总结：除了TNBC1外，其余3个数据集中Subtype与PreOrPost，Subtype与Response均不独立，故可以通过Subtype来判断Response状况，PreOrPost状况

#下一步转入对应分析，将PreOrPost与Response融合成一个变量，判断Subtype与融合变量之间的关联，找出能预测不同融合变量状态下的Subtype


# 添加：保存三维列联表建模分析和压缩分层分析的结论表格

# Create model comparison summary table
model_comparison_summary <- data.frame(
  CancerType = c("PC", "ESCC", "NSCLC", "TNBC1"),
  BestModel = c("To be determined", "To be determined", "Model 8 (Three-way interaction)", "Model 7 (Two-way interaction)"),
  ModelDescription = c("Need further analysis", "Need further analysis", 
                      "PreOrPost*Response + PreOrPost*Subtype + Response*Subtype", 
                      "PreOrPost*Subtype + Response*Subtype"),
  BiologicalInterpretation = c("Complex interactions among variables", "Complex interactions among variables",
                             "Complex interactions among three variables, relationships between any two variables are affected by the third variable",
                             "Given subtype condition, PreOrPost and Response are independent")
)

# Create association analysis summary table
association_analysis_summary <- data.frame(
  CancerType = c("PC", "ESCC", "NSCLC", "TNBC1"),
  PrePost_Subtype_Independent = c("No", "No", "No", "No"),
  Response_Subtype_Independent = c("No", "No", "No", "No"),
  PrePost_Subtype_PValue = c(PC_results$pre_post_subtype_test$p.value, 
                            ESCC_results$pre_post_subtype_test$p.value,
                            NSCLC_results$pre_post_subtype_test$p.value,
                            TNBC1_results$pre_post_subtype_test$p.value),
  Response_Subtype_PValue = c(PC_results$response_subtype_test$p.value,
                             ESCC_results$response_subtype_test$p.value,
                             NSCLC_results$response_subtype_test$p.value,
                             TNBC1_results$response_subtype_test$p.value),
  Stratified_Findings = c("All strata are dependent", "All strata are dependent", "All strata are dependent", 
                         "PreOrPost*Subtype: Partial stratum independent, Stable stratum dependent")
)

# Save tables to CSV files
write.csv(model_comparison_summary, 
          file = file.path(data_path, "model_comparison_summary.csv"), 
          row.names = FALSE)

write.csv(association_analysis_summary, 
          file = file.path(data_path, "association_analysis_summary.csv"), 
          row.names = FALSE)

cat("\n\n==================== Summary Tables Saved ====================\n")
cat("Model comparison summary table saved to:", file.path(data_path, "model_comparison_summary.csv"), "\n")
cat("Association analysis summary table saved to:", file.path(data_path, "association_analysis_summary.csv"), "\n")

# Create detailed stratified analysis results table
stratified_analysis_details <- NULL

# Create detailed results table for each cancer type and stratification analysis
for(i in 1:length(cancer_types)) {
  cancer_type <- cancer_types[i]
  cancer_result <- cancer_results[[i]]
  
  # Process results stratified by Response
  response_levels <- names(cancer_result$stratified_by_response)
  for(level in response_levels) {
    if(!is.null(cancer_result$stratified_by_response[[level]])) {
      temp_df <- data.frame(
        CancerType = cancer_type,
        AnalysisType = "Stratified_by_Response",
        StratificationLevel = level,
        ChiSquareStatistic = cancer_result$stratified_by_response[[level]]$chi_square_test$statistic,
        DegreesOfFreedom = cancer_result$stratified_by_response[[level]]$chi_square_test$parameter,
        PValue = cancer_result$stratified_by_response[[level]]$chi_square_test$p.value,
        stringsAsFactors = FALSE
      )
      stratified_analysis_details <- rbind(stratified_analysis_details, temp_df)
    }
  }
  
  # Process results stratified by PreOrPost
  prepost_levels <- names(cancer_result$stratified_by_pre_or_post)
  for(level in prepost_levels) {
    if(!is.null(cancer_result$stratified_by_pre_or_post[[level]])) {
      temp_df <- data.frame(
        CancerType = cancer_type,
        AnalysisType = "Stratified_by_PreOrPost",
        StratificationLevel = level,
        ChiSquareStatistic = cancer_result$stratified_by_pre_or_post[[level]]$chi_square_test$statistic,
        DegreesOfFreedom = cancer_result$stratified_by_pre_or_post[[level]]$chi_square_test$parameter,
        PValue = cancer_result$stratified_by_pre_or_post[[level]]$chi_square_test$p.value,
        stringsAsFactors = FALSE
      )
      stratified_analysis_details <- rbind(stratified_analysis_details, temp_df)
    }
  }
}

# Save stratified analysis detailed results
write.csv(stratified_analysis_details,
          file = file.path(data_path, "stratified_analysis_details.csv"),
          row.names = FALSE)

cat("Stratified analysis detailed results saved to:", file.path(data_path, "stratified_analysis_details.csv"), "\n")

# Final summary information
cat("\n\n==================== Final Summary ====================\n")
cat("All analysis results and summary tables have been saved to:", data_path, "\n")
cat("- Model comparison summary: model_comparison_summary.csv\n")
cat("- Association analysis summary: association_analysis_summary.csv\n")
cat("- Stratified analysis details: stratified_analysis_details.csv\n")
cat("- Three-dimensional contingency table data: high_dimensional_contingency_table.csv\n")
cat("- Each cancer type 3D table: contingency_table_[CancerType].csv\n")
