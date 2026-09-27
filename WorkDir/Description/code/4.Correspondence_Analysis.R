#2025/12/11
#本脚本用于Correspondence analysis对应分析

#input
#   1.contingency_table_ESCC.csv
#   2.contingency_table_PC.csv
#   3.contingency_table_NSCLC.csv
#   4.contingency_table_TNBC1.csv

#output    
#   1.CA_plot_PC.pdf
#   2.CA_plot_NSCLC.pdf
#   3.CA_plot_TNBC1.pdf
#   4.CA_plot_ESCC.pdf
#   5.Subtype_CancerType_Distribution.pdf 不同亚类中来自不同癌症类型的细胞占比



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#2026/1/7修改：将Minor_Cell_Type的亚类结果作为输入，重新跑该脚本代码。结果均保存在Minor_Cell_Type_Input文件夹之下。


#input:
#   1.

#output
#   1.



# 加载必要的库
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ca))

# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data/Minor_Cell_Type_Input", sep = "")
plot_path <- paste(file_path, "/plot/Minor_Cell_Type_Input", sep = "")





# PC分析
print("开始分析PC数据...")
PC_data <- read.csv(paste(data_path, "/contingency_table_PC.csv", sep = ""))

# 转换因子变量
PC_data$PreOrPost <- as.factor(PC_data$PreOrPost)
PC_data$Response <- as.factor(PC_data$Response)
PC_data$Subtype <- as.factor(PC_data$Subtype)

# 如果CancerType不为NSCLC，则在分析中排除Subtype为3或9的数据
PC_data <- PC_data[!PC_data$Subtype %in% c("Macro_MARCO", "Neutro_FCGR3B"), ]
print("在PC数据中排除了Subtype为Macro_MARCO和Neutro_FCGR3B的记录")

# 合并PreOrPost和Response变量
PC_data$PreOrPost_Response <- interaction(PC_data$PreOrPost, PC_data$Response, sep = "_")

# 创建二维列联表
PC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = PC_data)
print("PC二维列联表:")
print(PC_table_2d)

# 检查并删除零频数行(Subtype)
PC_subtype_sums <- margin.table(PC_table_2d, margin = 1)  # 计算每个Subtype的总频数
PC_zero_freq_subtypes <- names(PC_subtype_sums)[PC_subtype_sums == 0]  # 找出频数为0的Subtype

# 如果存在零频数Subtype，则输出信息并从数据中删除
if (length(PC_zero_freq_subtypes) > 0) {
  print(paste("删除PC中频数为0的Subtype:", paste(PC_zero_freq_subtypes, collapse = ", ")))
  PC_data <- PC_data[!PC_data$Subtype %in% PC_zero_freq_subtypes, ]
  
  # 重新创建二维列联表
  PC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = PC_data)
}

# 检查并删除零频数列(PreOrPost_Response)
PC_prepost_response_sums <- margin.table(PC_table_2d, margin = 2)  # 计算每个PreOrPost_Response的总频数
PC_zero_freq_columns <- names(PC_prepost_response_sums)[PC_prepost_response_sums == 0]  # 找出频数为0的列

# 如果存在零频数列，则输出信息并从数据中删除
if (length(PC_zero_freq_columns) > 0) {
  print(paste("删除PC中频数为0的列:", paste(PC_zero_freq_columns, collapse = ", ")))
  PC_data <- PC_data[!PC_data$PreOrPost_Response %in% PC_zero_freq_columns, ]
  
  # 重新创建二维列联表
  PC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = PC_data)
}

# 检查处理后是否还有数据
if (nrow(PC_data) == 0 || nrow(PC_table_2d) == 0 || sum(PC_table_2d) == 0) {
  stop("错误: PC数据中所有Subtype频数都为0，无法进行对应分析")
}

# 再次检查行和列是否还有零频数
PC_subtype_sums <- margin.table(PC_table_2d, margin = 1)
PC_prepost_response_sums <- margin.table(PC_table_2d, margin = 2)
PC_zero_freq_subtypes <- names(PC_subtype_sums)[PC_subtype_sums == 0]
PC_zero_freq_columns <- names(PC_prepost_response_sums)[PC_prepost_response_sums == 0]

# 清理零频数行
if (length(PC_zero_freq_subtypes) > 0) {
  print(paste("清理PC中频数为0的Subtype:", paste(PC_zero_freq_subtypes, collapse = ", ")))
  PC_table_2d <- PC_table_2d[!rownames(PC_table_2d) %in% PC_zero_freq_subtypes, , drop = FALSE]
}

# 清理零频数列
if (length(PC_zero_freq_columns) > 0) {
  print(paste("清理PC中频数为0的列:", paste(PC_zero_freq_columns, collapse = ", ")))
  PC_table_2d <- PC_table_2d[, !colnames(PC_table_2d) %in% PC_zero_freq_columns, drop = FALSE]
}

# 最终检查
if (nrow(PC_table_2d) == 0 || ncol(PC_table_2d) == 0 || sum(PC_table_2d) == 0) {
  stop("错误: PC数据清理后不满足对应分析要求，无法进行分析")
}

# 进行对应分析
PC_ca_result <- ca(PC_table_2d)

# 创建对应分析图并保存为PDF
PC_plot_title <- "Correspondence Analysis - PC"

# 提取行和列的坐标
PC_row_coords <- data.frame(
  x = PC_ca_result$rowcoord[,1],
  y = PC_ca_result$rowcoord[,2],
  label = rownames(PC_ca_result$rowcoord),
  type = "Row"
)

PC_col_coords <- data.frame(
  x = PC_ca_result$colcoord[,1],
  y = PC_ca_result$colcoord[,2],
  label = rownames(PC_ca_result$colcoord),
  type = "Column"
)

# 合并坐标数据
PC_all_coords <- rbind(PC_row_coords, PC_col_coords)

# 计算坐标轴限制，确保横纵坐标尺度一致
PC_x_range <- range(PC_all_coords$x)
PC_y_range <- range(PC_all_coords$y)
PC_axis_limit <- c(min(PC_x_range[1], PC_y_range[1]), max(PC_x_range[2], PC_y_range[2]))

# 创建ggplot图形
PC_ca_plot <- ggplot(PC_all_coords, aes(x = x, y = y, color = type)) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -1, size = 5) +
  coord_equal(xlim = PC_axis_limit, ylim = PC_axis_limit) +  # 保持坐标轴比例一致
  labs(
#    title = PC_plot_title,
    x = paste("Dimension 1 (Inertia:", round(PC_ca_result$sv[1]^2/sum(PC_ca_result$sv^2)*100, 2), "%)"),
    y = paste("Dimension 2 (Inertia:", round(PC_ca_result$sv[2]^2/sum(PC_ca_result$sv^2)*100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
#    panel.grid = element_blank(),  # 删除背景网格线
    axis.text = element_text(size = 12),  # 确保坐标轴刻度文本可见
    axis.title = element_text(size = 14)  # 确保坐标轴标题可见
  )

# 保存为PDF文件
PC_pdf_file_name <- "CA_plot_PC.pdf"
PC_pdf_file_path <- paste(plot_path, "/", PC_pdf_file_name, sep = "")
ggsave(PC_pdf_file_path, plot = PC_ca_plot, width = 10, height = 10, device = "pdf")
print(paste("已保存对应分析图到:", PC_pdf_file_path))

# 显示图形
print(PC_ca_plot)
print(PC_ca_result)





# ESCC分析
print("开始分析ESCC数据...")
ESCC_data <- read.csv(paste(data_path, "/contingency_table_ESCC.csv", sep = ""))

# 转换因子变量
ESCC_data$PreOrPost <- as.factor(ESCC_data$PreOrPost)
ESCC_data$Response <- as.factor(ESCC_data$Response)
ESCC_data$Subtype <- as.factor(ESCC_data$Subtype)

# 如果CancerType不为NSCLC，则在分析中排除Subtype为3或9的数据
ESCC_data <- ESCC_data[!ESCC_data$Subtype %in% c("Macro_MARCO", "Neutro_FCGR3B"), ]
print("在ESCC数据中排除了Subtype为Macro_MARCO和Neutro_FCGR3B的记录")

# 合并PreOrPost和Response变量
ESCC_data$PreOrPost_Response <- interaction(ESCC_data$PreOrPost, ESCC_data$Response, sep = "_")

# 创建二维列联表
ESCC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = ESCC_data)
print("ESCC二维列联表:")
print(ESCC_table_2d)

# 检查并删除零频数行(Subtype)
ESCC_subtype_sums <- margin.table(ESCC_table_2d, margin = 1)  # 计算每个Subtype的总频数
ESCC_zero_freq_subtypes <- names(ESCC_subtype_sums)[ESCC_subtype_sums == 0]  # 找出频数为0的Subtype

# 如果存在零频数Subtype，则输出信息并从数据中删除
if (length(ESCC_zero_freq_subtypes) > 0) {
  print(paste("删除ESCC中频数为0的Subtype:", paste(ESCC_zero_freq_subtypes, collapse = ", ")))
  ESCC_data <- ESCC_data[!ESCC_data$Subtype %in% ESCC_zero_freq_subtypes, ]
  
  # 重新创建二维列联表
  ESCC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = ESCC_data)
}

# 检查并删除零频数列(PreOrPost_Response)
ESCC_prepost_response_sums <- margin.table(ESCC_table_2d, margin = 2)  # 计算每个PreOrPost_Response的总频数
ESCC_zero_freq_columns <- names(ESCC_prepost_response_sums)[ESCC_prepost_response_sums == 0]  # 找出频数为0的列

# 如果存在零频数列，则输出信息并从数据中删除
if (length(ESCC_zero_freq_columns) > 0) {
  print(paste("删除ESCC中频数为0的列:", paste(ESCC_zero_freq_columns, collapse = ", ")))
  ESCC_data <- ESCC_data[!ESCC_data$PreOrPost_Response %in% ESCC_zero_freq_columns, ]
  
  # 重新创建二维列联表
  ESCC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = ESCC_data)
}

# 检查处理后是否还有数据
if (nrow(ESCC_data) == 0 || nrow(ESCC_table_2d) == 0 || sum(ESCC_table_2d) == 0) {
  stop("错误: ESCC数据中所有Subtype频数都为0，无法进行对应分析")
}

# 再次检查行和列是否还有零频数
ESCC_subtype_sums <- margin.table(ESCC_table_2d, margin = 1)
ESCC_prepost_response_sums <- margin.table(ESCC_table_2d, margin = 2)
ESCC_zero_freq_subtypes <- names(ESCC_subtype_sums)[ESCC_subtype_sums == 0]
ESCC_zero_freq_columns <- names(ESCC_prepost_response_sums)[ESCC_prepost_response_sums == 0]

# 清理零频数行
if (length(ESCC_zero_freq_subtypes) > 0) {
  print(paste("清理ESCC中频数为0的Subtype:", paste(ESCC_zero_freq_subtypes, collapse = ", ")))
  ESCC_table_2d <- ESCC_table_2d[!rownames(ESCC_table_2d) %in% ESCC_zero_freq_subtypes, , drop = FALSE]
}

# 清理零频数列
if (length(ESCC_zero_freq_columns) > 0) {
  print(paste("清理ESCC中频数为0的列:", paste(ESCC_zero_freq_columns, collapse = ", ")))
  ESCC_table_2d <- ESCC_table_2d[, !colnames(ESCC_table_2d) %in% ESCC_zero_freq_columns, drop = FALSE]
}

# 最终检查
if (nrow(ESCC_table_2d) == 0 || ncol(ESCC_table_2d) == 0 || sum(ESCC_table_2d) == 0) {
  stop("错误: ESCC数据清理后不满足对应分析要求，无法进行分析")
}

# 进行对应分析
ESCC_ca_result <- ca(ESCC_table_2d)

# 创建对应分析图并保存为PDF
ESCC_plot_title <- "Correspondence Analysis - ESCC"

# 提取行和列的坐标
ESCC_row_coords <- data.frame(
  x = ESCC_ca_result$rowcoord[,1],
  y = ESCC_ca_result$rowcoord[,2],
  label = rownames(ESCC_ca_result$rowcoord),
  type = "Row"
)

ESCC_col_coords <- data.frame(
  x = ESCC_ca_result$colcoord[,1],
  y = ESCC_ca_result$colcoord[,2],
  label = rownames(ESCC_ca_result$colcoord),
  type = "Column"
)

# 合并坐标数据
ESCC_all_coords <- rbind(ESCC_row_coords, ESCC_col_coords)

# 计算坐标轴限制，确保横纵坐标尺度一致
ESCC_x_range <- range(ESCC_all_coords$x)
ESCC_y_range <- range(ESCC_all_coords$y)
ESCC_axis_limit <- c(min(ESCC_x_range[1], ESCC_y_range[1]), max(ESCC_x_range[2], ESCC_y_range[2]))

# 创建ggplot图形
ESCC_ca_plot <- ggplot(ESCC_all_coords, aes(x = x, y = y, color = type)) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -1, size = 5) +
  coord_equal(xlim = ESCC_axis_limit, ylim = ESCC_axis_limit) +  # 保持坐标轴比例一致
  labs(
#    title = ESCC_plot_title,
    x = paste("Dimension 1 (Inertia:", round(ESCC_ca_result$sv[1]^2/sum(ESCC_ca_result$sv^2)*100, 2), "%)"),
    y = paste("Dimension 2 (Inertia:", round(ESCC_ca_result$sv[2]^2/sum(ESCC_ca_result$sv^2)*100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
#    panel.grid = element_blank(),  # 删除背景网格线
    axis.text = element_text(size = 12),  # 确保坐标轴刻度文本可见
    axis.title = element_text(size = 14)  # 确保坐标轴标题可见
  )

# 保存为PDF文件
ESCC_pdf_file_name <- "CA_plot_ESCC.pdf"
ESCC_pdf_file_path <- paste(plot_path, "/", ESCC_pdf_file_name, sep = "")
ggsave(ESCC_pdf_file_path, plot = ESCC_ca_plot, width = 10, height = 10, device = "pdf")
print(paste("已保存对应分析图到:", ESCC_pdf_file_path))

# 显示图形
print(ESCC_ca_plot)
print(ESCC_ca_result)





# NSCLC分析
print("开始分析NSCLC数据...")
NSCLC_data <- read.csv(paste(data_path, "/contingency_table_NSCLC.csv", sep = ""))

# 转换因子变量
NSCLC_data$PreOrPost <- as.factor(NSCLC_data$PreOrPost)
NSCLC_data$Response <- as.factor(NSCLC_data$Response)
NSCLC_data$Subtype <- as.factor(NSCLC_data$Subtype)



# 合并PreOrPost和Response变量
NSCLC_data$PreOrPost_Response <- interaction(NSCLC_data$PreOrPost, NSCLC_data$Response, sep = "_")

# 创建二维列联表
NSCLC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = NSCLC_data)
print("NSCLC二维列联表:")
print(NSCLC_table_2d)

# 检查并删除零频数行(Subtype)
NSCLC_subtype_sums <- margin.table(NSCLC_table_2d, margin = 1)  # 计算每个Subtype的总频数
NSCLC_zero_freq_subtypes <- names(NSCLC_subtype_sums)[NSCLC_subtype_sums == 0]  # 找出频数为0的Subtype

# 如果存在零频数Subtype，则输出信息并从数据中删除
if (length(NSCLC_zero_freq_subtypes) > 0) {
  print(paste("删除NSCLC中频数为0的Subtype:", paste(NSCLC_zero_freq_subtypes, collapse = ", ")))
  NSCLC_data <- NSCLC_data[!NSCLC_data$Subtype %in% NSCLC_zero_freq_subtypes, ]
  
  # 重新创建二维列联表
  NSCLC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = NSCLC_data)
}

# 检查并删除零频数列(PreOrPost_Response)
NSCLC_prepost_response_sums <- margin.table(NSCLC_table_2d, margin = 2)  # 计算每个PreOrPost_Response的总频数
NSCLC_zero_freq_columns <- names(NSCLC_prepost_response_sums)[NSCLC_prepost_response_sums == 0]  # 找出频数为0的列

# 如果存在零频数列，则输出信息并从数据中删除
if (length(NSCLC_zero_freq_columns) > 0) {
  print(paste("删除NSCLC中频数为0的列:", paste(NSCLC_zero_freq_columns, collapse = ", ")))
  NSCLC_data <- NSCLC_data[!NSCLC_data$PreOrPost_Response %in% NSCLC_zero_freq_columns, ]
  
  # 重新创建二维列联表
  NSCLC_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = NSCLC_data)
}

# 检查处理后是否还有数据
if (nrow(NSCLC_data) == 0 || nrow(NSCLC_table_2d) == 0 || sum(NSCLC_table_2d) == 0) {
  stop("错误: NSCLC数据中所有Subtype频数都为0，无法进行对应分析")
}

# 再次检查行和列是否还有零频数
NSCLC_subtype_sums <- margin.table(NSCLC_table_2d, margin = 1)
NSCLC_prepost_response_sums <- margin.table(NSCLC_table_2d, margin = 2)
NSCLC_zero_freq_subtypes <- names(NSCLC_subtype_sums)[NSCLC_subtype_sums == 0]
NSCLC_zero_freq_columns <- names(NSCLC_prepost_response_sums)[NSCLC_prepost_response_sums == 0]

# 清理零频数行
if (length(NSCLC_zero_freq_subtypes) > 0) {
  print(paste("清理NSCLC中频数为0的Subtype:", paste(NSCLC_zero_freq_subtypes, collapse = ", ")))
  NSCLC_table_2d <- NSCLC_table_2d[!rownames(NSCLC_table_2d) %in% NSCLC_zero_freq_subtypes, , drop = FALSE]
}

# 清理零频数列
if (length(NSCLC_zero_freq_columns) > 0) {
  print(paste("清理NSCLC中频数为0的列:", paste(NSCLC_zero_freq_columns, collapse = ", ")))
  NSCLC_table_2d <- NSCLC_table_2d[, !colnames(NSCLC_table_2d) %in% NSCLC_zero_freq_columns, drop = FALSE]
}

# 最终检查
if (nrow(NSCLC_table_2d) == 0 || ncol(NSCLC_table_2d) == 0 || sum(NSCLC_table_2d) == 0) {
  stop("错误: NSCLC数据清理后不满足对应分析要求，无法进行分析")
}

# 进行对应分析
NSCLC_ca_result <- ca(NSCLC_table_2d)

# 创建对应分析图并保存为PDF
NSCLC_plot_title <- "Correspondence Analysis - NSCLC"

# 提取行和列的坐标
NSCLC_row_coords <- data.frame(
  x = NSCLC_ca_result$rowcoord[,1],
  y = NSCLC_ca_result$rowcoord[,2],
  label = rownames(NSCLC_ca_result$rowcoord),
  type = "Row"
)

NSCLC_col_coords <- data.frame(
  x = NSCLC_ca_result$colcoord[,1],
  y = NSCLC_ca_result$colcoord[,2],
  label = rownames(NSCLC_ca_result$colcoord),
  type = "Column"
)

# 合并坐标数据
NSCLC_all_coords <- rbind(NSCLC_row_coords, NSCLC_col_coords)

# 计算坐标轴限制，确保横纵坐标尺度一致
NSCLC_x_range <- range(NSCLC_all_coords$x)
NSCLC_y_range <- range(NSCLC_all_coords$y)
NSCLC_axis_limit <- c(min(NSCLC_x_range[1], NSCLC_y_range[1]), max(NSCLC_x_range[2], NSCLC_y_range[2]))

# 创建ggplot图形
NSCLC_ca_plot <- ggplot(NSCLC_all_coords, aes(x = x, y = y, color = type)) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -1, size = 5) +
  coord_equal(xlim = NSCLC_axis_limit, ylim = NSCLC_axis_limit) +  # 保持坐标轴比例一致
  labs(
#    title = NSCLC_plot_title,
    x = paste("Dimension 1 (Inertia:", round(NSCLC_ca_result$sv[1]^2/sum(NSCLC_ca_result$sv^2)*100, 2), "%)"),
    y = paste("Dimension 2 (Inertia:", round(NSCLC_ca_result$sv[2]^2/sum(NSCLC_ca_result$sv^2)*100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
#    panel.grid = element_blank(),  # 删除背景网格线
    axis.text = element_text(size = 12),  # 确保坐标轴刻度文本可见
    axis.title = element_text(size = 14)  # 确保坐标轴标题可见
  )

# 保存为PDF文件
NSCLC_pdf_file_name <- "CA_plot_NSCLC.pdf"
NSCLC_pdf_file_path <- paste(plot_path, "/", NSCLC_pdf_file_name, sep = "")
ggsave(NSCLC_pdf_file_path, plot = NSCLC_ca_plot, width = 10, height = 10, device = "pdf")
print(paste("已保存对应分析图到:", NSCLC_pdf_file_path))

# 显示图形
print(NSCLC_ca_plot)
print(NSCLC_ca_result)





# TNBC1分析
print("开始分析TNBC1数据...")
TNBC1_data <- read.csv(paste(data_path, "/contingency_table_TNBC1.csv", sep = ""))

# 转换因子变量
TNBC1_data$PreOrPost <- as.factor(TNBC1_data$PreOrPost)
TNBC1_data$Response <- as.factor(TNBC1_data$Response)
TNBC1_data$Subtype <- as.factor(TNBC1_data$Subtype)

# 如果CancerType不为NSCLC，则在分析中排除Subtype为3或9的数据
TNBC1_data <- TNBC1_data[!TNBC1_data$Subtype %in% c("Macro_MARCO", "Neutro_FCGR3B"), ]
print("在TNBC1数据中排除了Subtype为Macro_MARCO和Neutro_FCGR3B的记录")

# 对于TNBC1，使用并列形式的列属性而不是组合形式
# 创建三维列联表
TNBC1_table_3d <- xtabs(Freq ~ Subtype + PreOrPost + Response, data = TNBC1_data)
print("TNBC1三维列联表:")
print(TNBC1_table_3d)

# 将三维表转换为二维表用于对应分析
TNBC1_data$PreOrPost_Response <- interaction(TNBC1_data$PreOrPost, TNBC1_data$Response, sep = "_")
TNBC1_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = TNBC1_data)

print("TNBC1二维列联表:")
print(TNBC1_table_2d)

# 检查并删除零频数行(Subtype)
TNBC1_subtype_sums <- margin.table(TNBC1_table_2d, margin = 1)  # 计算每个Subtype的总频数
TNBC1_zero_freq_subtypes <- names(TNBC1_subtype_sums)[TNBC1_subtype_sums == 0]  # 找出频数为0的Subtype

# 如果存在零频数Subtype，则输出信息并从数据中删除
if (length(TNBC1_zero_freq_subtypes) > 0) {
  print(paste("删除TNBC1中频数为0的Subtype:", paste(TNBC1_zero_freq_subtypes, collapse = ", ")))
  TNBC1_data <- TNBC1_data[!TNBC1_data$Subtype %in% TNBC1_zero_freq_subtypes, ]
  
  # 重新创建二维列联表
  TNBC1_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = TNBC1_data)
}

# 检查并删除零频数列(PreOrPost_Response)
TNBC1_prepost_response_sums <- margin.table(TNBC1_table_2d, margin = 2)  # 计算每个PreOrPost_Response的总频数
TNBC1_zero_freq_columns <- names(TNBC1_prepost_response_sums)[TNBC1_prepost_response_sums == 0]  # 找出频数为0的列

# 如果存在零频数列，则输出信息并从数据中删除
if (length(TNBC1_zero_freq_columns) > 0) {
  print(paste("删除TNBC1中频数为0的列:", paste(TNBC1_zero_freq_columns, collapse = ", ")))
  TNBC1_data <- TNBC1_data[!TNBC1_data$PreOrPost_Response %in% TNBC1_zero_freq_columns, ]
  
  # 重新创建二维列联表
  TNBC1_table_2d <- xtabs(Freq ~ Subtype + PreOrPost_Response, data = TNBC1_data)
}

# 检查处理后是否还有数据
if (nrow(TNBC1_data) == 0 || nrow(TNBC1_table_2d) == 0 || sum(TNBC1_table_2d) == 0) {
  stop("错误: TNBC1数据中所有Subtype频数都为0，无法进行对应分析")
}

# 再次检查行和列是否还有零频数
TNBC1_subtype_sums <- margin.table(TNBC1_table_2d, margin = 1)
TNBC1_prepost_response_sums <- margin.table(TNBC1_table_2d, margin = 2)
TNBC1_zero_freq_subtypes <- names(TNBC1_subtype_sums)[TNBC1_subtype_sums == 0]
TNBC1_zero_freq_columns <- names(TNBC1_prepost_response_sums)[TNBC1_prepost_response_sums == 0]

# 清理零频数行
if (length(TNBC1_zero_freq_subtypes) > 0) {
  print(paste("清理TNBC1中频数为0的Subtype:", paste(TNBC1_zero_freq_subtypes, collapse = ", ")))
  TNBC1_table_2d <- TNBC1_table_2d[!rownames(TNBC1_table_2d) %in% TNBC1_zero_freq_subtypes, , drop = FALSE]
}

# 清理零频数列
if (length(TNBC1_zero_freq_columns) > 0) {
  print(paste("清理TNBC1中频数为0的列:", paste(TNBC1_zero_freq_columns, collapse = ", ")))
  TNBC1_table_2d <- TNBC1_table_2d[, !colnames(TNBC1_table_2d) %in% TNBC1_zero_freq_columns, drop = FALSE]
}

# 最终检查
if (nrow(TNBC1_table_2d) == 0 || ncol(TNBC1_table_2d) == 0 || sum(TNBC1_table_2d) == 0) {
  stop("错误: TNBC1数据清理后不满足对应分析要求，无法进行分析")
}

# 进行对应分析
TNBC1_ca_result <- ca(TNBC1_table_2d)

# 创建对应分析图并保存为PDF
TNBC1_plot_title <- "Correspondence Analysis - TNBC1"

# 提取行和列的坐标
TNBC1_row_coords <- data.frame(
  x = TNBC1_ca_result$rowcoord[,1],
  y = TNBC1_ca_result$rowcoord[,2],
  label = rownames(TNBC1_ca_result$rowcoord),
  type = "Row"
)

TNBC1_col_coords <- data.frame(
  x = TNBC1_ca_result$colcoord[,1],
  y = TNBC1_ca_result$colcoord[,2],
  label = rownames(TNBC1_ca_result$colcoord),
  type = "Column"
)

# 合并坐标数据
TNBC1_all_coords <- rbind(TNBC1_row_coords, TNBC1_col_coords)

# 计算坐标轴限制，确保横纵坐标尺度一致
TNBC1_x_range <- range(TNBC1_all_coords$x)
TNBC1_y_range <- range(TNBC1_all_coords$y)
TNBC1_axis_limit <- c(min(TNBC1_x_range[1], TNBC1_y_range[1]), max(TNBC1_x_range[2], TNBC1_y_range[2]))

# 创建ggplot图形
TNBC1_ca_plot <- ggplot(TNBC1_all_coords, aes(x = x, y = y, color = type)) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -1, size = 5) +
  coord_equal(xlim = TNBC1_axis_limit, ylim = TNBC1_axis_limit) +  # 保持坐标轴比例一致
  labs(
#    title = TNBC1_plot_title,
    x = paste("Dimension 1 (Inertia:", round(TNBC1_ca_result$sv[1]^2/sum(TNBC1_ca_result$sv^2)*100, 2), "%)"),
    y = paste("Dimension 2 (Inertia:", round(TNBC1_ca_result$sv[2]^2/sum(TNBC1_ca_result$sv^2)*100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
#    panel.grid = element_blank(),  # 删除背景网格线
    axis.text = element_text(size = 12),  # 确保坐标轴刻度文本可见
    axis.title = element_text(size = 14)  # 确保坐标轴标题可见
  )

# 保存为PDF文件
TNBC1_pdf_file_name <- "CA_plot_TNBC1.pdf"
TNBC1_pdf_file_path <- paste(plot_path, "/", TNBC1_pdf_file_name, sep = "")
ggsave(TNBC1_pdf_file_path, plot = TNBC1_ca_plot, width = 10, height = 10, device = "pdf")
print(paste("已保存对应分析图到:", TNBC1_pdf_file_path))

# 显示图形
print(TNBC1_ca_plot)
print(TNBC1_ca_result)

# 新增功能：汇总不同CancerType的细胞，计算不同Subtype下各CancerType所占比例，并输出为堆积柱状图
# 读取所有癌症类型的数据并合并
all_data <- data.frame()

cancer_types <- c("PC", "ESCC", "NSCLC", "TNBC1")
for (cancer_type in cancer_types) {
  file_name <- paste("contingency_table_", cancer_type, ".csv", sep = "")
  file_path <- paste(data_path, "/", file_name, sep = "")
  data <- read.csv(file_path)
  
  # 添加癌症类型列
  data$CancerType <- cancer_type
  
  
  # 合并数据
  all_data <- rbind(all_data, data)
}

# 转换因子变量
all_data$PreOrPost <- as.factor(all_data$PreOrPost)
all_data$Response <- as.factor(all_data$Response)
all_data$Subtype <- as.factor(all_data$Subtype)
all_data$CancerType <- as.factor(all_data$CancerType)

# 创建按Subtype和CancerType汇总的列联表
subtype_cancer_table <- xtabs(Freq ~ Subtype + CancerType, data = all_data)

# 计算每个Subtype中各CancerType的比例
subtype_cancer_prop <- prop.table(subtype_cancer_table, margin = 1)

# 转换为数据框格式用于绘图
subtype_cancer_df <- as.data.frame(subtype_cancer_prop)
colnames(subtype_cancer_df) <- c("Subtype", "CancerType", "Proportion")

# 创建堆积柱状图
stacked_bar_plot <- ggplot(subtype_cancer_df, aes(x = Subtype, y = Proportion, fill = CancerType)) +
  geom_col(position = "stack") +
  labs(
#    title = "Distribution of Cancer Types within Each Subtype",
    x = "Subtype",
    y = "Proportion",
    fill = "Cancer Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# 保存堆积柱状图
stacked_bar_pdf_path <- paste(plot_path, "/Subtype_CancerType_Distribution.pdf", sep = "")
ggsave(stacked_bar_pdf_path, plot = stacked_bar_plot, width = 12, height = 8, device = "pdf")
print(paste("已保存Subtype中CancerType分布的堆积柱状图到:", stacked_bar_pdf_path))

# 显示图形
print(stacked_bar_plot)





#解释结果：
#Subtype——Macro_5和Neutro_1中，来自NSCLC的细胞占比超过85%，故认为该两种Subtype是NSCLC专属的亚型（组织驻留型）



#总结所有结果：
#PC：列变量的四个取值互相分离，Pre_Partial有较为独立的Neighborhood（和DC_1,DC_2,Mast）
#ESCC：列变量有六个取值，彼此离得比较近，不易互相区分。Pre_pCR和Macro_1高度接近
#NSCLC：列变量四个取值，Post_pCR和Post_MPR离得比较近，与其他两个变量离得比较远。
#       Pre_NMPR周围无Subtype编号（注意：可能是由于本研究数据未匹配导致，Pre中只有两例NMPR的结果，缺少关键对比）
#       Post_NMPR和Macro_2高度接近
#TNBC1：列变量四个取值，彼此分离。Pre_Partial和Mast高度接近，Pre_Stable和Macro_1高度接近



#由于组织驻留型细胞的存在，怀疑Macro_5和Neutro_1两类细胞是NSCLC专属存在的亚型（在各自亚型中来自NSCLC的细胞占比超过90%）
#故在分析其他三种癌症类型的时候应该删除这两类






#
# 新增功能：为指定亚类绘制比例柱状图并添加折线图


# 分别为每种癌症类型创建图表

# 为ESCC数据中的Macro_MS4A6A和Macro_APOE亚类创建比例柱状图并添加折线图
print("正在为ESCC中的Macro_MS4A6A和Macro_APOE亚类绘制比例柱状图...")
ESCC_table_2d <- get("ESCC_table_2d")
target_subtypes_ESCC <- c("Macro_MS4A6A", "Macro_APOE")

# 筛选出目标亚类的数据
selected_rows_ESCC <- rownames(ESCC_table_2d)[rownames(ESCC_table_2d) %in% target_subtypes_ESCC]
subtype_table_ESCC <- ESCC_table_2d[selected_rows_ESCC, , drop = FALSE]

if(length(selected_rows_ESCC) > 0){
  # 将表格转换为长格式数据框
  plot_data_ESCC <- as.data.frame(as.table(subtype_table_ESCC))
  colnames(plot_data_ESCC) <- c("Subtype", "TreatmentStatus", "Freq")
  plot_data_ESCC$Freq <- as.numeric(plot_data_ESCC$Freq)
  
  # 计算每个治疗状态下所有亚类的总和
  total_by_treatmentstatus <- as.data.frame(as.table(ESCC_table_2d))
  colnames(total_by_treatmentstatus) <- c("Subtype", "TreatmentStatus", "TotalFreq")
  total_by_treatmentstatus$TotalFreq <- as.numeric(total_by_treatmentstatus$TotalFreq)
  total_by_treatmentstatus <- aggregate(TotalFreq ~ TreatmentStatus, data = total_by_treatmentstatus, FUN = sum)
  
  # 合并数据以计算比例
  plot_data_ESCC <- merge(plot_data_ESCC, total_by_treatmentstatus, by = "TreatmentStatus")
  plot_data_ESCC$Prop <- plot_data_ESCC$Freq / plot_data_ESCC$TotalFreq
  
  # 移除Prop为0的行
  plot_data_ESCC <- plot_data_ESCC[plot_data_ESCC$Prop > 0, ]
  
  # 按Subtype和TreatmentStatus排序，确保折线图能正确连接点
  plot_data_ESCC <- plot_data_ESCC[order(plot_data_ESCC$Subtype, plot_data_ESCC$TreatmentStatus),]
  
  # 为亚类指定特定颜色
  subtype_colors <- c(
    "Macro_MS4A6A" = "#A6CEE3",
    "Macro_APOE" = "#FB9A99"
  )
  
  # 为当前数据中实际存在的亚类选择颜色
  available_colors <- subtype_colors[unique(plot_data_ESCC$Subtype)]
  
  p_ESCC <- ggplot(plot_data_ESCC, aes(x = TreatmentStatus, y = Prop, fill = Subtype, color = Subtype)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.7, width = 0.6) +
    geom_line(aes(group = Subtype), position = position_dodge(width = 0.8), size = 1) +
    geom_point(aes(group = Subtype), position = position_dodge(width = 0.8), size = 2) +
    scale_fill_manual(values = available_colors) +
    scale_color_manual(values = available_colors) +
    labs(
#      title = paste("Proportion Distribution for Macro_MS4A6A and Macro_APOE in ESCC"),
      x = "Treatment Status (PreOrPost_Response)",
      y = "Proportion",
      fill = "Subtype",
      color = "Subtype",
      caption = "Bar chart with overlay line graph for two subtypes in ESCC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.05))
  
  # 保存图表
  pdf_filename_ESCC_detailed <- "Prop_Bar_Line_Plot_ESCC_Macro_MS4A6A_APOE.pdf"
  pdf_filepath_ESCC_detailed <- paste(plot_path, "/", pdf_filename_ESCC_detailed, sep = "")
  ggsave(pdf_filepath_ESCC_detailed, plot = p_ESCC, width = 12, height = 6, device = "pdf")
  print(paste("已保存ESCC中Macro_MS4A6A和Macro_APOE亚类的比例柱状图(带折线):", pdf_filepath_ESCC_detailed))
  
  # 显示图形
  print(p_ESCC)
}

# 为NSCLC数据中的Macro_APOE亚类创建比例柱状图并添加折线图
print("正在为NSCLC中的Macro_APOE亚类绘制比例柱状图...")
NSCLC_table_2d <- get("NSCLC_table_2d")
target_subtypes_NSCLC <- c("Macro_APOE")

# 筛选出目标亚类的数据
selected_rows_NSCLC <- rownames(NSCLC_table_2d)[rownames(NSCLC_table_2d) %in% target_subtypes_NSCLC]
subtype_table_NSCLC <- NSCLC_table_2d[selected_rows_NSCLC, , drop = FALSE]

# 如果是NSCLC，过滤掉Pre_MPR和Pre_pCR列
cols_to_keep_NSCLC <- colnames(NSCLC_table_2d)[!grepl("Pre_MPR|Pre_pCR", colnames(NSCLC_table_2d))]
subtype_table_NSCLC <- subtype_table_NSCLC[, cols_to_keep_NSCLC, drop = FALSE]
NSCLC_table_2d_filtered <- NSCLC_table_2d[, cols_to_keep_NSCLC, drop = FALSE]

if(length(selected_rows_NSCLC) > 0){
  # 将表格转换为长格式数据框
  plot_data_NSCLC <- as.data.frame(as.table(subtype_table_NSCLC))
  colnames(plot_data_NSCLC) <- c("Subtype", "TreatmentStatus", "Freq")
  plot_data_NSCLC$Freq <- as.numeric(plot_data_NSCLC$Freq)
  
  # 计算每个治疗状态下所有亚类的总和
  total_by_treatmentstatus_NSCLC <- as.data.frame(as.table(NSCLC_table_2d_filtered))
  colnames(total_by_treatmentstatus_NSCLC) <- c("Subtype", "TreatmentStatus", "TotalFreq")
  total_by_treatmentstatus_NSCLC$TotalFreq <- as.numeric(total_by_treatmentstatus_NSCLC$TotalFreq)
  total_by_treatmentstatus_NSCLC <- aggregate(TotalFreq ~ TreatmentStatus, data = total_by_treatmentstatus_NSCLC, FUN = sum)
  
  # 合并数据以计算比例
  plot_data_NSCLC <- merge(plot_data_NSCLC, total_by_treatmentstatus_NSCLC, by = "TreatmentStatus")
  plot_data_NSCLC$Prop <- plot_data_NSCLC$Freq / plot_data_NSCLC$TotalFreq
  
  # 移除Prop为0的行
  plot_data_NSCLC <- plot_data_NSCLC[plot_data_NSCLC$Prop > 0, ]
  
  # 按Subtype和TreatmentStatus排序，确保折线图能正确连接点
  plot_data_NSCLC <- plot_data_NSCLC[order(plot_data_NSCLC$Subtype, plot_data_NSCLC$TreatmentStatus),]
  
  # 为亚类指定特定颜色
  subtype_colors_NSCLC <- c(
    "Macro_APOE" = "#FB9A99"
  )
  
  # 为当前数据中实际存在的亚类选择颜色
  available_colors_NSCLC <- subtype_colors_NSCLC[unique(plot_data_NSCLC$Subtype)]
  
  p_NSCLC <- ggplot(plot_data_NSCLC, aes(x = TreatmentStatus, y = Prop, fill = Subtype, color = Subtype)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.7, width = 0.6) +
    geom_line(aes(group = Subtype), position = position_dodge(width = 0.8), size = 1) +
    geom_point(aes(group = Subtype), position = position_dodge(width = 0.8), size = 2) +
    scale_fill_manual(values = available_colors_NSCLC) +
    scale_color_manual(values = available_colors_NSCLC) +
    labs(
#      title = paste("Proportion Distribution for Macro_APOE in NSCLC"),
      x = "Treatment Status (PreOrPost_Response)",
      y = "Proportion",
      fill = "Subtype",
      color = "Subtype",
      caption = "Bar chart with overlay line graph for Macro_APOE in NSCLC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.05))
  
  # 保存图表
  pdf_filename_NSCLC_detailed <- "Prop_Bar_Line_Plot_NSCLC_Macro_APOE.pdf"
  pdf_filepath_NSCLC_detailed <- paste(plot_path, "/", pdf_filename_NSCLC_detailed, sep = "")
  ggsave(pdf_filepath_NSCLC_detailed, plot = p_NSCLC, width = 12, height = 6, device = "pdf")
  print(paste("已保存NSCLC中Macro_APOE亚类的比例柱状图(带折线):", pdf_filepath_NSCLC_detailed))
  
  # 显示图形
  print(p_NSCLC)
}

# 为TNBC1数据中的Mast和Macro_MS4A6A亚类创建比例柱状图并添加折线图
print("正在为TNBC1中的Mast和Macro_MS4A6A亚类绘制比例柱状图...")
TNBC1_table_2d <- get("TNBC1_table_2d")
target_subtypes_TNBC1 <- c("Mast", "Macro_MS4A6A")

# 筛选出目标亚类的数据
selected_rows_TNBC1 <- rownames(TNBC1_table_2d)[rownames(TNBC1_table_2d) %in% target_subtypes_TNBC1]
subtype_table_TNBC1 <- TNBC1_table_2d[selected_rows_TNBC1, , drop = FALSE]

if(length(selected_rows_TNBC1) > 0){
  # 将表格转换为长格式数据框
  plot_data_TNBC1 <- as.data.frame(as.table(subtype_table_TNBC1))
  colnames(plot_data_TNBC1) <- c("Subtype", "TreatmentStatus", "Freq")
  plot_data_TNBC1$Freq <- as.numeric(plot_data_TNBC1$Freq)
  
  # 计算每个治疗状态下所有亚类的总和
  total_by_treatmentstatus_TNBC1 <- as.data.frame(as.table(TNBC1_table_2d))
  colnames(total_by_treatmentstatus_TNBC1) <- c("Subtype", "TreatmentStatus", "TotalFreq")
  total_by_treatmentstatus_TNBC1$TotalFreq <- as.numeric(total_by_treatmentstatus_TNBC1$TotalFreq)
  total_by_treatmentstatus_TNBC1 <- aggregate(TotalFreq ~ TreatmentStatus, data = total_by_treatmentstatus_TNBC1, FUN = sum)
  
  # 合并数据以计算比例
  plot_data_TNBC1 <- merge(plot_data_TNBC1, total_by_treatmentstatus_TNBC1, by = "TreatmentStatus")
  plot_data_TNBC1$Prop <- plot_data_TNBC1$Freq / plot_data_TNBC1$TotalFreq
  
  # 移除Prop为0的行
  plot_data_TNBC1 <- plot_data_TNBC1[plot_data_TNBC1$Prop > 0, ]
  
  # 按Subtype和TreatmentStatus排序，确保折线图能正确连接点
  plot_data_TNBC1 <- plot_data_TNBC1[order(plot_data_TNBC1$Subtype, plot_data_TNBC1$TreatmentStatus),]
  
  # 为亚类指定特定颜色
  subtype_colors_TNBC1 <- c(
    "Mast" = "#B2DF8A",
    "Macro_MS4A6A" = "#A6CEE3"
  )
  
  # 为当前数据中实际存在的亚类选择颜色
  available_colors_TNBC1 <- subtype_colors_TNBC1[unique(plot_data_TNBC1$Subtype)]
  
  p_TNBC1 <- ggplot(plot_data_TNBC1, aes(x = TreatmentStatus, y = Prop, fill = Subtype, color = Subtype)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.7, width = 0.6) +
    geom_line(aes(group = Subtype), position = position_dodge(width = 0.8), size = 1) +
    geom_point(aes(group = Subtype), position = position_dodge(width = 0.8), size = 2) +
    scale_fill_manual(values = available_colors_TNBC1) +
    scale_color_manual(values = available_colors_TNBC1) +
    labs(
#      title = paste("Proportion Distribution for Mast and Macro_MS4A6A in TNBC"),
      x = "Treatment Status (PreOrPost_Response)",
      y = "Proportion",
      fill = "Subtype",
      color = "Subtype",
      caption = "Bar chart with overlay line graph for two subtypes in TNBC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.05))
  
  # 保存图表
  pdf_filename_TNBC1_detailed <- "Prop_Bar_Line_Plot_TNBC1_Mast_Macro_MS4A6A.pdf"
  pdf_filepath_TNBC1_detailed <- paste(plot_path, "/", pdf_filename_TNBC1_detailed, sep = "")
  ggsave(pdf_filepath_TNBC1_detailed, plot = p_TNBC1, width = 12, height = 6, device = "pdf")
  print(paste("已保存TNBC1中Mast和Macro_MS4A6A亚类的比例柱状图(带折线):", pdf_filepath_TNBC1_detailed))
  
  # 显示图形
  print(p_TNBC1)
}

# 为PC数据中的Mast亚类创建比例柱状图并添加折线图
print("正在为PC中的Mast亚类绘制比例柱状图...")
PC_table_2d <- get("PC_table_2d")
target_subtypes_PC <- c("Mast")

# 筛选出目标亚类的数据
selected_rows_PC <- rownames(PC_table_2d)[rownames(PC_table_2d) %in% target_subtypes_PC]
subtype_table_PC <- PC_table_2d[selected_rows_PC, , drop = FALSE]

if(length(selected_rows_PC) > 0){
  # 将表格转换为长格式数据框
  plot_data_PC <- as.data.frame(as.table(subtype_table_PC))
  colnames(plot_data_PC) <- c("Subtype", "TreatmentStatus", "Freq")
  plot_data_PC$Freq <- as.numeric(plot_data_PC$Freq)
  
  # 计算每个治疗状态下所有亚类的总和
  total_by_treatmentstatus_PC <- as.data.frame(as.table(PC_table_2d))
  colnames(total_by_treatmentstatus_PC) <- c("Subtype", "TreatmentStatus", "TotalFreq")
  total_by_treatmentstatus_PC$TotalFreq <- as.numeric(total_by_treatmentstatus_PC$TotalFreq)
  total_by_treatmentstatus_PC <- aggregate(TotalFreq ~ TreatmentStatus, data = total_by_treatmentstatus_PC, FUN = sum)
  
  # 合并数据以计算比例
  plot_data_PC <- merge(plot_data_PC, total_by_treatmentstatus_PC, by = "TreatmentStatus")
  plot_data_PC$Prop <- plot_data_PC$Freq / plot_data_PC$TotalFreq
  
  # 移除Prop为0的行
  plot_data_PC <- plot_data_PC[plot_data_PC$Prop > 0, ]
  
  # 按Subtype和TreatmentStatus排序，确保折线图能正确连接点
  plot_data_PC <- plot_data_PC[order(plot_data_PC$Subtype, plot_data_PC$TreatmentStatus),]
  
  # 为亚类指定特定颜色
  subtype_colors_PC <- c(
    "Mast" = "#B2DF8A"
  )
  
  # 为当前数据中实际存在的亚类选择颜色
  available_colors_PC <- subtype_colors_PC[unique(plot_data_PC$Subtype)]
  
  p_PC <- ggplot(plot_data_PC, aes(x = TreatmentStatus, y = Prop, fill = Subtype, color = Subtype)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.7, width = 0.6) +
    geom_line(aes(group = Subtype), position = position_dodge(width = 0.8), size = 1) +
    geom_point(aes(group = Subtype), position = position_dodge(width = 0.8), size = 2) +
    scale_fill_manual(values = available_colors_PC) +
    scale_color_manual(values = available_colors_PC) +
    labs(
#      title = paste("Proportion Distribution for Mast in PC"),
      x = "Treatment Status (PreOrPost_Response)",
      y = "Proportion",
      fill = "Subtype",
      color = "Subtype",
      caption = "Bar chart with overlay line graph for Mast in PC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.05))
  
  # 保存图表
  pdf_filename_PC_detailed <- "Prop_Bar_Line_Plot_PC_Mast.pdf"
  pdf_filepath_PC_detailed <- paste(plot_path, "/", pdf_filename_PC_detailed, sep = "")
  ggsave(pdf_filepath_PC_detailed, plot = p_PC, width = 12, height = 6, device = "pdf")
  print(paste("已保存PC中Mast亚类的比例柱状图(带折线):", pdf_filepath_PC_detailed))
  
  # 显示图形
  print(p_PC)
}
