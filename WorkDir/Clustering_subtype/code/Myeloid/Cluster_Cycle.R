#2025/11/11
#本脚本用与循环多种聚类结果，从中挑选有亚群Response.Comprehensive比例显著不同的聚类结果



#input：
#   1.SData_Integrated_rpca.RDS     整合后未聚类的数据

#output：
#   1.SData_Integrated_ClusterCycle.RDS     循环聚类后的RDS文件
#   2.Cluster_Cycle文件夹里的全部图片



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
SData <- readRDS(file.path(data_path, "SData_Integrated_rpca.RDS"))
#DemoCellName <- sample(colnames(SData), size = 5000)  
#Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)       #2000 5000
#str(Demo)
#saveRDS(Demo, paste(data_path, "SData_Integrated_rpca(Demo).RDS", sep = "/"))



# 循环聚类
response_col_name <- "Response.Comprehensive" # 指定 Response 列名
cluster_prefix <- "SCT_snn_res."

# 定义“不均衡”阈值 (例如，认为比例 > 0.8 或 < 0.2 是不均衡的)
imbalance_threshold_high <- 0.65
imbalance_threshold_low <- 0.35

# 存储结果的向量
unbalanced_resolutions <- c()
SData <- FindNeighbors(SData, reduction = "integrated.dr", verbose = FALSE)
resolution_sequence <- seq(0.1, 2, by = 0.1)
for (res in resolution_sequence) {
    cat("Processing resolution:", res, "\n") # 添加进度提示
    SData <- FindClusters(SData, resolution = res, verbose = FALSE)
    # 构造当前分辨率下的聚类列名
    cluster_col_name <- paste0(cluster_prefix, res)
    if (!cluster_col_name %in% colnames(SData@meta.data)) {
        warning(paste0("Expected cluster column '", cluster_col_name, "' not found in metadata. Skipping resolution ", res))
        next
    }
    # 提取聚类和响应信息
    meta_tmp <- SData@meta.data[, c(cluster_col_name, response_col_name)]
    colnames(meta_tmp) <- c("cluster", "response") # 简化列名以便处理
    # 检查 Response 列是否存在且是 Yes/No
    if (!response_col_name %in% colnames(SData@meta.data)) {
        stop(paste0("Response column '", response_col_name, "' not found in metadata."))
    }
}

# 计算每个 cluster 的 Response 比例
cluster_response_summary <- meta_tmp %>%
    group_by(cluster) %>%
    summarise(
        total_cells = n(),
        yes_count = sum(response == "Yes", na.rm = TRUE),
        no_count = sum(response == "No", na.rm = TRUE),
        prop_yes = yes_count / total_cells,
        .groups = "drop"
    )

# 检查是否有任何 cluster 的比例超出不平衡阈值
# 我们检查 Yes 的比例是否过高或过低
is_unbalanced <- any(
    !is.na(cluster_response_summary$prop_yes) &
        (cluster_response_summary$prop_yes > imbalance_threshold_high |
            cluster_response_summary$prop_yes < imbalance_threshold_low)
)

# 如果发现不均衡，记录分辨率
if (is_unbalanced) {
    unbalanced_resolutions <- c(unbalanced_resolutions, res)
    cat("  -> Unbalanced cluster found at resolution", res, "\n")
}


# --- 输出结果 ---
cat("\n--- Summary ---\n")
if (length(unbalanced_resolutions) > 0) {
    cat("Resolutions with unbalanced 'Response.Comprehensive' in at least one cluster:\n")
    print(unbalanced_resolutions)
} else {
    cat("No resolutions between 0.1 and 2 resulted in clusters with unbalanced 'Response.Comprehensive' (using thresholds <", imbalance_threshold_low, "or >", imbalance_threshold_high, ").\n")
}



saveRDS(SData, file = file.path(data_path, "SData_Integrated_ClusterCycle.RDS"))



#检查结果
SData_Result <- readRDS(file.path(data_path, "SData_Integrated_ClusterCycle.RDS"))      #6GB

cluster_prefix <- "SCT_snn_res." # 聚类列的前缀
response_col <- "Response.Comprehensive" # 响应列名
plot_output_dir <- file.path(plot_path, "Cluster_Cycle") # 图片输出目录


# 获取所有匹配的聚类列名
cluster_cols <- grep(cluster_prefix, colnames(SData_Result@meta.data), value = TRUE)
# 过滤掉可能不是数值分辨率的列（如果有必要）
# cluster_cols <- cluster_cols[grep("\\.\\d+", cluster_cols)] # 例如，只保留包含小数点的

# 检查响应列是否存在
if (!(response_col %in% colnames(SData_Result@meta.data))) {
  stop(glue::glue("Column '{response_col}' not found in SData_Result@meta.data"))
}

# 检查响应列的唯一值 (警告非 Yes/No 值)
unique_responses <- unique(SData_Result@meta.data[[response_col]])
cat("Unique values in", response_col, ":", unique_responses, "\n")
if (!all(unique_responses[!is.na(unique_responses)] %in% c("Yes", "No"))) {
  warning(glue::glue("'{response_col}' contains values other than 'Yes' or 'No'. Plotting might be affected."))
}

# --- 3. 循环处理每个分辨率并绘图 ---
all_plots <- list() # 可选：存储所有图对象用于后续操作

for (clust_col in cluster_cols) {
  # 从列名中提取分辨率数值 (假设格式为 "prefixX.X")
  # 这里使用正则表达式提取最后一个数字部分
  res_value_str <- stringr::str_extract(clust_col, "\\d+\\.?\\d*$")
  if (is.na(res_value_str) || res_value_str == "") {
    warning(glue::glue("Could not extract resolution value from column name '{clust_col}'. Skipping."))
    next
  }
  res_value <- as.numeric(res_value_str)
  if (is.na(res_value)) {
     warning(glue::glue("Extracted resolution value '{res_value_str}' from '{clust_col}' is not numeric. Skipping."))
     next
  }


  cat("Processing column:", clust_col, "(Resolution:", res_value, ")\n")

  # --- 4. 数据准备 ---
  # 提取当前聚类和响应数据
  plot_data <- SData_Result@meta.data %>%
    select(!!sym(clust_col), !!sym(response_col)) %>%
    rename(cluster_id = !!sym(clust_col), response = !!sym(response_col))

  # 检查数据
  if (nrow(plot_data) == 0) {
    warning(glue::glue("No data found for column '{clust_col}'. Skipping."))
    next
  }

  # 计算每个 cluster 的 Yes/No 计数和比例
  summary_data <- plot_data %>%
    filter(!is.na(cluster_id), !is.na(response)) %>% # 移除 NA 值
    count(cluster_id, response) %>% # 计算每个 cluster-response 组合的数量
    group_by(cluster_id) %>%
    mutate(total = sum(n), prop = n / total) %>% # 计算总数和比例
    ungroup() %>%
    select(cluster_id, response, prop) # 只保留 cluster_id, response, prop

  # --- 5. 数据重塑 (Pivot) ---
  # 将数据从长格式转换为宽格式，使每个 response 类别成为一列
  # 然后再转回长格式，以便 ggplot 正确堆叠
  plot_data_wide <- summary_data %>%
    tidyr::pivot_wider(names_from = response, values_from = prop, values_fill = 0) # 用0填充缺失的组合

  # 重新整理为 ggplot 堆积柱状图所需的长格式
  # 确保列名是 "Yes" 和 "No" (根据你的实际数据调整)
  # 如果某些响应类别不存在，pivot_wider 会创建值为 0 的列
  expected_responses <- c("Yes", "No")
  for(resp in expected_responses) {
    if(!(resp %in% colnames(plot_data_wide))) {
      plot_data_wide[[resp]] <- 0 # 添加缺失的响应列并填充0
    }
  }

  # 最终用于绘图的长格式数据
  plot_data_final <- plot_data_wide %>%
    tidyr::pivot_longer(cols = all_of(expected_responses), names_to = "response", values_to = "proportion")


  # --- 6. 绘图 ---
  p <- ggplot(plot_data_final, aes(x = factor(cluster_id), y = proportion, fill = response)) +
    geom_col(position = "stack") + # 堆积柱状图
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + # Y轴为百分比，范围 0-100%
    scale_fill_manual(values = c("Yes" = "#66C2A5", "No" = "#FC8D62")) + # 自定义颜色 (可选)
    labs(
      title = glue::glue("Response Distribution per Cluster (Resolution = {res_value})"),
      x = "Cluster ID",
      y = "Proportion (%)",
      fill = "Response"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # X轴标签倾斜
      panel.grid.minor.y = element_blank() # 移除 Y 轴次网格线使图更清晰
    )

  # --- 7. 保存图形 ---
  filename <- file.path(plot_output_dir, glue::glue("Cluster_Response_Res_", stringr::str_replace_all(res_value_str, '\\.', '_'), ".pdf"))
  # 使用 str_replace_all 将文件名中的点替换为下划线，避免潜在问题
  ggsave(filename, plot = p, width = 12, height = 8, dpi = 300) # 调整尺寸和 DPI
  cat("Saved plot to:", filename, "\n")

  # --- 8. (可选) 存储图形对象 ---
  # all_plots[[length(all_plots) + 1]] <- p

  # --- 9. (可选) 在 RStudio 中显示图形 ---
  # print(p) # 注意：大量图形可能导致 RStudio 响应缓慢

}

cat("Finished generating plots.\n")
# 如果存储了所有图形对象，可以用 patchwork 或 gridExtra 等包组合它们
# 例如 (组合前几个):
# combined_plot <- wrap_plots(all_plots[1:min(4, length(all_plots))], ncol = 2)
# print(combined_plot)




#图片结果
#resolution = 0.1时，有5类
#resolution = 0.5时，有16类
#resolution = 1时，有27类
#resolution = 2时， 有44类

#resolution = 0.1时，类2， No占比70%    总共有（6）类
#0.2，2，70                     （9）
#0.3，2，70                 （总共11类）
#0.4，2/12，70/70       （总共13类）
#0.5，11，75       （总共16类）
#0.6，11，75      （总共19类）      此处聚类个数有一个显著增加，故resolution暂定为0.6
#0.7，1/14，70/75     （总共24类）
#0.8，13，75            （24）
#0.9，13，75            （25）
#1，5/24/26，75/70/70       （27）
#1.1，                          （27）
#1.2，                        （31）





#检查Pre Or Post
cluster_prefix <- "SCT_snn_res." # 聚类列的前缀
response_col <- "PreOrPost" # 响应列名
plot_output_dir <- file.path(plot_path, "Cluster_Cycle") # 图片输出目录


# 获取所有匹配的聚类列名
cluster_cols <- grep(cluster_prefix, colnames(SData_Result@meta.data), value = TRUE)
# 过滤掉可能不是数值分辨率的列（如果有必要）
# cluster_cols <- cluster_cols[grep("\\.\\d+", cluster_cols)] # 例如，只保留包含小数点的

# 检查响应列是否存在
if (!(response_col %in% colnames(SData_Result@meta.data))) {
  stop(glue::glue("Column '{response_col}' not found in SData_Result@meta.data"))
}

# 检查响应列的唯一值 (警告非 Yes/No 值)
unique_responses <- unique(SData_Result@meta.data[[response_col]])
cat("Unique values in", response_col, ":", unique_responses, "\n")
if (!all(unique_responses[!is.na(unique_responses)] %in% c("Pre", "Post", "Prog"))) {
  warning(glue::glue("'{response_col}' contains values other than 'Yes' or 'No'. Plotting might be affected."))
}

# --- 3. 循环处理每个分辨率并绘图 ---
all_plots <- list() # 可选：存储所有图对象用于后续操作

for (clust_col in cluster_cols) {
  # 从列名中提取分辨率数值 (假设格式为 "prefixX.X")
  # 这里使用正则表达式提取最后一个数字部分
  res_value_str <- stringr::str_extract(clust_col, "\\d+\\.?\\d*$")
  if (is.na(res_value_str) || res_value_str == "") {
    warning(glue::glue("Could not extract resolution value from column name '{clust_col}'. Skipping."))
    next
  }
  res_value <- as.numeric(res_value_str)
  if (is.na(res_value)) {
     warning(glue::glue("Extracted resolution value '{res_value_str}' from '{clust_col}' is not numeric. Skipping."))
     next
  }


  cat("Processing column:", clust_col, "(Resolution:", res_value, ")\n")

  # --- 4. 数据准备 ---
  # 提取当前聚类和响应数据
  plot_data <- SData_Result@meta.data %>%
    select(!!sym(clust_col), !!sym(response_col)) %>%
    rename(cluster_id = !!sym(clust_col), response = !!sym(response_col))

  # 检查数据
  if (nrow(plot_data) == 0) {
    warning(glue::glue("No data found for column '{clust_col}'. Skipping."))
    next
  }

  # 计算每个 cluster 的 Yes/No 计数和比例
  summary_data <- plot_data %>%
    filter(!is.na(cluster_id), !is.na(response)) %>% # 移除 NA 值
    count(cluster_id, response) %>% # 计算每个 cluster-response 组合的数量
    group_by(cluster_id) %>%
    mutate(total = sum(n), prop = n / total) %>% # 计算总数和比例
    ungroup() %>%
    select(cluster_id, response, prop) # 只保留 cluster_id, response, prop

  # --- 5. 数据重塑 (Pivot) ---
  # 将数据从长格式转换为宽格式，使每个 response 类别成为一列
  # 然后再转回长格式，以便 ggplot 正确堆叠
  plot_data_wide <- summary_data %>%
    tidyr::pivot_wider(names_from = response, values_from = prop, values_fill = 0) # 用0填充缺失的组合

  # 重新整理为 ggplot 堆积柱状图所需的长格式
  # 确保列名是 "Yes" 和 "No" (根据你的实际数据调整)
  # 如果某些响应类别不存在，pivot_wider 会创建值为 0 的列
  expected_responses <- c("Pre", "Post", "Prog")
  for(resp in expected_responses) {
    if(!(resp %in% colnames(plot_data_wide))) {
      plot_data_wide[[resp]] <- 0 # 添加缺失的响应列并填充0
    }
  }

  # 最终用于绘图的长格式数据
  plot_data_final <- plot_data_wide %>%
    tidyr::pivot_longer(cols = all_of(expected_responses), names_to = "PreOrPost", values_to = "proportion")


  # --- 6. 绘图 ---
  p <- ggplot(plot_data_final, aes(x = factor(cluster_id), y = proportion, fill = PreOrPost)) +
    geom_col(position = "stack") + # 堆积柱状图
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + # Y轴为百分比，范围 0-100%
    scale_fill_manual(values = c("Pre" = "#66C2A5", "Post" = "#FC8D62", "Prog" = "#50bcdf")) + # 自定义颜色 (可选)
    labs(
      title = glue::glue("PreOrPost Distribution per Cluster (Resolution = {res_value})"),
      x = "Cluster ID",
      y = "Proportion (%)",
      fill = "PreOrPost"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # X轴标签倾斜
      panel.grid.minor.y = element_blank() # 移除 Y 轴次网格线使图更清晰
    )

  # --- 7. 保存图形 ---
  filename <- file.path(plot_output_dir, glue::glue("Cluster_PreOrPost_Res_", stringr::str_replace_all(res_value_str, '\\.', '_'), ".pdf"))
  # 使用 str_replace_all 将文件名中的点替换为下划线，避免潜在问题
  ggsave(filename, plot = p, width = 12, height = 8, dpi = 300) # 调整尺寸和 DPI
  cat("Saved plot to:", filename, "\n")

  # --- 8. (可选) 存储图形对象 ---
  # all_plots[[length(all_plots) + 1]] <- p

  # --- 9. (可选) 在 RStudio 中显示图形 ---
  # print(p) # 注意：大量图形可能导致 RStudio 响应缓慢

}

cat("Finished generating plots.\n")


#Reseponse图片结果
#resolution = 0.1时，有5类
#resolution = 0.5时，有16类
#resolution = 1时，有27类
#resolution = 2时， 有44类

#resolution = 0.1时，类2， No占比70%    总共有（6）类
#0.2，2，70                     （9）
#0.3，2，70                 （总共11类）
#0.4，2/12，70/70       （总共13类）
#0.5，11，75       （总共16类）
#0.6，11，75      （总共19类）      此处聚类个数有一个显著增加，故resolution暂定为0.6或者0.5
#0.7，1/14，70/75     （总共24类）
#0.8，13，75            （24）
#0.9，13，75            （25）
#1，5/24/26，75/70/70       （27）
#1.1，                          （27）
#1.2，                        （31）



#PreOrPost图片结果
#resolution = 0.1时，类2， Post占比70%
#0.2，3，70
#0.3，2，70                 与Reseponse同异
#0.4，1，70
#0.5，5/11/14，70/78/30     与Reseponse同异
#0.6，11/16，75/30          与Reseponse同异
#0.7，2/14/21，70/75/30     与Reseponse同异
#0.8，3/13/21，70/85/30     与Reseponse同异
#0.9，3/13/21，70/75/30     与Reseponse同异
#1，0/5/22，70/70/30        与Reseponse同异
