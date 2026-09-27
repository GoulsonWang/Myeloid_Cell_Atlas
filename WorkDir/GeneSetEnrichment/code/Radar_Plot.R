#2025/12/27     本脚本用与汇总Macro和Mono的富集结果



#input
#   1.KEGG_GSEA_Result.rds


#output
#   1.Macro_Mono_Top5_Pathways.csv

suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment")  

# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")

#载入数据——富集结果
gsea_result <- readRDS(file.path(data_path, "KEGG_GSEA_Result.rds"))

# 定义cluster名称
cluster_names <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", 
                   "Mono_1", "Mono_2", 
                   "DC_1", "DC_2", "DC_3", 
                   "Mast", 
                   "Neutro_1")



#合并富集结果
all_results_list <- list()
for(i in 1:length(gsea_result)) {
  # 检查当前 GSEA 结果是否有效
  if(!is.null(gsea_result[[i]]) && nrow(gsea_result[[i]]@result) > 0) { # 假设结果在 @result slot
    
    # 提取 S4 对象内部的 data.frame
    current_result_df <- gsea_result[[i]]@result 
    
    # 添加 cluster 列
    current_result_df$cluster <- cluster_names[i]  
    
    # 将处理好的 data.frame 添加到列表中
    all_results_list[[i]] <- current_result_df
    
  }
}
# 使用 do.call 和 rbind 来合并列表中的所有 data.frame
if(length(all_results_list) > 0) {
  combined_gsea_results <- do.call(rbind, all_results_list)
  
  # 重置行名（可选）
  rownames(combined_gsea_results) <- NULL
  
  print("GSEA 结果已成功整合到一个 data.frame 中。")
  print(paste("总共有", nrow(combined_gsea_results), "行数据。"))   #675个通路

} else {
  print("警告：没有有效的 GSEA 结果可以合并。")
  combined_gsea_results <- data.frame() # 创建一个空的 data.frame
}
#筛选p.adjust < 0.05
combined_gsea_results <- combined_gsea_results[combined_gsea_results$p.adjust < 0.05, ] #444条数据



#挑选出共有的通路（至少三个cluster共有）
combined_gsea_results_MacroMono <- subset(combined_gsea_results, cluster %in% c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2"))
shared_pathway_result_MacroMono <- combined_gsea_results_MacroMono %>%
   distinct(Description, cluster) %>%  # 确保每个通路-聚类组合只计算一次
   count(Description) %>%             # 计算每个通路出现的 cluster 数量
   filter(n >= 3) %>%                 # 筛选至少出现在 3 个 cluster 的通路
   pull(Description)                  # 提取通路名称


print(shared_pathway_result_MacroMono)

 # 获取这些通路的详细结果
if(length(shared_pathway_result_MacroMono) > 0) {
 filtered_results_dplyr <- combined_gsea_results %>%
     filter(Description %in% shared_pathway_result_MacroMono) %>%
     filter(cluster %in% c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2"))
  print(paste("共有", nrow(filtered_results_dplyr), "条记录属于这些共有通路。"))
  print("这些共有通路及其所在 cluster 的部分 GSEA 统计信息 (示例前10行):")
  print(filtered_results_dplyr[, c("Description", "cluster", "NES", "p.adjust")]) # 只显示部分列作为示例
} else {
  print("使用 dplyr 方法没有找到至少在 3 个 cluster 中富集的通路。")
}

#能用来描述巨噬细胞/单核细胞的通路：
#hsa04146: Autophagy - other **(自噬 - 其他)  ——区分稳态or应激
#hsa04625: C-type lectin receptor signaling pathway **(C型凝集素受体信号通路)
#hsa04110: Cell cycle （反映是否处于增殖状态/区分Mono和Macro）
#hsa00981: Cholesterol metabolism（区分M1 M2）
#hsa04060: Cytokine-cytokine receptor interaction（细胞因子受体相互作用，免疫交互、极化、激活的基础）
#hsa03030: DNA replication（反映是否处于增殖状态/区分Mono和Macro）
#hsa04657: IL-17 signaling pathway 应（髓系细胞可响应IL-17，增强炎症反应）
#hsa04064: NF-kappa B signaling pathway(核心)
#hsa04621: NOD-like receptor signaling pathway(核心)
#hsa04161: PD-L1 expression and PD-1 checkpoint pathway in cancer(PD-L1免疫治疗相关)
#hsa04668: TNF signaling pathway（核心）



# --- 定义核心通路名称和细胞类型 ---
core_pathways_desc <- c(
  "Autophagy - other",                  # hsa04146
  "C-type lectin receptor signaling pathway", # hsa04625
  "Cell cycle",                         # hsa04110
  "Cholesterol metabolism",             # hsa00981
  "Cytokine-cytokine receptor interaction", # hsa04060
  "DNA replication",                    # hsa03030
  "IL-17 signaling pathway",            # hsa04657 (假设 GSEA 结果中是这个格式，如果包含 "应" 字则需调整)
  "NF-kappa B signaling pathway",       # hsa04064
  "NOD-like receptor signaling pathway",# hsa04621
  "PD-L1 expression and PD-1 checkpoint pathway in cancer", # hsa04161
  "TNF signaling pathway"               # hsa04668
)

cell_types <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2")

# --- 筛选并重塑数据 ---
# 1. 筛选出感兴趣的通路和细胞类型
# 确保 combined_gsea_results$Description 中的名称与 core_pathways_desc 完全一致
filtered_data <- combined_gsea_results %>%
  filter(
    cluster %in% cell_types,
    Description %in% core_pathways_desc
  )

# 2. 检查筛选结果
if(nrow(filtered_data) == 0) {
  stop("没有找到符合条件的数据 (指定的通路和细胞类型)。")
  # 你可以打印 unique(filtered_data$Description) 来检查实际的通路名称
  # print(unique(combined_gsea_results$Description))
}

# 3. 计算每个通路在所有 cluster 中的 NES 最小值 (用于填充缺失值)
# 计算 *所有* cluster 中每个通路的最小 NES (基于 filtered_data 或原始 combined_gsea_results)
pathway_min_nes <- filtered_data %>%
  filter(Description %in% core_pathways_desc) %>% # 确保只考虑目标通路
  group_by(Description) %>%
  summarise(min_nes = min(NES, na.rm = TRUE), .groups = 'drop')

# 4. 创建一个包含所有 cluster-pathway 组合的完整框架
all_combinations <- expand.grid(
  cluster = cell_types,
  Description = core_pathways_desc,
  stringsAsFactors = FALSE
)

# 5. 将筛选后的数据与完整框架合并，并填充缺失值
plot_data_joined <- all_combinations %>%
  left_join(filtered_data, by = c("cluster", "Description")) %>%
  # 左连接后，不存在的组合 NES 会是 NA
  # 将 NA 值替换为该通路在所有 cluster 中的最小 NES
  left_join(pathway_min_nes, by = "Description") %>% # 将最小值表连接过来
  mutate(
    NES = ifelse(is.na(NES), min_nes, NES) # 如果 NES 是 NA，则用 min_nes 替换
  ) %>%
  select(-min_nes) # 移除辅助列 min_nes
plot_data_joined <- plot_data_joined[, 1:8]
colnames(plot_data_joined)
#[1] "cluster"         "Description"     "ID"              "setSize"        
#[5] "enrichmentScore" "NES"             "pvalue"          "p.adjust"       
#[9] "qvalue"          "rank"            "leading_edge"    "core_enrichment"
plot_data_joined




#2025/12/27     本脚本用与汇总Macro和Mono的富集结果



#input
#   1.KEGG_GSEA_Result.rds


#output
#   1.Macro_Mono_Top5_Pathways.csv

suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment")  

# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")

#载入数据——富集结果
gsea_result <- readRDS(file.path(data_path, "KEGG_GSEA_Result.rds"))

# 定义cluster名称
cluster_names <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", 
                   "Mono_1", "Mono_2", 
                   "DC_1", "DC_2", "DC_3", 
                   "Mast", 
                   "Neutro_1")



#合并富集结果
all_results_list <- list()
for(i in 1:length(gsea_result)) {
  # 检查当前 GSEA 结果是否有效
  if(!is.null(gsea_result[[i]]) && nrow(gsea_result[[i]]@result) > 0) { # 假设结果在 @result slot
    
    # 提取 S4 对象内部的 data.frame
    current_result_df <- gsea_result[[i]]@result 
    
    # 添加 cluster 列
    current_result_df$cluster <- cluster_names[i]  
    
    # 将处理好的 data.frame 添加到列表中
    all_results_list[[i]] <- current_result_df
    
  }
}
# 使用 do.call 和 rbind 来合并列表中的所有 data.frame
if(length(all_results_list) > 0) {
  combined_gsea_results <- do.call(rbind, all_results_list)
  
  # 重置行名（可选）
  rownames(combined_gsea_results) <- NULL
  
  print("GSEA 结果已成功整合到一个 data.frame 中。")
  print(paste("总共有", nrow(combined_gsea_results), "行数据。"))   #675个通路

} else {
  print("警告：没有有效的 GSEA 结果可以合并。")
  combined_gsea_results <- data.frame() # 创建一个空的 data.frame
}
#筛选p.adjust < 0.05
combined_gsea_results <- combined_gsea_results[combined_gsea_results$p.adjust < 0.05, ] #444条数据



#挑选出共有的通路（至少三个cluster共有）
combined_gsea_results_MacroMono <- subset(combined_gsea_results, cluster %in% c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2"))
shared_pathway_result_MacroMono <- combined_gsea_results_MacroMono %>%
   distinct(Description, cluster) %>%  # 确保每个通路-聚类组合只计算一次
   count(Description) %>%             # 计算每个通路出现的 cluster 数量
   filter(n >= 3) %>%                 # 筛选至少出现在 3 个 cluster 的通路
   pull(Description)                  # 提取通路名称


print(shared_pathway_result_MacroMono)

 # 获取这些通路的详细结果
if(length(shared_pathway_result_MacroMono) > 0) {
 filtered_results_dplyr <- combined_gsea_results %>%
     filter(Description %in% shared_pathway_result_MacroMono) %>%
     filter(cluster %in% c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2"))
  print(paste("共有", nrow(filtered_results_dplyr), "条记录属于这些共有通路。"))
  print("这些共有通路及其所在 cluster 的部分 GSEA 统计信息 (示例前10行):")
  print(filtered_results_dplyr[, c("Description", "cluster", "NES", "p.adjust")]) # 只显示部分列作为示例
} else {
  print("使用 dplyr 方法没有找到至少在 3 个 cluster 中富集的通路。")
}

#能用来描述巨噬细胞/单核细胞的通路：
#hsa04146: Autophagy - other **(自噬 - 其他)  ——区分稳态or应激
#hsa04625: C-type lectin receptor signaling pathway **(C型凝集素受体信号通路)
#hsa04110: Cell cycle （反映是否处于增殖状态/区分Mono和Macro）
#hsa00981: Cholesterol metabolism（区分M1 M2）
#hsa04060: Cytokine-cytokine receptor interaction（细胞因子受体相互作用，免疫交互、极化、激活的基础）
#hsa03030: DNA replication（反映是否处于增殖状态/区分Mono和Macro）
#hsa04657: IL-17 signaling pathway 应（髓系细胞可响应IL-17，增强炎症反应）
#hsa04064: NF-kappa B signaling pathway(核心)
#hsa04621: NOD-like receptor signaling pathway(核心)
#hsa04161: PD-L1 expression and PD-1 checkpoint pathway in cancer(PD-L1免疫治疗相关)
#hsa04668: TNF signaling pathway（核心）



# --- 定义核心通路名称和细胞类型 ---
core_pathways_desc <- c(
  "Autophagy - other",                  # hsa04146
  "C-type lectin receptor signaling pathway", # hsa04625
  "Cell cycle",                         # hsa04110
  "Cholesterol metabolism",             # hsa00981
  "Cytokine-cytokine receptor interaction", # hsa04060
  "DNA replication",                    # hsa03030
  "IL-17 signaling pathway",            # hsa04657 (假设 GSEA 结果中是这个格式，如果包含 "应" 字则需调整)
  "NF-kappa B signaling pathway",       # hsa04064
  "NOD-like receptor signaling pathway",# hsa04621
  "PD-L1 expression and PD-1 checkpoint pathway in cancer", # hsa04161
  "TNF signaling pathway"               # hsa04668
)

cell_types <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2")

# --- 筛选并重塑数据 ---
# 1. 筛选出感兴趣的通路和细胞类型
# 确保 combined_gsea_results$Description 中的名称与 core_pathways_desc 完全一致
filtered_data <- combined_gsea_results %>%
  filter(
    cluster %in% cell_types,
    Description %in% core_pathways_desc
  )

# 2. 检查筛选结果
if(nrow(filtered_data) == 0) {
  stop("没有找到符合条件的数据 (指定的通路和细胞类型)。")
  # 你可以打印 unique(filtered_data$Description) 来检查实际的通路名称
  # print(unique(combined_gsea_results$Description))
}

# 3. 计算每个通路在所有 cluster 中的 NES 最小值 (用于填充缺失值)
# 计算 *所有* cluster 中每个通路的最小 NES (基于 filtered_data 或原始 combined_gsea_results)
pathway_min_nes <- filtered_data %>%
  filter(Description %in% core_pathways_desc) %>% # 确保只考虑目标通路
  group_by(Description) %>%
  summarise(min_nes = min(NES, na.rm = TRUE), .groups = 'drop')

# 4. 创建一个包含所有 cluster-pathway 组合的完整框架
all_combinations <- expand.grid(
  cluster = cell_types,
  Description = core_pathways_desc,
  stringsAsFactors = FALSE
)

# 5. 将筛选后的数据与完整框架合并，并填充缺失值
plot_data_joined <- all_combinations %>%
  left_join(filtered_data, by = c("cluster", "Description")) %>%
  # 左连接后，不存在的组合 NES 会是 NA
  # 将 NA 值替换为该通路在所有 cluster 中的最小 NES
  left_join(pathway_min_nes, by = "Description") %>% # 将最小值表连接过来
  mutate(
    NES = ifelse(is.na(NES), min_nes, NES) # 如果 NES 是 NA，则用 min_nes 替换
  ) %>%
  select(-min_nes) # 移除辅助列 min_nes
plot_data_joined <- plot_data_joined[, 1:8]
colnames(plot_data_joined)
#[1] "cluster"         "Description"     "ID"              "setSize"        
#[5] "enrichmentScore" "NES"             "pvalue"          "p.adjust"       
#[9] "qvalue"          "rank"            "leading_edge"    "core_enrichment"
plot_data_joined




#画雷达图
# --- 雷达图绘制部分 ---

# 1. 准备绘图数据 (长格式)
# 确保 Description 的顺序一致 (可选，但推荐)
unique_descriptions <- unique(plot_data_joined$Description) # 获取所有核心通路的顺序

plot_data_for_plotting <- plot_data_joined %>%
  select(cluster, Description, NES) %>% # 只保留需要的列
  arrange(cluster, match(Description, unique_descriptions)) # 按 cluster 分组，并按 unique_descriptions 的顺序排列

# 2. 为了闭合多边形，选取每个 cluster 的第一行并添加到数据末尾
first_points <- plot_data_for_plotting %>%
  group_by(cluster) %>%
  slice(1) %>% # 选取每个 cluster 的第一行
  ungroup()

# 3. 合并原始数据和第一点，形成闭合图形的数据
plot_data_long <- bind_rows(plot_data_for_plotting, first_points) %>%
  arrange(cluster) # 确保顺序正确，每个 cluster 的数据是连续的，且以第一个点结尾

# 4. 检查 plot_data_long 的结构（可选，用于调试）
# print(head(plot_data_long))
# print(tail(plot_data_long))

# 5. 获取唯一的 Description 用于 x 轴和标签
# unique_descriptions 已经在步骤 1 中定义过了

# 6. 绘制雷达图
p <- ggplot(plot_data_long, aes(x = Description, y = NES, group = cluster, color = cluster)) +
  geom_polygon(fill = NA, size = 1, alpha = 0.7) + # 绘制多边形轮廓
  geom_point(size = 2) + # 绘制数据点
  geom_line(size = 0.8) + # 连接点形成线
  # 添加轴标签
  geom_text(
    data = data.frame(Description = unique_descriptions, NES = max(plot_data_long$NES) * 1.1), # 将标签放在最大值稍外的位置
    aes(x = Description, y = NES, label = Description),
    inherit.aes = FALSE, # 不继承上面的 aes，避免颜色干扰
    angle = 360 * (1:length(unique_descriptions)) / length(unique_descriptions) - 90, # 计算角度，使标签朝向正确
    hjust = 1.2, # 水平对齐
    vjust = 0.5  # 垂直对齐
  ) +
  scale_x_discrete(limits = unique_descriptions, expand = c(0, 0)) + # 设置 x 轴为离散变量，并设置范围
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # 扩展 y 轴范围，为标签留出空间
  coord_polar(theta = "x") + # 转换为极坐标
  labs(
    title = "Radar Chart of NES by Cluster and Pathway",
    color = "Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), # 移除轴标题
    axis.text.y = element_text(angle = 0, hjust = 1), # y 轴文本角度
    # axis.text.x = element_text(angle = 45, hjust = 1), # x 轴文本角度（在极坐标下效果可能不佳）
    panel.grid.major.x = element_line(color = "grey", linetype = "dashed", linewidth = 0.5), # 主网格线（径向）
    panel.grid.major.y = element_line(color = "grey", linetype = "solid", linewidth = 0.5),  # 主网格线（圆周）
    panel.grid.minor = element_blank(), # 移除次网格线
    legend.position = "right" # 图例位置
  )

# 7. 显示图形
print(p)

# 8. 如果你想保存图片
# ggsave(file.path(plot_path, "Macro_Mono_Pathway_RadarChart.png"), plot = p, width = 10, height = 8, dpi = 300)