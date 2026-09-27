#2025/12/23
#本脚本用与构建用与GSEA的input



#input
#   1.FindMarker_Prep_Res030.RDS

#output
#   1.FCList_MacroMono.rds
#   2.FCList_DC.rds
#   3.FCList_Mono.rds
#   4.FCList_MastNeut.rds
#   5.Minor_Cell_Type_UMAP.pdf
#   6.SData_Minor_Cell_Type.RDS         Minor_Cell_Type结果



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))



# 设置工作目录
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/code/")

# 定义路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/FindMarker_Prep_Res030.RDS")
SData@active.assay      #"SCT"
#SData中SCT_snn_res.0.3为4的为Mast细胞，5，15，17为DC细胞，其余为Macro/Mono细胞（共计14类）

table(SData@meta.data$SCT_snn_res.0.3)
#    0     1     2     3     4     5     6     7     8     9    10    11    12 
#14174  9925  9061  8362  6019  5770  5548  5264  3108  2526  2002  1360  1072 
#   13    14    15    16    17 
#  910   831   735   613   526 

Marker <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/Macro_Marker_Result_Res030_Params.RDS")
Marker$min_pct_0.45_min_diff_pct_0.1["8"][[1]] %>% rownames()
#尝试对Macro/Mono的亚类进行合并
#   结合FindMerker6_Supplement4.R中的结果：
#       0大概率为Macro，激活状态偏向M2型或混合型
#       1和13归为一类，大概率为Macro，激活状态偏向M2型
#       2大概率为活化的Macro，激活状态为M1型
#       3大概率为浸润后的中性粒细胞，且促炎，肺癌特有的
#       4为Mast
#       5为DC
#       6Mono混合M1型
#       7大概率为Mono
#       8缺乏特异性基因表达
#       9疑似肺泡驻留巨噬细胞，肺癌特有的
#       10强烈的干扰素应答，偏M1（IFN-γ激活）
#       11活化的单核细胞，或单核细胞——巨噬细胞的中间态
#       12激活的单核细胞
#       13巨噬细胞，偏向M2型
#       14促肿瘤 M2巨噬细胞
#       15DC
#       16是热休克蛋白高表达簇，无意义
#       17DC
#   合并规则：将7，11合为一类(Mono);1和13,14归为一类(M2Macro)，0和10归为一类，16删除，8删除。
#   Macro/Mono共计9个子类，其中2个是肺癌独有的，其余7类中2类Mono，5类Macro。

# 删除cluster 8和cluster 16的细胞
# 首先找出要删除的细胞的名称
cells_to_remove <- colnames(SData)[SData@meta.data$SCT_snn_res.0.3 %in% c("8", "16")]
cat("将删除", length(cells_to_remove), "个细胞 (cluster 8和16)\n")

# 然后从SData中剔除这些细胞
SData_clean <- SData[, !colnames(SData) %in% cells_to_remove]

# 使用ifelse函数的嵌套结构为每个细胞添加Minor_Cell_Type注释（使用字符串命名格式）
SData_clean$Minor_Cell_Type <- ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 4, "Mast",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 5, "DC_1",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 15, "DC_2",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 17, "DC_3",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 %in%c(0, 10), "Macro_1",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 %in% c(1, 13, 14), "Macro_2",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 2, "Macro_3",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 3, "Neutro_1",  
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 6, "Macro_4",                                                
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 %in% c(7, 11), "Mono_1",
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 9, "Macro_5",                                                    
                         ifelse(SData_clean@meta.data$SCT_snn_res.0.3 == 12, "Mono_2",  
                         as.character(SData@meta.data$SCT_snn_res.0.3)))))))))))))
# 查看新添加的元数据分布
cat("细胞类型分布:\n")
print(table(SData_clean$Minor_Cell_Type))

# 生成并保存UMAP图
umapdata <- Embeddings(SData_clean, reduction = "umap")
MinorInfo <- SData_clean@meta.data$Minor_Cell_Type
ggplotDataframe <- data.frame(umapdata, MinorInfo)

# 定义高对比度颜色调色板（使用基础R函数）
n_clusters <- length(unique(MinorInfo))
# 手动定义一组高对比度颜色
base_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00", 
                 "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", 
                 "#FFFF99", "#B15928", "#F0F0F0", "#A0522D", "#FF69B4",
                 "#8B4513", "#40E0D0", "#6B8E23", "#FF1493", "#00008B")

if(n_clusters <= length(base_colors)) {
  cluster_colors <- base_colors[1:n_clusters]
} else {
  # 如果聚类数量超过预定义颜色数量，使用rainbow生成
  cluster_colors <- rainbow(n_clusters)
}

# 为每个细胞类型分配颜色
names(cluster_colors) <- unique(MinorInfo)

# 总图 - 所有癌症类型细胞的UMAP图
p_total <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(fill = MinorInfo), size = 3.5, shape = 21, stroke = 0.3, colour = "black", alpha = 0.7) +  # 使用fill进行填充，colour作为边缘
        scale_fill_manual(values = cluster_colors) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)) + 
        labs(title = "UMAP visualization of Minor Cell Types")  # 添加居中标题

# 计算每个cluster的中心位置（最密集区域）
total_cluster_centers <- ggplotDataframe %>%
  group_by(MinorInfo) %>%
  summarise(
    center_x = median(umap_1, na.rm = TRUE),
    center_y = median(umap_2, na.rm = TRUE),
    .groups = 'drop'
  )

# 添加标签，使用白色背景
p_total <- p_total +
  geom_label(
    data = total_cluster_centers,
    aes(x = center_x, y = center_y, label = MinorInfo),
    inherit.aes = FALSE,
    colour = "white",            # 白色边框
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    label.padding = unit(0.1, "lines"),  # 标签内边距
    label.r = unit(0.15, "lines"),       # 圆角
    label.size = 0.25,                   # 标签边框大小
    show.legend = FALSE                  # 不显示标签图例
  ) +
  geom_text(
    data = total_cluster_centers,
    aes(x = center_x, y = center_y, label = MinorInfo),
    inherit.aes = FALSE,
    colour = "black",            # 黑色文字
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    show.legend = FALSE          # 不显示标签图例
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3, stroke = 0.5, colour = "black")))  # 调整图例点大小和边框

# 保存总图
total_filename <- paste(plot_path, "/Minor_Cell_Type_UMAP.pdf", sep = "")
ggsave(total_filename, plot = p_total, width = 12, height = 10)

print("Minor_Cell_Type元数据已成功添加到SData对象中")
saveRDS(SData_clean, file = paste0(data_path, "/SData_Minor_Cell_Type.RDS"))



#在Minor_Cell_Type基础上进行计算FoldChange
SData_clean$Minor_Cell_Type <- as.factor(SData_clean$Minor_Cell_Type)
SData_clean@active.ident <- SData_clean$Minor_Cell_Type
table(SData_clean$Minor_Cell_Type)



# 分别对Macro Mono和DC各群体计算FoldChange
SData_Macro <- subset(SData_clean, idents = c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5"))
SData_Macro@active.ident <- SData_Macro$Minor_Cell_Type
SData_Mono <- subset(SData_clean, idents = c("Mono_1", "Mono_2"))
SData_Mono@active.ident <- SData_Mono$Minor_Cell_Type
SData_DC <- subset(SData_clean, idents = c("DC_1", "DC_2", "DC_3"))
SData_DC@active.ident <- SData_DC$Minor_Cell_Type
# 计算FoldChange，并将结果存储为列表
#1.Macro
macro_levels <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5")
FCList_Macro <- list()
cat("开始计算Macro/Mono细胞亚群的FoldChange...\n")
for(level in macro_levels){
  cat("正在计算", level, "的FoldChange...\n")
  FCList_Macro[[level]] <- FoldChange(SData_Macro, ident.1 = level)
}
cat("Macro细胞亚群的FoldChange计算完成\n")

#2.DC
dc_levels <- c("DC_1", "DC_2", "DC_3")
FCList_DC <- list()
cat("开始计算DC细胞亚群的FoldChange...\n")
for(level in dc_levels){
  cat("正在计算", level, "的FoldChange...\n")
  FCList_DC[[level]] <- FoldChange(SData_DC, ident.1 = level)
}
cat("DC细胞亚群的FoldChange计算完成\n")

#3.Mono
mono_levels <- c("Mono_1", "Mono_2")
FCList_Mono <- list()
cat("开始计算Mono细胞亚群的FoldChange...\n")
for(level in mono_levels){
    cat("正在计算", level, "的FoldChange...\n")
    FCList_Mono[[level]] <- FoldChange(SData_Mono, ident.1 = level)
}
cat("Mono细胞亚群的FoldChange计算完成\n")

#4.Mast和Neutrophils        此时以所有其他髓系细胞作为对照
MastNuet_levels <- c("Mast", "Neutro_1")
FCList_MastNeut <- list()
for(level in MastNuet_levels){
    cat("开始计算", level, "的FoldChange...\n")
    FCList_MastNeut[[level]] <- FoldChange(SData_clean, ident.1 = level)
}
cat("Mast和Neutrophils的FoldChange计算完成\n")

# 输出结果信息
print("Macro/Mono细胞的FoldChange列表长度:")
print(length(FCList_Macro))
print("DC细胞的FoldChange列表长度:")
print(length(FCList_DC))
print("Mono细胞的FoldChange列表长度:")
print(length(FCList_Mono))
print("Mast和Neutrophils细胞的FoldChange列表长度:")
print(length(FCList_MastNeut))

# 保存结果
saveRDS(FCList_Macro, file = paste0(data_path, "/FCList_Macro.rds"))
saveRDS(FCList_DC, file = paste0(data_path, "/FCList_DC.rds"))
saveRDS(FCList_Mono, file = paste0(data_path, "/FCList_Mono.rds"))
saveRDS(FCList_MastNeut, file = paste0(data_path, "/FCList_MastNeut.rds"))

print("FoldChange结果已保存到数据文件夹中")
