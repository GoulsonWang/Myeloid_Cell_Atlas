#2025/12/17
#本脚本用于美化不同癌症类型的细胞所构成的umap图，以及所有癌症类型细胞的总umap图


#input
#   1.SData030.RDS
#   2.

#output
#   1. umap_By_CancerType_ALL.pdf
#   2. umap_By_CancerType_ESCC.pdf
#   3. umap_By_CancerType_NSCLC.pdf
#   4. umap_By_CancerType_PC.pdf
#   5. umap_By_CancerType_TNBC1.R



# 加载必要的库
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code/")

# 载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/SData030.RDS")    #6GB



#提取画图所需的数据
umapdata <- Embeddings(SData, reduction = "umap")
tsnedata <- Embeddings(SData, reduction = "tsne")
clusterdata <- SData@active.ident
cancertype <- SData@meta.data$CancerType
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata, cancertype)

# 获取所有唯一的CancerType
cancer_types <- unique(ggplotDataframe$cancertype)
cat("找到的CancerType有:", cancer_types, "\n")

# 为每个CancerType生成UMAP图
for (cancer_type in cancer_types) {
  # 筛选特定CancerType的数据
  subset_data <- ggplotDataframe[which(ggplotDataframe$cancertype == cancer_type), ]
  
  # 检查是否有数据
  if (nrow(subset_data) == 0) {
    cat("警告：未找到属于", cancer_type, "的细胞，跳过绘图。\n")
    next
  }
  
  # 创建UMAP图
  p <- ggplot(subset_data, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 21, aes(fill = clusterdata), size = 3.0, stroke = 0.1, colour = "black") +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(fill = guide_legend(override.aes = list(size = 3.5, stroke = 0.3))) + # 调整图例点大小
        labs(title = paste("UMAP visualization of", cancer_type, "cancer cells"))  # 添加居中标题

  # 计算每个cluster的中心位置（最密集区域）
  cluster_centers <- subset_data %>%
    group_by(clusterdata) %>%
    summarise(
      center_x = median(umap_1, na.rm = TRUE),
      center_y = median(umap_2, na.rm = TRUE),
      .groups = 'drop'
    )

  # 添加标签，使用双层标签结构确保可读性
  p <- p +
    geom_label(
      data = cluster_centers,
      aes(x = center_x, y = center_y, label = clusterdata, fill = clusterdata),
      inherit.aes = FALSE,
      colour = "white",            # 白色文字
      size = 4,                    # 标签字体大小
      fontface = "bold",           # 粗体
      vjust = 0.5,                 # 垂直居中
      hjust = 0.5,                 # 水平居中
      check_overlap = TRUE,        # 避免标签重叠
      label.padding = unit(0.3, "lines"),  # 标签内边距
      label.r = unit(0.15, "lines"),       # 圆角
      label.size = 0.8,                    # 标签边框大小
      show.legend = FALSE                  # 不显示标签图例
    ) +
    geom_text(
      data = cluster_centers,
      aes(x = center_x, y = center_y, label = clusterdata),
      inherit.aes = FALSE,
      colour = "white",            # 白色文字
      size = 4,                    # 标签字体大小
      fontface = "bold",           # 粗体
      vjust = 0.5,                 # 垂直居中
      hjust = 0.5,                 # 水平居中
      check_overlap = TRUE         # 避免标签重叠
    ) +
    guides(fill = "none")  # 完全移除图例
  
  # 显示图形
  print(p)
  
  # 保存图形
  filename <- paste(plot_path, "/umap_By_CancerType_", cancer_type, ".pdf", sep = "")
  ggsave(filename, plot = p, width = 10, height = 10)
  cat("已保存", cancer_type, "的UMAP图到", filename, "\n")
}

# 总图 - 所有癌症类型细胞的UMAP图
p_total <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 21, aes(fill = clusterdata), size = 3.0, stroke = 0.1, colour = "black") +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(fill = guide_legend(override.aes = list(size = 3.5, stroke = 0.3))) + # 调整图例点大小
        labs(title = "UMAP visualization of all cancer cells")  # 添加居中标题

# 计算每个cluster的中心位置（最密集区域）
total_cluster_centers <- ggplotDataframe %>%
  group_by(clusterdata) %>%
  summarise(
    center_x = median(umap_1, na.rm = TRUE),
    center_y = median(umap_2, na.rm = TRUE),
    .groups = 'drop'
  )

# 添加标签，使用双层标签结构确保可读性
p_total <- p_total +
  geom_label(
    data = total_cluster_centers,
    aes(x = center_x, y = center_y, label = clusterdata, fill = clusterdata),
    inherit.aes = FALSE,
    colour = "white",            # 白色文字
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    check_overlap = TRUE,        # 避免标签重叠
    label.padding = unit(0.3, "lines"),  # 标签内边距
    label.r = unit(0.15, "lines"),       # 圆角
    label.size = 0.8,                    # 标签边框大小
    show.legend = FALSE                  # 不显示标签图例
  ) +
  geom_text(
    data = total_cluster_centers,
    aes(x = center_x, y = center_y, label = clusterdata),
    inherit.aes = FALSE,
    colour = "white",            # 白色文字
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    check_overlap = TRUE         # 避免标签重叠
  ) +
  guides(fill = "none")  # 完全移除图例

# 显示总图
print(p_total)

# 保存总图
total_filename <- paste(plot_path, "/umap_By_CancerType_ALL.pdf", sep = "")
ggsave(total_filename, plot = p_total, width = 10, height = 10)
cat("已保存所有癌症类型的UMAP图到", total_filename, "\n")



#绘制Major_Cell_Type的umap图
major_annotations <- rep("Macro/Mono", ncol(SData))
current_clusters <- Idents(SData)
major_annotations[which(current_clusters == "4")] <- "Mast"
major_annotations[which(current_clusters %in% c("5", "15", "17"))] <- "DCs"
names(major_annotations) <- colnames(SData)
SData <- AddMetaData(SData, metadata = major_annotations, col.name = "Major_Cell_Type")
MajorCellTypeData <- SData@meta.data$Major_Cell_Type
ggplotDataFrame2 <- data.frame(umapdata, MajorCellTypeData)

p_total_MajorCellType <- ggplot(ggplotDataFrame2, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 21, aes(fill = MajorCellTypeData), size = 3.0, stroke = 0.1, colour = "black") +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(fill = guide_legend(override.aes = list(size = 3.5, stroke = 0.3))) + # 调整图例点大小
        labs(title = "UMAP visualization of Major Cell Types") +  # 修改标题以准确描述图表内容
        guides(fill = "none")  # 移除图例

# 计算每个Major Cell Type的中心位置（最密集区域）
major_cell_type_centers <- ggplotDataFrame2 %>%
  group_by(MajorCellTypeData) %>%
  summarise(
    center_x = median(umap_1, na.rm = TRUE),
    center_y = median(umap_2, na.rm = TRUE),
    .groups = 'drop'
  )

# 添加标签，使用双层标签结构确保可读性
p_total_MajorCellType <- p_total_MajorCellType +
  geom_label(
    data = major_cell_type_centers,
    aes(x = center_x, y = center_y, label = MajorCellTypeData, fill = MajorCellTypeData),
    inherit.aes = FALSE,
    colour = "white",            # 白色文字
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    check_overlap = TRUE,        # 避免标签重叠
    label.padding = unit(0.3, "lines"),  # 标签内边距
    label.r = unit(0.15, "lines"),       # 圆角
    label.size = 0.8,                    # 标签边框大小
    show.legend = FALSE                  # 不显示标签图例
  ) +
  geom_text(
    data = major_cell_type_centers,
    aes(x = center_x, y = center_y, label = MajorCellTypeData),
    inherit.aes = FALSE,
    colour = "white",            # 白色文字
    size = 4,                    # 标签字体大小
    fontface = "bold",           # 粗体
    vjust = 0.5,                 # 垂直居中
    hjust = 0.5,                 # 水平居中
    check_overlap = TRUE         # 避免标签重叠
  )

# 显示Major Cell Type图
print(p_total_MajorCellType)

# 保存Major Cell Type图
major_filename <- paste(plot_path, "/umap_By_Major_Cell_Type.pdf", sep = "")
ggsave(major_filename, plot = p_total_MajorCellType, width = 10, height = 10)
cat("已保存Major Cell Type的UMAP图到", major_filename, "\n")
