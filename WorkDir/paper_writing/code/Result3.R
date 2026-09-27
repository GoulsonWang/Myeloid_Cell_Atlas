#本脚本用于生成髓系umap图



suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#最终参与了Minor_Cell_Annotation的细胞
SData <- readRDS(file.path(data_path, "1PrepSCT_SData_Rename.RDS"))
table(SData$Minor_Cell_Type)

# 生成并保存UMAP图
umapdata <- Embeddings(SData, reduction = "umap")
MinorInfo <- SData@meta.data$Minor_Cell_Type
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
        geom_point(aes(fill = MinorInfo), size = 3.5, shape = 21, stroke = 0.2, colour = "black", alpha = 0.7) +  # 使用fill进行填充，colour作为边缘
        scale_fill_manual(values = cluster_colors) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)) #+ 
#        labs(title = "UMAP visualization of Minor Cell Types")  # 添加居中标题

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
