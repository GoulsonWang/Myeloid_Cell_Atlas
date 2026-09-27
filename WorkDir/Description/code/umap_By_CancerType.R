#2025/12/9
#本脚本用于画不同癌症类型的细胞所构成的umap图


#input
#   1.SData_Cluster030.RDS
#   2.

#output
#   1. umap_By_CancerType_ESCC.pdf    （每个CancerType的单独UMAP图）
#   2. umap_By_CancerType_NSCLC.pdf
#   3. umap_By_CancerType_PC.pdf
#   4. umap_By_CancerType_TNBC1.pdf



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
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/SData_Cluster030.RDS")    #6GB

# 获取所有唯一的CancerType
cancer_types <- unique(SData@meta.data$CancerType)
cat("找到的CancerType有:", cancer_types, "\n")

# 为每个CancerType创建单独的UMAP图，保持原始分组着色
for (cancer_type in cancer_types) {
  # 子集化数据，只包含当前指定的CancerType
  cells_to_plot <- colnames(SData)[SData@meta.data$CancerType == cancer_type]
  if (length(cells_to_plot) == 0) {
    cat("警告：未找到属于", cancer_type, "的细胞，跳过绘图。\n")
    next
  }
  
  SData_subset <- subset(SData, cells = cells_to_plot)
  
  # 绘制UMAP图，保留原始的聚类分组着色（例如 seurat_clusters）
  p <- DimPlot(SData_subset, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
       ggtitle(paste("UMAP plot for CancerType:", cancer_type))
  
  # 保存图像
  filename <- paste(plot_path, "/umap_By_CancerType_", cancer_type, ".pdf", sep = "")
  ggsave(filename, plot = p, width = 10, height = 10)
  cat("已保存", cancer_type, "的UMAP图\n")
  
  # 清理临时数据
  rm(SData_subset)
}
