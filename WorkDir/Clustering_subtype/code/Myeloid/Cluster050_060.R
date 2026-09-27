#2025/11/11         暂定亚型聚类resolution取0.6或0.5



#input：
#   1.SData_Integrated_rpca.RDS     整合后未聚类的数据

#output：
#   1.SData_Integrated_Res050.RDS
#   2.SData_Integrated_Res060.RDS
#   3.Integration_Cluster_Res050.pdf
#   4.Integration_Cluster_Res060.pdf



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
SData <-  FindNeighbors(SData, reduction = "integrated.dr", verbose = F) 
SData_Res050 <- FindClusters(SData, resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:30, reduction = "integrated.dr", verbose = F) %>%
    RunTSNE(dims = 1:30, reduction = "integrated.dr", verbose = F)
SData_Res060 <- FindClusters(SData, resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:30, reduction = "integrated.dr", verbose = F) %>%
    RunTSNE(dims = 1:30, reduction = "integrated.dr", verbose = F)
saveRDS(SData_Res050, file = file.path(data_path, "SData_Integrated_Res050.RDS"))
saveRDS(SData_Res060, file = file.path(data_path, "SData_Integrated_Res060.RDS"))



#聚类umap图 Res050
umapdata <- Embeddings(SData_Res050, reduction = "umap")
ClusterData <- SData_Res050$SCT_snn_res.0.5
ggplotdataframe1 <- data.frame(umapdata, ClusterData)
p1 <- ggplot(ggplotdataframe1, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = ClusterData), size = 0.3, stroke = 0.2) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke  = 1.5))) 
cluster_centers <- ggplotdataframe1 %>%
  group_by(ClusterData) %>%
  summarise(
    center_x = mean(umap_1, na.rm = TRUE),
    center_y = mean(umap_2, na.rm = TRUE),
    .groups = 'drop' # 可选，控制分组行为
  )
p1 <- p1 +
  # 使用 geom_text 添加标签
  geom_text(
    data = cluster_centers, 
    aes(x = center_x, y = center_y, label = ClusterData),
    inherit.aes = FALSE, # 不继承主图的 aes，使用 data 和 aes 中指定的内容
    # 可选的美化参数：
    colour = "black",      # 标签颜色
    size = 3,              # 标签字体大小
    # check_overlap = TRUE, # 如果标签重叠，可以尝试隐藏部分标签
    vjust = -1,            # 垂直调整，使标签在点上方 (-1 在上方, 0.5 居中, 1 在下方)
    hjust = 0.5            # 水平调整，使标签居中 (0 左对齐, 0.5 居中, 1 右对齐)
  )
ggsave(
    filename = paste(plot_path, "/Integration_Cluster_Res050.pdf", sep = ""),
    plot = p1, 
    width = 9, 
    height = 7
)
#Res 060
umapdata <- Embeddings(SData_Res060, reduction = "umap")
ClusterData <- SData_Res060$SCT_snn_res.0.6
ggplotdataframe2 <- data.frame(umapdata, ClusterData)
p2 <- ggplot(ggplotdataframe2, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = ClusterData), size = 0.3, stroke = 0.2) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke  = 1.5))) 
cluster_centers <- ggplotdataframe2 %>%
  group_by(ClusterData) %>%
  summarise(
    center_x = mean(umap_1, na.rm = TRUE),
    center_y = mean(umap_2, na.rm = TRUE),
    .groups = 'drop' # 可选，控制分组行为
  )
p2 <- p2 +
  # 使用 geom_text 添加标签
  geom_text(
    data = cluster_centers, 
    aes(x = center_x, y = center_y, label = ClusterData),
    inherit.aes = FALSE, # 不继承主图的 aes，使用 data 和 aes 中指定的内容
    # 可选的美化参数：
    colour = "black",      # 标签颜色
    size = 3,              # 标签字体大小
    # check_overlap = TRUE, # 如果标签重叠，可以尝试隐藏部分标签
    vjust = -1,            # 垂直调整，使标签在点上方 (-1 在上方, 0.5 居中, 1 在下方)
    hjust = 0.5            # 水平调整，使标签居中 (0 左对齐, 0.5 居中, 1 右对齐)
  )
ggsave(
    filename = paste(plot_path, "/Integration_Cluster_Res060.pdf", sep = ""),
    plot = p2, 
    width = 9, 
    height = 7
)



#根据聚类结果来看，选择resolution为0.6的