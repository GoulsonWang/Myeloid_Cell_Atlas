#2025/11/19     本脚本用与生成res0.25~0.3之间的聚类树状图



#input：
#   1.SData_Integrated_rpca.RDS

#output：
#   1.SData_Cluster_Cycle_Res025_030.RDS
#   2.ClusterTree_Res025_030.pdf   res以0.005为间隔，0.25~0.3之间



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid/FindMarker_Series")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid/FindMarker_Series", sep = "")
data_path <- paste(file_path, "/data/Myeloid/HSP_Exclude", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid/HSP_Exclude", sep = "")



#提交时使用
plan(multicore, workers = 30)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 10 * 1024^3)   



#载入数据
#SData <- readRDS(file.path(data_path, "SData_Integrated_rpca.RDS"))



#聚类
#SData_Cluster <- FindNeighbors(SData, reduction = "integrated.dr", dims = 1:40, verbose = F)
#resolution_sequence <- seq(0.25, 0.3, by = 0.005)
#for (res in resolution_sequence) {
#    cat("Processing resolution:", res, "\n") # 添加进度提示
3    SData_Cluster <- FindClusters(SData_Cluster, resolution = res, verbose = FALSE)
#}
#saveRDS(SData_Cluster, file = file.path(data_path, "SData_Cluster_Cycle_Res025_030.RDS"))

#library(clustree)           #不知道什么原因，必须library，不能用clustree：：clustree的方式。
#p1 <- clustree(SData_Cluster, prefix = "SCT_snn_res.") + coord_flip()
#p = p1 + patchwork::plot_layout(widths = c(3,1))
#ggsave(filename = file.path(plot_path, "ClusterTree_Res025_030.pdf"), p, width = 30, height = 14)




#聚类umap图
SData_Cluster <- readRDS(file.path(data_path, "SData_Cluster_Cycle_Res025_030.RDS"))
SData_Cluster <- RunUMAP(SData_Cluster, dims = 1:40, reduction = "integrated.dr", verbose = F)
umapdata <- Embeddings(SData_Cluster, reduction = "umap")
ClusterData <- SData_Cluster$SCT_snn_res.0.3
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
    vjust = -1,            # 垂直调整，使标签在点上方 (-1 在上方, 0.5 居中, 1 在下方)
    hjust = 0.5            # 水平调整，使标签居中 (0 左对齐, 0.5 居中, 1 右对齐)
  )
ggsave(
    filename = paste(plot_path, "/SData_Cluster03.pdf", sep = ""),
    plot = p1, 
    width = 9, 
    height = 7
)



#Res025时，Cluster0有一小簇离群，共16子类
#Res0255时，Cluster0依然未变，但总类别为15个
#Res03时看的最顺眼，一共18个子类