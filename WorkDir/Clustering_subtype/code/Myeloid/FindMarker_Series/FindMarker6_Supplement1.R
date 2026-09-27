#2025/11/18    本脚本用于对FindMarker6.R做补充，主要画聚类树状图



#input：
#   1.SData_Integrated_rpca.RDS

#output：
#   1.SData_Cluster_Cycle.RDS   res以0.5为间隔
#   2.ClusterTree.pdf



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
SData <- readRDS(file.path(data_path, "SData_Integrated_rpca.RDS"))



#聚类
#SData_Cluster <- FindNeighbors(SData, reduction = "integrated.dr", dims = 1:40, verbose = FALSE)
#resolution_sequence <- seq(0.05, 1, by = 0.05)
#for (res in resolution_sequence) {
#    cat("Processing resolution:", res, "\n") # 添加进度提示
#    SData_Cluster <- FindClusters(SData_Cluster, resolution = res, verbose = FALSE)
#}
#saveRDS(SData_Cluster, file = file.path(data_path, "SData_Cluster_Cycle.RDS"))

SData_Cluster <- readRDS(file.path(data_path, "SData_Cluster_Cycle.RDS"))
library(clustree)           #不知道什么原因，必须library，不能用clustree：：clustree的方式。
p1 <- clustree(SData_Cluster, prefix = "SCT_snn_res.") + coord_flip()
p = p1 + patchwork::plot_layout(widths = c(3,1))
ggsave(filename = file.path(plot_path, "ClusterTree.pdf"), p, width = 30, height = 14)



#由图决定，将res暂定为0.2，（13个子类）     
