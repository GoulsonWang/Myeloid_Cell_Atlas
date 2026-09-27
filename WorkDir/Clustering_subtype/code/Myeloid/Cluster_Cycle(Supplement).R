#2025/11/16   本脚本用与循环聚类的补充，主要为了输出聚类树状图。



#input：
#   1.SData_Integrated_ClusterCycle.RDS     

#output：
#   1.ClusterTree.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(clustree))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "SData_Integrated_ClusterCycle.RDS"))
p1 = clustree(SData,prefix = "SCT_snn_res.")+coord_flip()
p = p1 + patchwork::plot_layout(widths = c(3,1))
ggsave(p, filename = file.path(data_path,"ClusterTree.pdf"), p, width = 30, height = 14)
