#2025/11/6
#本脚本用与处理Integration1.R的报错：忘记在SCT后接RunPCA了



#input：
#   1.SData_Processed.RDS（未RunPCA（））

#output：
#   1.SData_Processed.RDS（RunPCA()后）



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



SData <- readRDS(file.path(data_path, "SData_Processed.RDS"))
SData <- RunPCA(SData, verbose = F)
saveRDS(SData, file = file.path(data_path, "SData_Processed.RDS"))