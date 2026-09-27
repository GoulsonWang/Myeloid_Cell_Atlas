#2025/9/16      计算percent.MT并添加至meta.data中

# Input: SeuratData.RDS(TNBC2数据)   地址为/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC2/SeuratData.RDS
# Output: 
#   1.SeuratData_AddPercentMT.RDS   运行失败！



suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")

SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC2/SeuratData.RDS")    #2GB

SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")      ##为什么全是NA值？？？
str(SData)

saveRDS(SData, file = paste(data_path, "SeuratData_AddPercentMT.RDS", sep = "/"))
