#2025/9/14      计算percent.MT并添加至meta.data中

# Input: SeuratDataAML(Rename).RDS   (AML数据)
# Output: 
#   1.SeuratDataAML(Rename)_AddPercentMT.RDS



suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/AML", sep = "")
data_path <- paste(file_path, "/data/AML", sep = "")
plot_path <- paste(file_path, "/plot/AML", sep = "")

SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/AML/SeuratDataAML(Rename).RDS")    #2GB

SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")
str(SData)

saveRDS(SData, file = paste(data_path, "SeuratDataAML(Rename)_AddPercentMT.RDS", sep = "/"))
