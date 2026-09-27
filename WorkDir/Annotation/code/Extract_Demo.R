#本脚本用与提取测试脚本时所用的数据Demo，利用Demo跑通脚本后再输入数据提交作业。

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#对细胞随机取样
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/TotalData/SData_Integrated_Filtered_Processed.RDS")       #5GB
DemoCellName <- sample(colnames(SData), size = 5000)    #此处的SData已经进行过了FindVariableFeatures()
Demo <- subset(SData, cells = DemoCellName)
dim(Demo)       #2000 3000
str(Demo)
saveRDS(Demo, paste(data_path, "SData_Integrated_Filtered_Processed(Demo).RDS", sep = "/"))