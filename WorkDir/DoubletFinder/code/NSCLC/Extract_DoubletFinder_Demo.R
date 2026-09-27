suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(future))
suppressMessages(library(DoubletFinder))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/code/")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder"
code_path <- paste(file_path, "/code/NSCLC", sep = "")
data_path <- paste(file_path, "/data/NSCLC", sep = "")
plot_path <- paste(file_path, "/plot/NSCLC", sep = "")

#载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/NSCLC/SeuratData(Rename).RDS")   
SDataList <- SplitObject(SData, split.by = "orig.ident")
Demo <- list()
Demo[[1]] <- SDataList[[1]]
Demo[[2]] <- SDataList[[2]]
Demo[[3]] <- SDataList[[3]]

saveRDS(Demo, file = file.path(data_path, "DoubletFinder_Demo.RDS"))