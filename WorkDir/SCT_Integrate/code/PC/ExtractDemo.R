#2025/9/15      从Integrate的数据中提取demo

# Input: SData_Integrated.RDS   #整合后的数据
# Output: 
#   1.SData_Integrated(Demo).RDS
#   2.
    

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")

#载入数据
SData <- readRDS(paste(data_path, "SData_Integrated.RDS", sep = "/"))       #5GB
dim(SData)      #2000 62750
str(SData)

DemoCellName <- sample(colnames(SData), size = 3000)
Demo <- subset(SData, cells = DemoCellName)
dim(Demo)       #2000 3000
str(Demo)

saveRDS(Demo, paste(data_path, "SData_Integrated(Demo).RDS"), sep = "/")
