#2025/11/9
#本脚本用与检查转移灶和原发灶的细胞比例



# Input: 
#   1.Description1.RDS

# Output: 



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "Description1.RDS"))   #21GB
table(SData$PrimaryOrMet)
#结果：转移灶细胞占总细胞比例大于1/4