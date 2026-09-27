#2025/9/15      从Integrate的数据中提取demo

# Input: SDataList_Processed.RDS    #整合后的数据
# Output: 
#   1.SDataList_Processed(Demo).RDS
#   2.
    

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")

#载入数据
SData <- readRDS(paste(data_path, "SDataList_Processed.RDS", sep = "/"))       #6GB

#Demo选取规则：只保留列表的前五个单元
SDataDemo <- list()
for (i in 1:10){
    SDataDemo[[i]] <- SData[[i]]
}

saveRDS(SDataDemo, paste(data_path, "SDataList_Processed(Demo).RDS", sep = "/"))
