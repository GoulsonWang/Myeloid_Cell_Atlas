#本脚本用与提取测试脚本时所用的数据Demo，利用Demo跑通脚本后再输入数据提交作业。

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")


#dim(SData)     
#str(SData)



#对细胞随机取样
#SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/ESCC/SeuratDataList(QC).RDS")       #5GB
#DemoCellName <- sample(colnames(SData), size = 3000)    #此处的SData已经进行过了FindVariableFeatures()
#Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)       #2000 3000
#str(Demo)
#saveRDS(Demo, paste(data_path, "SData_Integrated(Demo).RDS", sep = "/"))



#取前俩样本
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/ESCC/SeuratDataList(QC).RDS")
SDataListDemo <- list()
for(i in 1:2){
    SDataListDemo[[i]] <- SDataList[[i]]
}
str(SDataListDemo)  #结果保存为str(SDataList_Demo).txt
saveRDS(SDataListDemo, paste(data_path, "SDataList_Demo.RDS", sep = "/"))
