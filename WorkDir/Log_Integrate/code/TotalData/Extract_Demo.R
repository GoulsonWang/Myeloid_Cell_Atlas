#本脚本用与提取测试脚本时所用的数据Demo，利用Demo跑通脚本后再输入数据提交作业。



# Input: 
#   1.PC-----SData_Integrated_Processed.RDS
#   2.TNBC1----SData_Integrated_Processed.RDS
#   3.ESCC----SData_Integrated_Processed.RDS
#   4.NSCLC----SData_Integrated_Processed.RDS

# Output: 四个数据集的Demo
#   1.ESCC(Demo).RDS
#   2.NSCLC(Demo).RDS
#   3.PC(Demo).RDS
#   4.TNBC1(Demo).RDS



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")



#对细胞随机取样
DataSet_Path <- c(
    "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/ESCC/SData_Integrated_Processed.RDS", 
    "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/NSCLC/SData_Integrated_Processed.RDS", 
    "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/PC/SData_Integrated_Processed.RDS", 
    "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/TNBC1/SData_Integrated_Processed.RDS"
)   #依次是ESCC， NSCLC， PC， TNBC1
DataSet_Name <- c("ESCC", "NSCLC", "PC", "TNBC1")
SDataList <- list()
for(i in 1:length(DataSet_Path)){
    SDataList[[i]] <- readRDS(DataSet_Path[i]) 
    DefaultAssay(SDataList[[i]]) <- "integrated"
    DemoCellName <- sample(colnames(SDataList[[i]]), size = 1000)
    Demo <- subset(SDataList[[i]], cells = DemoCellName)
    dim(Demo)
    print("输出str(Demo)")
    str(Demo)
    saveRDS(Demo, paste(data_path, "/", DataSet_Name[i], "(Demo).RDS", sep = ""))
}




#取前俩样本
#SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TotalData/SeuratDataList(QC).RDS")
#SDataListDemo <- list()
#for(i in 1:2){
#    SDataDemo[[i]] <- SDataList[[i]]
#}
#str(SDataListDemo)  #结果保存为str(SDataList_Demo).txt
#saveRDS(SDataListDemo, paste(data_path, "SDataList_Demo.RDS", sep = "/"))
