#本脚本用与提取测试脚本时所用的数据Demo，利用Demo跑通脚本后再输入数据提交作业。


suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")



SData <- readRDS(paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))       #5GB
#dim(SData)      #2000 167525
str(SData)      #assay中有27组counts

#core code
DemoCellName <- sample(colnames(SData), size = 3000)    #此处的SData已经进行过了FindVariableFeatures()
Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)       #2000 3000
str(Demo)

saveRDS(Demo, paste(data_path, "SData_Integrated_Processed(Demo).RDS", sep = "/"))
