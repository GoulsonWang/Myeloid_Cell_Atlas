#2025/9/15   利用Seurat-CCA-Anchor 方法正对数据集内部进行整合

# Input: SDataList_Processed(QC3).RDS   #2GB大小（压缩后）
# Output: 
#   1.SData_Integrated.RDS

    

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
SDataList <- readRDS(paste(data_path, "SDataList_Processed(QC3).RDS", sep = "/"))   #约5GB

options(future.globals.maxSize = 10 * 1024^3)
anchors <- FindIntegrationAnchors(
    SDataList,    #默认assay是SCT
    reduction = "cca",
    verbose = F)
IntegratedSData <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "SCT", 
    verbose = F)
str(IntegratedSData)
saveRDS(IntegratedSData, paste(data_path, "SData_Integrated.RDS", sep = "/"))

#注：本脚本运行6w细胞量时间约为50分钟。