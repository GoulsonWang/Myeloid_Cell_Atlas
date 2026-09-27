# 2025/9/12  QC3 预处理画UMAP图之后，发现8，9，10，13号样本效果较差，故重新进行QC，主要收紧percent.MT和HB的阈值

# Input: QC2SeuratDataList.RDS 
# Output: 
#   1.
#   2.
    



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/PC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")



#载入数据
SDataList <- readRDS(paste(data_path, "/QC2SeuratDataList.RDS", sep = ""))

#设置QC阈值
maxMT <- 15
maxHB <- 3


QCSDataList <- list()
for(i in 1:length(SDataList)){
    QCSDataList[[i]] <- subset(SDataList[[i]],
                               subset = percent.MT < maxMT & percent.HB < maxHB)
}

saveRDS(QCSDataList, file = paste(data_path, "/QC3SDataList.RDS", sep = ""))
