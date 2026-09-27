#2025/9/18      本脚本用与检查processing_modified.R的报错原因

#报错原因：内存不足（申请了160G）。猜测原因：存在某些样本的细胞量太小。

#input：SDataList_Processed.RDS 
#output:
#   
#   1.
#   2.
#   3.
#   4.



suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")

SDataListTransformed <- readRDS(paste(data_path, "SDataList_Processed.RDS", sep = "/"))  #12GB

#检查dim
for ( i in 1:length(SDataListTransformed)){
    cat("第", i, "个样本的dim值为：", dim(SDataListTransformed[[i]]), "\n")
}
#74个样本， 5个样本细胞量<100,  11个样本细胞量<200
