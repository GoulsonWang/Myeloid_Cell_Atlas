#2025/10/13     本脚本用与补充Integration_Total.R脚本的输出。主要将meta.data$Is.Doublet中的NA全部修改为Singlet



#input
#   1.SData_Integrated_Processed.RDS



#output
#   1.SData_Integrated_Processed_Supplement.RDS
#   2.Doublet_Distribution_Supplement.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/code/")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")

#载入数据
SData <- readRDS(file.path(data_path, "SData_Integrated_Processed.RDS"))



#将meta.data$Is_Doublet中的"NA"改为"Singlet"。
SData@meta.data$Is_Doublet[is.na(SData$Is_Doublet)] <- "Singlet"
p <- DimPlot(SData, reduction = "umap", group.by = "Is_Doublet")
ggsave(
    filename = paste(plot_path, "/Doublet_Distribution_Supplement.pdf", sep = ""),
    plot = p
)