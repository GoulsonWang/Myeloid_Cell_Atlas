#2025/10/27     本脚本用与将SData_Annotation1.RDS（Annotation1.R的output）中TNBC1部分提取出来，用与Validation_Annotation.R



# Input: 
#   1.SData_Annotation1.RDS      Annotation1.R的output

# Output: 
#   1.Annotated_TNBC1.RDS
#   



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code/Validation", sep = "")
data_path <- paste(file_path, "/data/Validation", sep = "")
plot_path <- paste(file_path, "/plot/Validation", sep = "")




SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/data/SData_Annotation1.RDS")
#TNBC1 <- subset(SData, subset = orig.ident %in% c(14:29))
#saveRDS(TNBC1, file = file.path(data_path, "Annotated_TNBC1.RDS"))

NSCLC <- subset(SData, subset = orig.ident %in% c(79:93))
table(NSCLC$orig.ident)
saveRDS(NSCLC, file = file.path(data_path, "Annotated_NSCLC.RDS"))