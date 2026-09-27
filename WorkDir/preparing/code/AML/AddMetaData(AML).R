# 2025/8/25  添加了MetaData

# Input: SeuratData(Rename).RDS
# Output:
#   1.
#   



suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/AML")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/AML", sep = "")
data_path <- paste(file_path, "/data/AML", sep = "")
plot_path <- paste(file_path, "/plot/AML", sep = "")


#载入数据
SData <- readRDS(paste(data_path, "SeuratDataAML(Rename).RDS", sep = "/"))
Allmetadata <- data.table::fread("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/AML/GSE198052_RAW/GSM5936941_metadata_AML_Cells.txt")
#str(SData)  只有orig.ident\nCount_RNA\nFeature_RNA

#添加percent.MT
SData <- AddMetaData(SData, Allmetadata$percent.mt, col.name = "percent.MT")  #注意，如此添加metadata时需要保证txt文档中细胞的排列顺序与SData中一致
#添加PatientID
SData <- AddMetaData(SData, Allmetadata$patient, col.name = "PatientID")
#添加CancerType
SData@meta.data$CancerType <- rep("AML", ncol(SData))
#添加PreOrPost
