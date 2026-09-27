#本脚本用于提取出SData_Annatation1.RDS中的B,Plasma细胞的数据发给师弟
#2026/5/20


library(Seurat)
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/data/SData_Annotation1.RDS")
colnames(SData@meta.data)
Demo <- subset(SData, subset = Lineage  %in% c("B", "Plasma"))
dim(Demo)
str(Demo)
saveRDS(Demo, file = "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/data/B_Plasma.RDS")