#本脚本用于查看hdWGCNA的input的文件结构



suppressMessages(library(Seurat))


SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_Rename.RDS")


str(SData)
table(SData$Minor_Cell_Type)
