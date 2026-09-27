#Merge_Data.R脚本运行的补充，主要解决counts名称混乱


Total_SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/TotalData/Merge_Data.RDS")

Total_SData <- JoinLayers(Total_SData, assay = "RNA")
str(Total_SData)

Total_SData[["RNA"]] <- split(Total_SData[["RNA"]], f = Total_SData$orig.ident)     #count命名格式为count.[orig.ident]
str(Total_SData)        #结果保存为str(Merge_Data_Supplement).txt
saveRDS(Total_SData, file = "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/TotalData/Merge_Data_Supplement.RDS")

