#2025/10/11     本脚本用与将四个数据集都经过DoubletFinder检测后再Merge、Integration



# Input: 
#   1.PC        运行了DoubletFinder的output
#   2.TNBC1     运行了DoubletFinder的output
#   3.ESCC      运行了DoubletFinder的output
#   4.NSCLC     运行了DoubletFinder的output
#   注意：此处的input都是DoubletFinder的结果

# Output: 
#   1.Merge_SData2.RDS     合并后的结果



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))



#载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")

#载入数据
ESCC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/ESCC/SDataList.RDS")   
NSCLC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/NSCLC/SDataList.RDS")   
PC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/PC/SDataList.RDS")
#注意：PC数据集中orig.ident=4的样本细胞只有177个， orig.ident=3的样本细胞只有553
TNBC1 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/TNBC1/SDataList.RDS") 



#添加meta.data信息：肿瘤类别
# 合并所有 list
CancerData <- c(ESCC, NSCLC, PC, TNBC1)
cancer_types <- rep(c("ESCC", "NSCLC", "PC", "TNBC1"), 
                    times = c(length(ESCC), length(NSCLC), length(PC), length(TNBC1)))
for (i in seq_along(CancerData)) {       #seq_along(CancerData)等同于1：length（CancerData），优点在于当CancerData为0时前者返回integer(0)，后者返回c(1, 0)
    obj <- CancerData[[i]]
    cancer <- cancer_types[i]
    cell_count <- ncol(obj)
    cancer_metadata <- rep(cancer, times = cell_count)
    CancerData[[i]] <- AddMetaData(obj, metadata = cancer_metadata, col.name = "CancerType")
    #cat(table(CancerData[[i]]@meta.data$CancerType))
}



# 合并所有对象
Total_SData <- merge(CancerData[[1]], CancerData[2:length(CancerData)], merge.data = F)
cat("merge成功！")



#重新整合，处理counts名称混乱问题
Total_SData <- JoinLayers(Total_SData, assay = "RNA")
Merged_Count <- LayerData(Total_SData, layer = "counts", assay = "RNA")
Merged_Metadata <- Total_SData@meta.data
SData_Merged <- CreateSeuratObject(
    counts = Merged_Count, 
    meta.data = Merged_Metadata
)
cells_to_keep <- colnames(SData_Merged)[!SData_Merged$orig.ident %in% c(4, 16, 29)]
SData_Merged <- subset(SData_Merged, cells = cells_to_keep)       #此处删除小于200的样本, orig,ident序号为4，另外16，29已经在DoubletFinder时剔除
SData_Merged[["RNA"]] <- split(SData_Merged[["RNA"]], f = SData_Merged$orig.ident)     #count命名格式为count.[orig.ident]
cat("count名称混乱问题已解决！")



# 保存
cat("输出str(Merge_Data)")
str(SData_Merged)    #结果保存为str(Merge_Data).txt
saveRDS(SData_Merged, file = file.path(data_path, "Merge_Data2.RDS"))
cat("Merge的输出数据已保存！")