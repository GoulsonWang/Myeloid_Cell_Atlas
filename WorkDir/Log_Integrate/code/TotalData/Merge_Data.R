#2025/10/8  合并多个Seurat对象为一个Seurat对象（不校正）



# Input: 
#   1.PC    
#   2.TNBC1
#   3.ESCC
#   4.NSCLC
#   注意：此处的input是preparing的结果，而非数据集整合的结果

# Output: 
#   1.Merge_SData.RDS     合并后的结果



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))



#载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")

#载入数据
ESCC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/ESCC/SeuratDataList(QC).RDS")   
NSCLC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/NSCLC/SeuratData(Rename).RDS")   
NSCLC <- SplitObject(NSCLC, split.by = "orig.ident")
PC <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/PC/QC3SDataList.RDS")
#注意：PC数据集中orig.ident=4的样本细胞只有177个， orig.ident=3的样本细胞只有553
TNBC1 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC1/SeuratData(Rename).RDS") 
TNBC1 <- subset(TNBC1, orig.ident != 16 & orig.ident != 29)     #此处筛除两个样本（原因：样本含细胞太少，分别为8和147）
TNBC1 <- SplitObject(TNBC1, split.by = "orig.ident")



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
Total_SData <- CancerData[[1]]

for (i in 2:length(CancerData)) {
  Total_SData <- merge(
    x = Total_SData,
    y = CancerData[[i]]
  )
}



# 保存
str(Total_SData)    #结果保存为str(Merge_Data).txt
saveRDS(Total_SData, file = file.path(data_path, "Merge_Data.RDS"))
cat("Merge complete. Total cells:", ncol(Total_SData), "\n")
