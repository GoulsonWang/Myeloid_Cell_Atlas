#2025/9/22      本脚本用与对比Seurat@assay$integrate@assay.orig下的值分别为RNA和SCT时的不同

# Input: 
# Output: 
#   1.
#   2.
#   3.
#   4.
#   5.


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")

#取Demo
#SDataListTransformed <- readRDS(paste(data_path, "SDataList_Processed.RDS", sep = "/"))
#nlist <- length(SDataListTransformed)
#SDataListTransformed_Demo <- list()
#for(i in 1:3){      #取列表前三个样本为Demo， 约1.5w个细胞
#    SDataListTransformed_Demo[[i]] <- SDataListTransformed[[i]]
#}
#saveRDS(SDataListTransformed_Demo, file = paste(data_path, "SDataList_Processed(Demo).RDS", sep = "/"))
#str(SDataListTransformed_Demo)  #结果保存为str(SDataListTransformed(Demo)).txt



#载入数据
SDataListTransformed <- readRDS(paste(data_path, "SDataList_Processed(Demo).RDS", sep = "/"))   #1GB
class(SDataListTransformed)

#利用旧代码来整合
DefaultAssay(SDataListTransformed[[1]])
options(future.globals.maxSize = 2 * 1024^3)
anchors <- FindIntegrationAnchors(          #运行2分钟
    SDataListTransformed,    #默认assay是SCT
    reduction = "cca",
    verbose = F)
IntegratedSData <- IntegrateData(           #有Warning：Layer counts isn't present in the assay object; returning NULL
    anchorset = anchors, 
    normalization.method = "SCT", 
    verbose = F)
IntegratedSData@assays$integrate@assay.orig    #此处返回的是RNA
str(IntegratedSData)    #结果保存为str(test1_old)




#利用SCT来整合
nlist <- length(SDataListTransformed)
for ( i in 1:nlist){
    DefaultAssay(SDataListTransformed[[i]]) <- "SCT"
}
anchors <- FindIntegrationAnchors(      #运行时间约3min
    SDataListTransformed,    #默认assay是SCT
    reduction = "cca",
    verbose = T)
IntegratedSData <- IntegrateData(       
    anchorset = anchors,  
    verbose = T)
IntegratedSData@assays$integrate@assay.orig    #此处返回的是NULL
str(IntegratedSData)    #结果保存为str(test1_SCT).txt

#processing again
DefaultAssay(IntegratedSData) <- "integrated"   
#注意：此时IntegratedSData@assay@integrated$scale.data为空，直接运行RunPCA()会报错。Error in data.use[features, ] : no 'dimnames' attribute for array
#为什么呢？？
IntegratedSData_Processed <- RunPCA(IntegratedSData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)

pSCT <- DimPlot(IntegratedSData_Processed, reduction = )



#总结：IntegrateData()函数中，normalization.method = "SCT"参数不可或缺，其显式的说明了标准化所用方法，如果缺失这个参数，则可能会导致判断错误，在RunPCA()时候出错。