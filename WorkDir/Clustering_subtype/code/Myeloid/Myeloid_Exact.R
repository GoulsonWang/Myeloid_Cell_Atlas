#2025/10/30     
#本脚本用与提取注释结果中的Macrophage和Mast两大类，用与后续区分亚型
#注意：将所有数据Merge后进行预处理、聚类，未进行批次效应矫正



#input：
#   1.Description1.RDS

#output：
#   1.Myeloid_SData.RDS                 在注释文件种被注释为Macrophage和Mast的细胞  (所有count被merge得到)
#   2.SData_Merged_Processed.RDS        未经过批次效应矫正、处理后的RDS文件 （注意，对整体进行预处理，而非个体样本）
#   3.Batch_Identification.pdf          对6种潜在的批次效应来源变量的umap图



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
#SData_Overall <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/data/Description1.RDS")
#countsData <- LayerData(SData, layer = "counts")
#metaData <- SData@meta.data
#New_SData <- CreateSeuratObject(countsData, meta.data = metaData)
#str(New_SData)
#SData <- subset(New_SData, subset = Lineage %in% c("Macrophage", "Mast"))
#saveRDS(SData, file = file.path(data_path, "Myeloid_SData.RDS"))
#cat("输出str(Myeloid_SData)")
#str(SData)
#上述代码是提交至服务器运行的，直接输出了Myeloid_SData.RDS



SData <- readRDS(file.path(data_path, "Myeloid_SData.RDS"))       #2GB
table(SData$Lineage)
#Macrophage       Mast 
#     72148       5658 
#提取Demo
#Demo <- subset(SData, subset = orig.ident %in% c(1, 2))
#saveRDS(Demo, file = file.path(data_path, "Myeloid_SData(Demo).RDS"))



#预处理，先不处理批次效应，后续查看批次效应的主要来源是orig.ident还是CancerType,还是Patient.ID
dim(SData)  #44296   208
SData <- SCTransform(SData, vars.to.regress = "percent.MT", verbose = F) %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 1, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    RunTSNE(dims = 1:30, verbose = F)
#str(SData)
saveRDS(SData, file = file.path(data_path, "SData_Merged_Processed.RDS"))



#绘图
p1 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T) 
p2 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T, group.by = "CancerType") 
p3 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T, group.by = "orig.ident")
p4 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T, group.by = "Patient.ID")  
p5 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T, group.by = "PreOrPost") 
p6 = DimPlot(SData, reduction = "umap", pt.size = 0.1, label=T, group.by = "PrimaryOrMet") 
p <- p1 + p2 + p3 + p4 + p5 + p6 +
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(
        title = "Identify the origin of batch effect", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Batch_Identification.pdf", sep = ""),
    plot = p, 
    width = 18, 
    height = 12
)


