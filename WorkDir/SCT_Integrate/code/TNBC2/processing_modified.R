#2025/9/16    处理processing.R的报错

#处理TNBC2 data，该data中包含了TNBC1和TNBC2数据中所有的Myeloid细胞，共56180。
#有详细的meta.data信息


#input：SDataList_Processed.RDS 
#output:
#   
#   1.SData_Integrated.RDS
#   2.SData_Integrated_Processed.RDS
#   3.DimPlot_Orig.ident-i.pdf(单个样本聚类后的UMAP图)
#   4.Integrate_Data_Cluster_Plot(MyeloidCell).pdf   多样本整合后聚类的UMAP、TSNE图



suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")


SDataListTransformed <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/data/TNBC2/SDataList_Processed.RDS")      #12GB
#class(SDataListTransformed)     #"list"


#Integrate
#此处尝试使用IntegrateLayers()函数，而非FindIntegrationAnchors+IntegrateData
options(future.globals.maxSize = 145 * 1024^3)
#此处SDataListTransformed中的每个单元的defaultAssays = "SCT"
#IntegratedSData <- IntegrateLayers(     #注：该函数的输入为一个Seurat对象，其中每个批次分别保存在一个layers下方。
#    SDataListTransformed, 
#    method = CCAIntegration,  
#    orig.reduction = "pca",
#    new.reduction = "integrated.cca",   #该reduction结果可直接用于聚类和可视化
#    verbose = FALSE
#)
#IntegratedSData[["RNA"]] <- JoinLayers(IntegratedSData[["RNA"]])

#去除SDataListTransformed中为Null的值
nlist <- length(SDataListTransformed)
SDataListTransformed <- SDataListTransformed[c(1 : 6, 8 : 21, 23 : 37, 39 : 62, 64 : nlist)]
class(SDataListTransformed[[7]])
class(SDataListTransformed[[22]])
class(SDataListTransformed[[38]])
class(SDataListTransformed[[63]])

anchors <- FindIntegrationAnchors(      #k.anchor参数越大，矫正力度越大
    object.list = SDataListTransformed,    #默认assay是SCT
    reduction = "cca",      #此处用rpca会更快一点
    verbose = F)
IntegratedSData <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "SCT", 
    k.weight = 50,          #默认值是100
    verbose = F)
IntegratedSData <- JoinLayers(IntegratedSData)      #内存不足，猜测原因是某几个样本细胞量太小。提示k.weight要小于56。
print("输出str(SData_Integrated)")
str(IntegratedSData)        #结果保存为str(SData_Integrated).txt
saveRDS(IntegratedSData, paste(data_path, "SData_Integrated.RDS", sep = "/"))



#processing again
#用IntegratedLayers()时，运行下列代码。
#DefaultAssay(IntegratedSData) <- "integrated"
#IntegratedSData_Processed <- FindNeighbors(IntegratedSData, reduction = "integrated.cca", verbose = F) %>%
#    FindClusters(resolution = 0.6, verbose = F) %>%
#    RunUMAP(dims = 1:30, reduction = "integrated.cca", verbose = F) %>%
#    RunTSNE(dims = 1:30, reduction = "integrated.cca", verbose = F)

#用FindIntegrationAnchors()+IntegratedSData()时，运行下列代码。
DefaultAssay(IntegratedSData) <- "integrated"
IntegratedSData_Processed <- RunPCA(IntegratedSData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    RunTSNE(dims = 1:30, verbose = F)
print("输出str(IntegratedSData_Processed)")
str(IntegratedSData_Processed)     #结果保存为str(IntegratedSData_Processed).txt
saveRDS(IntegratedSData_Processed, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))



#画图 umap+tsne
umapdata <- Embeddings(IntegratedSData_Processed, reduction = "umap")
tsnedata <- Embeddings(IntegratedSData_Processed, reduction = "tsne")
clusterdata <- IntegratedSData_Processed@meta.data$SubCluster
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata)

p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe, aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "TNBC2 Myeloid Cells Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot(MyeloidCell).pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)

