#2025/9/14      processing,Integrate,and processing again.


# Input: SeuratDataAML(Rename)_AddPercentMT.RDS   (AML数据)
# Output: 
#   1.SDataList_Processed.RDS   (以样本为单位进行预处理)
#   2.DimPlot_Orig.ident-i.pdf (每个样本的聚类图)
#   3.SData_Integrated.RDS  (整合后的数据)
#   4.Integrate_Data_Cluster_Plot.pdf(整合后数据，再处理后的umap和tsne图)
#   5.SData_Integrated_Processed.RDS  (整合后、再处理后的数据)



#关键修改参数：RunUMAP(dims = 1:50， FindClusters(resolution = 0.6， RunUMAP(dims = 1:50
#RunTSNE(dims = 1:50， FindClusters(resolution = 0.6， RunUMAP(dims = 1:30, n.neighbors = 50.

suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/AML", sep = "")
data_path <- paste(file_path, "/data/AML", sep = "")
plot_path <- paste(file_path, "/plot/AML", sep = "")

#载入数据
SData <- readRDS(paste(data_path, "SeuratDataAML(Rename)_AddPercentMT.RDS", sep = "/"))    #2GB


#按样本拆分成列表
SDataList <- SplitObject(SData, split.by = "orig.ident")

#SCTransform
SDataListTransformed <- list()
nlist <- length(SDataList)
options(future.globals.maxSize = 60 * 1024^3)   #提高内存保护机制阈值
for( i in 1:nlist){
    SDataListTransformed[[i]] <- SCTransform(SDataList[[i]], vars.to.regress = "percent.MT", verbose = FALSE)
}

#聚类画图
for(i in 1:nlist){
    SDataListTransformed[[i]] <- RunPCA(SDataListTransformed[[i]], verbose = F) %>% 
        RunUMAP(dims = 1:50, verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(verbose = F)
    p1 <- DimPlot(SDataListTransformed[[i]], label = T)
    ggsave(paste(plot_path, "/DimPlot_Orig.ident-", i, ".pdf", sep = ""), plot = p1)
}

print("输出str(SDataListTransformed)")
str(SDataListTransformed)      #结果保存为str(SDataListTransformed).txt
saveRDS(SDataListTransformed, file = paste(data_path, "/SDataList_Processed.RDS", sep = ""))



#Integrate
options(future.globals.maxSize = 60 * 1024^3)
anchors <- FindIntegrationAnchors(
    SDataListTransformed,    #默认assay是SCT
    reduction = "cca",
    verbose = F)
IntegratedSData <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "SCT", 
    verbose = F)
#print("输出str(SData_Integrated)")
#str(IntegratedSData)        #结果保存为str(SData_Integrated).txt
saveRDS(IntegratedSData, paste(data_path, "SData_Integrated.RDS", sep = "/"))



#processing again
IntegratedSData_Processed <- RunPCA(IntegratedSData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:50, verbose = F) %>%
    RunTSNE(dims = 1:50, verbose = F)
#print("输出str(IntegratedSData_Processed)")
#str(IntegratedSData_Processed)     #结果保存为str(IntegratedSData_Processed).txt
saveRDS(IntegratedSData_Processed, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))

#画图 umap+tsne
umapdata <- Embeddings(IntegratedSData_Processed, reduction = "umap")
tsnedata <- Embeddings(IntegratedSData_Processed, reduction = "tsne")
clusterdata <- IntegratedSData_Processed@active.ident
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
        title = "AML Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
