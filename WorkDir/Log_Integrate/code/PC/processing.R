#2025/9/25      processing,Integrate(RPCA),and processing again.



# Input: QC3SDataList.RDS   (PC数据)
# Output: 
#   1.SDataList_Transformed.RDS   (以样本为单位进行预处理)
#   2.DimPlot_Orig.ident-i.pdf (每个样本的聚类图)
#   3.SData_Integrated.RDS  (整合后的数据)
#   4.Integrate_Data_Cluster_Plot.pdf(整合后数据，再处理后的umap和tsne图)
#   5.SData_Integrated_Processed.RDS  (整合后、再处理后的数据)



#提交命令时候需要修改的地方：
#1.plan()
#2.SDataList <- readRDS()
#3.options()





suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")
#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 95 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")


#载入数据
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/PC/QC3SDataList.RDS")    #2GB



#processing
SDataListTransformed <- list()
nlist <- length(SDataList)
for(i in 1:nlist){
    SDataListTransformed[[i]] <- NormalizeData(SDataList[[i]], verbose = F) %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
        RunPCA(verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(resolution = 0.2, verbose = F) %>%
        RunUMAP(dims = 1:30, verbose = F)
    p1 <- DimPlot(SDataListTransformed[[i]], label = T)
    ggsave(paste(plot_path, "/DimPlot_Orig.ident-", i, ".pdf", sep = ""), plot = p1)
}
cat("\n", "输出str(SDataList_Transformed)")
#str(SDataListTransformed)   #结果保存为str(SDataList_Transformed).txt

saveRDS(SDataListTransformed, paste(data_path, "SDataList_Transformed.RDS", sep = "/"))



#Integration(RPCA)
anchors <- FindIntegrationAnchors(
    SDataListTransformed,    
    reduction = "rpca",
    verbose = F)
IntegratedSData <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "LogNormalize", 
    verbose = F)
DefaultAssay(IntegratedSData) <- "RNA"
IntegratedSData <- JoinLayers(IntegratedSData, assay = "RNA")
cat("\n", "输出str(SData_Integrated)")
str(IntegratedSData)        #结果保存为str(SData_Integrated).txt
saveRDS(IntegratedSData, paste(data_path, "SData_Integrated.RDS", sep = "/"))



#processing again
DefaultAssay(IntegratedSData) <- "integrated"
IntegratedSData_Processed <- ScaleData(IntegratedSData, verbose = F) %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    RunTSNE(dims = 1:30, verbose = F)
cat("\n", "输出str(SData_Integrated_Processed)")
str(IntegratedSData_Processed)     #结果保存为str(IntegratedSData_Processed).txt
saveRDS(IntegratedSData_Processed, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))



#画图 umap+tsne(整合后)
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
        title = "PC Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)



#整合前后对比
DefaultAssay(IntegratedSData) <- "RNA"
UnintegratedSData <- GetAssayData(IntegratedSData, assay = "RNA", layer = "counts")
metadata <- IntegratedSData@meta.data
UnintegratedSData <- CreateSeuratObject(counts = UnintegratedSData, 
    meta.data = metadata)
UnintegratedSData <- NormalizeData(UnintegratedSData, verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)
umapdata2 <- Embeddings(UnintegratedSData, reduction = "umap")
groupdata2 <- UnintegratedSData@meta.data$orig.ident
ggplotDataframe2 <- data.frame(umapdata, umapdata2, groupdata2)
p1 <- ggplot(ggplotDataframe2, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = groupdata2), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe2, aes(x = umap_1.1, y = umap_2.1)) +
        geom_point(aes(color = groupdata2), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "PC Integration after vs before ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Evaluate_Integration_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
