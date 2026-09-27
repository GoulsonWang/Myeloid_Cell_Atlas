#2025/10/10     本脚本用与运行DoubletFinder



#input:     TNBC1数据在QC时发现数据已经经过了清洗，为进一步确定是否是双细胞去除后的数据，检测
#   1.SeuratData(Rename).RDS    TNBC1数据集经过preparing的output

#output:    
#   1.Doublet_Distribution-origi-Rate006.pdf
#   2.SDataList.RDS


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(future))
suppressMessages(library(DoubletFinder))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/code/")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder"
code_path <- paste(file_path, "/code/TNBC1", sep = "")
data_path <- paste(file_path, "/data/TNBC1", sep = "")
plot_path <- paste(file_path, "/plot/TNBC1", sep = "")



#载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC1/SeuratData(Rename).RDS")
#table(SData@meta.data$orig.ident)
cells_to_keep <- colnames(SData)[!SData$orig.ident %in% c(16, 29)]
SData_Subset <- subset(SData, cells = cells_to_keep)  #此处筛除两个样本16和29，原因：细胞量小于200
#table(SData_Subset@meta.data$orig.ident)
SDataList <- SplitObject(SData_Subset, split.by = "orig.ident")
#preparing
for( i in 1:length(SDataList)){
    SDataList[[i]] <- SDataList[[i]] %>%
        NormalizeData(verbose = F) %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
        RunPCA(verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(resolution = 0.1, verbose = F) %>%
        RunUMAP(dims = 1:30, reduction = "pca", verbose = F) %>%
        RunTSNE(dims = 1:30, reduction = "pca", verbose = F)

}



#DoubletFinder      注意：DoubletFinder函数运行要求一个样本为一个Seurat对象，多样本处理时利用list合并，再循环运行
sweep.res.list_kidney <- list()
sweep.stats_kidney <- list()
bcmvn_kidney <- list()
pk.best <- list()
for(i in 1:length(SDataList)){
    #估算在PC空间中要考虑的最近邻居个数pK
    sweep.res.list_kidney[[i]] <- paramSweep(
        SDataList[[i]], 
        PCs = 1:10, 
        sct = FALSE, 
#        num.cores = 10
        )
    sweep.stats_kidney[[i]] <- summarizeSweep(
        sweep.res.list_kidney[[i]], 
        GT = FALSE
        )
    bcmvn_kidney[[i]] <- find.pK(sweep.stats_kidney[[i]])
    pk.best[[i]] <- bcmvn_kidney[[i]] %>%
        dplyr::arrange(desc(BCmetric)) %>%
        dplyr::pull(pK) %>%
        .[1] %>% as.character() %>% as.numeric()
    #同源双细胞率估计（最终这部分细胞不纳入计算）
    annotations <- SDataList[[i]]@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.06*nrow(SDataList[[i]]@meta.data))       ## 此处假设双细胞率为0.02，此值一般为0.01-0.1之间，建议多次尝试
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    #Run DoubletFinder
    SDataList[[i]] <- doubletFinder(
        SDataList[[i]], 
        PCs = 1:10, 
        pN = 0.25,      #产生人工双细胞比例
        pK = pk.best[[i]],      
        nExp = nExp_poi.adj,        #预期双细胞个数
        reuse.pANN = F, 
        sct = FALSE)
    #修改meta.data列名
    colnames(SDataList[[i]]@meta.data)[length(colnames(SDataList[[i]]@meta.data))-1] <- "Doublet_Score"
    colnames(SDataList[[i]]@meta.data)[length(colnames(SDataList[[i]]@meta.data))] <- "Is_Doublet"
    #绘制Doublet分布UMAP图
    p1 <- DimPlot(SDataList[[i]], reduction = "umap", group.by = "Is_Doublet")
    p2 <- DimPlot(SDataList[[i]], reduction = "umap", group.by = "seurat_clusters", label = T)
    p <- p1 + p2 +
        patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(
            title = "Doublet Distribution ", 
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/Doublet_Distribution-orig", i, "-Rate006", ".pdf", sep = ""),
        plot = p, 
        width = 16, 
        height = 7
    )
    #print(p)
}
saveRDS(SDataList, file = file.path(data_path, "SDataList.RDS"))


#对数据进行筛选。筛选原则：观察Doublet是否高度集中于某一簇中、或表现为离群值，如果是，则删除该簇或离群值。


#对Demo（含3个样本）分析后发现，2各样本的Doublet存在成簇现象，应当过滤会比较好。故对所有TNBC1样本进行DoubletFinder
#担忧：同时进行DoubletFinder和Scrublet是否合理？