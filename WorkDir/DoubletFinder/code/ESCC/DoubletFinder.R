#2025/10/10     本脚本用与运行DoubletFinder



#input:     本步骤用与单样本QC、聚类、umap之后，多样本整合之前
#   1.SeuratDataList(QC).RDS    preparing的output

#output:    
#   1.Doublet_Distribution-origi.pdf    


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(future))
suppressMessages(library(DoubletFinder))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/code/")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")

#载入数据
#SData_Merged <- readRDS(file.path(data_path, "Merge_Data_Supplement.RDS"))
#注意：PC数据集中orig.ident=4的样本细胞只有177个，故筛除。剩余68个样本
#SData_Merged <- SplitObject(SData_Merged, split.by = "orig.ident")
#Demo <- list()
#Demo[[1]] <- SData_Merged[[1]]
#Demo[[2]] <- SData_Merged[[2]]



#载入数据
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/ESCC/SeuratDataList(QC).RDS")
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
        num.cores = 10
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
    nExp_poi <- round(0.06*nrow(SDataList[[i]]@meta.data))       ## 此处假设双细胞率为0.06，此值一般为0.01-0.1之间，建议多次尝试
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
        filename = paste(plot_path, "/Doublet_Distribution-orig", i, ".pdf", sep = ""),
        plot = p, 
        width = 16, 
        height = 7
    )
}



saveRDS(SDataList, file = file.path(data_path, "SDataList.RDS"))
#对数据进行筛选。筛选原则：观察Doublet是否高度集中于某一簇中、或表现为离群值，如果是，则删除该簇或离群值。
