#2025/8/26  按照文章中写到的标准来QC

# Input: SeuratDatalist(Raw_Tumor).RDS   
# Output:
#       1.SeuratDataList(QC).RDS


suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/ESCC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")


# 载入数据
SDataList <- readRDS(paste(data_path, "SeuratDatalist(Raw_Tumor).RDS", sep = "/"))  #约3GB


#QC
QCList <- list()
cat("所有样本的最小feature数阈值为：200")
for ( i in 1 :length(SDataList)){
    #QC标准（细胞角度）
    min_feature <- 200
    max_feature <- quantile(SDataList[[i]]@meta.data$nFeature_RNA, 0.98)
    max_percent_MT <- 50
    max_count <- quantile(SDataList[[i]]@meta.data$nCount_RNA, 0.98)
    #QC
    QCList[[i]] <- subset(SDataList[[i]], 
        nFeature_RNA > min_feature & nFeature_RNA < max_feature &
        percent.MT < max_percent_MT & nCount_RNA < max_count)
    #QC标准（基因角度）
    expressed_cells_per_gene <- rowSums(GetAssayData(QCList[[i]], layer = "count"))
    QCList[[i]] <- subset(QCList[[i]], 
        features = names(expressed_cells_per_gene[expressed_cells_per_gene >= 5]))
    #检查结果
    minfeature_before <- min(SDataList[[i]]@meta.data$nFeature_RNA)
    minfeature_after <- min(QCList[[i]]@meta.data$nFeature_RNA)
    cat("QC前最小feature数为：", minfeature_before, "QC后为：", minfeature_after, "\n")
    maxfeature_before <- max(SDataList[[i]]@meta.data$nFeature_RNA)
    maxfeature_after <- max(QCList[[i]]@meta.data$nFeature_RNA)
    cat("第", i, "个样本的最大feature数阈值为：", max_feature, "\n")
    cat("QC前最大feature数为：", maxfeature_before, "QC后为：", maxfeature_after, "\n")
    maxcount_before <- max(SDataList[[i]]@meta.data$nCount_RNA)
    maxcount_after <- max(QCList[[i]]@meta.data$nCount_RNA)
    cat("第", i, "个样本的最大count数阈值为：", max_count, "\n")
    cat("QC前最大count数为：", maxcount_before, "QC后为：", maxcount_after, "\n")
    n_expressed_cells_per_gene_before <- min(expressed_cells_per_gene)
    n_expressed_cells_per_gene_after <- min(rowSums(GetAssayData(QCList[[i]], layer = "count"))) 
    cat("第", i, "个样本QC前基因在最少细胞中表达的数量为", n_expressed_cells_per_gene_before, 
        "QC后为", n_expressed_cells_per_gene_after, "\n")
}
#str(QCList)
saveRDS(QCList, paste(data_path, "SeuratDataList(QC).RDS", sep = "/"))



#可视化
ggplotdatalist1 <- list()
for (i in 1:length(QCList)) {
    # 构建data.frame
    ggplotdatalist1[[i]] <- data.frame(
        QC = rep(1, times = ncol(QCList[[i]])), # QC=0表示是RawData，QC=1表示经过了QC
        nCount = QCList[[i]]@meta.data$nCount_RNA,
        nFeature = QCList[[i]]@meta.data$nFeature_RNA,
        percent.MT = QCList[[i]]@meta.data$percent.MT
    )

    # metric distribution可视化（ggplot编号为01）
    p011 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = nCount, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC")) +
        guides(fill = FALSE) +
        theme_light()
    p012 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = nFeature, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC")) +
        guides(fill = FALSE) +
        theme_light()
    p013 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = percent.MT, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC")) +
        guides(fill = FALSE) +
        theme_light()
    p01 <- p011 + p012 + p013 +
        patchwork::plot_layout(ncol = 3) +
        patchwork::plot_annotation(
            title = paste("QC Data Metrics Distribution(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/QC Data Metric Distribution (orig.ident:", i, ").pdf", sep = ""),
        plot = p01, width = 18
    )
}