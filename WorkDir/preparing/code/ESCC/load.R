# 2025/8/25  load脚本

# Input: ESCC数据  
# Output: 
    #1.SeuratDataList(Raw).RDS   包含癌组织和癌旁组织两部分
    #2.Raw Data Metrics Distribution(orig.ident: i).pdf  样本的nCount，nFeature，percent.MT的小提琴图


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/ESCC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")

#载入数据(此处先载入前两个数据)
RawDatalist <- list()
RawPath <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/ESCC/Raw/RawData(decompression)/data"
SDatalist <- list()
for (i in 1:46){
    RawDatalist[[i]] <- Read10X(paste(RawPath, i, sep = ""))
    SDatalist[[i]] <- CreateSeuratObject(RawDatalist[[i]])
}

#class(RawDatalist[[2]])     #"Matrix"
#dim(RawDatalist[[2]])     #25026  8612
#head(rownames(RawDatalist[[2]]))    #"TSPAN6"   "DPM1"     "SCYL3"    "C1orf112" "FGR"      "CFH"  
#head(colnames(RawDatalist[[2]]))    #"AAACATCGAAACATCGCAAGACTA"   此处细胞名称为24位，比CellRanger的长！！
#class(SDatalist[[2]])       #"SeuratObject"
#dim(RawDatalist[[2]])       #25026  8612
#head(rownames(RawDatalist[[2]]))    #"TSPAN6"   "DPM1"     "SCYL3"    "C1orf112" "FGR"      "CFH" 
#head(colnames(RawDatalist[[2]]))    #"AAACATCGAAACATCGCAAGACTA"




#质量控制

#线粒体基因比例
for ( i in 1 : length(SDatalist)){
    SDatalist[[i]][["percent.MT"]] <- PercentageFeatureSet(SDatalist[[i]], pattern = "^MT-")
}

# 准备ggplot的data，还是将所有的data.frame收到一个list中
ggplotdatalist1 <- list()
for (i in 1:length(SDatalist)) {
    # 构建data.frame
    ggplotdatalist1[[i]] <- data.frame(
        QC = rep(0, times = ncol(SDatalist[[i]])), # QC=0表示是RawData，QC=1表示经过了QC
        nCount = SDatalist[[i]]@meta.data$nCount_RNA,
        nFeature = SDatalist[[i]]@meta.data$nFeature_RNA,
        percent.MT = SDatalist[[i]]@meta.data$percent.MT
    )

    # metric distribution可视化（ggplot编号为01）
    p011 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = nCount, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("0" = "Raw Data")) +
        guides(fill = FALSE) +
        theme_light()
    p012 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = nFeature, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("0" = "Raw Data")) +
        guides(fill = FALSE) +
        theme_light()
    p013 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = percent.MT, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("0" = "Raw Data")) +
        guides(fill = FALSE) +
        theme_light()
    p01 <- p011 + p012 + p013 +
        patchwork::plot_layout(ncol = 3) +
        patchwork::plot_annotation(
            title = paste("Raw Data Metrics Distribution(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/Raw Data Metric Distribution (orig.ident:", i, ").pdf", sep = ""),
        plot = p01, width = 18
    )
}

for (i in 1:46) {
    minfeature <- min(SDatalist[[i]]@meta.data$nFeature_RNA)
    cat("data", i, "的最小feature数为", minfeature, "\n")
}

for (i in 1:46){
    dimdata <- dim(SDatalist[[i]])
    cat("data", i, "的dim为", dimdata, "\n")
}

##此处的数据细胞最小feature个数均小于文章中给出的标准——200，故推断此数据未进行QC

saveRDS(SDatalist, file = paste(data_path, "SeuratDataList(Raw).RDS", sep = "/"))



