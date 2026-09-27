# 2025/8/13  load脚本

# Input: TNBC1数据 
# Output: 
    #1.QC by Original author Metric Distribution.pdf(包含nCount，nFeature，percent.HB，percent.MT四个metric)
    #2.SeuratData(Unrename).RDS



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/TNBC1")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/TNBC1", sep = "")
data_path <- paste(file_path, "/data/TNBC1", sep = "")
plot_path <- paste(file_path, "/plot/TNBC1", sep = "")



# 载入数据  （此处载入TNBC1数据）
RawData <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC1/10X", gene.column = 1)  #约8GB

SData <- CreateSeuratObject(RawData)


#检查线粒体基因、红细胞基因占比（判断该份数据是否已经QC）
SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")
RedCell.genes <- c("Hbb", "Hba", "Gypa", "Klf1", "Epor", "SLC4A1", "Hbb-bh1", "Pbx1") # 红细胞特异表达的基因
RedCell.genes <- CaseMatch(RedCell.genes, rownames(SData))
SData[["percent.HB"]] <- PercentageFeatureSet(SData, features = RedCell.genes)

ggplotdataframe <- data.frame(
    QC = rep(0, times = ncol(SData)), # QC=0表示还不确定是否原数据是否进行了QC
    nCount = SData@meta.data$nCount_RNA,
    nFeature = SData@meta.data$nFeature_RNA,
    percent.MT = SData@meta.data$percent.MT,
    percent.HB = SData@meta.data$percent.HB
)
#metric distribution可视化
p011 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nCount, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p012 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nFeature, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p013 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.MT, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p014 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.HB, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p01 <- p011 + p012 + p013 + p014 +
    patchwork::plot_layout(ncol = 4) +
    patchwork::plot_annotation(
        title = "QC by Original author Metrics Distribution", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/QC by Original author Metric Distribution.pdf", sep = ""),
    plot = p01, width = 18
)

#原作者已经进行过QC，分别针对count（600-120000）、feature（400-8000）和percent.MT（<=10%）三个指标来QC。

saveRDS(SData, file = paste(data_path, "/SeuratData(Unrename).RDS", sep = ""))