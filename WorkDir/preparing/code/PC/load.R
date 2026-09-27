# 2025/8/13  load脚本

# Input: PC数据 
# Output: 
    #1.Raw Data Metric Distribution (orig.ident:--).pdf(包含nCount，nFeature，percent.HB，percent.MT四个metric)
    #2.Raw Data Metric Correlation (orig.ident:--).pdf（包含nCount分别与nFeature，percent.HB，percent.MT三个变量的散点图）
    #3.SeuratDataList.RDS（包含所有样本的Seurat对象）(未合并)



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/PC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")



# 载入数据  （此处先载入PC数据，sample ID 07 08 结尾）
RawDatalist <- list()
RawDatalist[[1]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279204_01-01B")
RawDatalist[[2]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279205_01-01T")
RawDatalist[[3]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279206_01-04B")
RawDatalist[[4]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279207_01-05B")
RawDatalist[[5]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279208_01-06T")
RawDatalist[[6]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279209_02-01B")
RawDatalist[[7]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279210_02-01TL")
RawDatalist[[8]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279211_02-01TP")
RawDatalist[[9]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279212_02-02B")
RawDatalist[[10]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279213_02-02T")
RawDatalist[[11]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279214_02-06B")
RawDatalist[[12]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279215_02-07B")
RawDatalist[[13]] <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279216_02-07T")


# 创建Seurat对象，收到一个list中
SDatalist <- list()
for (i in 1:length(RawDatalist)) {
    SDatalist[[i]] <- CreateSeuratObject(RawDatalist[[i]])
    colnames(SDatalist[[i]]) <- gsub("1", i, colnames(SDatalist[[i]])) # 改细胞名，以orig.ident结尾
    SDatalist[[i]] <- AddMetaData(SDatalist[[i]],   #添加orig.ident至meta.data中
        metadata = stringr::str_split(colnames(SDatalist[[i]]), "-", simplify = T)[, 2], # 用str_split将colnames分为两部分，选择后半段为orig.ident
        col.name = "orig.ident"
    )
}


# 检查
#table(SDatalist[[1]]@meta.data$orig.ident)  #无样本名称，需要在metadata 中加orig.ident一列
#table(SDatalist[[2]]@meta.data$orig.ident)
#dim(SDatalist[[1]])   #36601   270   只有270个细胞
#dim(SDatalist[[2]])   #36601   725   
#str(SData)
#head(rownames(SDatalist[[2]]))   #"MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"
#head(colnames(SDatalist[[2]]))    #"AAAGCAAAGAAACGAG-1" "AAAGCAACACATTCGA-1" "AAAGCAAGTCGAGTTT-1"
#table(SDatalist[[2]]@meta.data$orig.ident)



# 质量控制

# 线粒体和红细胞基因比例
RedCell.genes <- c("Hbb", "Hba", "Gypa", "Klf1", "Epor", "SLC4A1", "Hbb-bh1", "Pbx1") # 红细胞特异表达的基因
for ( i in 1 : length(SDatalist)){
    SDatalist[[i]][["percent.MT"]] <- PercentageFeatureSet(SDatalist[[i]], pattern = "^MT-")
    RedCell.genes <- CaseMatch(RedCell.genes, rownames(SDatalist[[i]]))
    SDatalist[[i]][["percent.HB"]] <- PercentageFeatureSet(SDatalist[[i]], features = RedCell.genes)
}


# 准备ggplot的data，还是将所有的data.frame收到一个list中
ggplotdatalist1 <- list()
for (i in 1:length(SDatalist)) {
    # 构建data.frame
    ggplotdatalist1[[i]] <- data.frame(
        QC = rep(0, times = ncol(SDatalist[[i]])), # QC=0表示是RawData，QC=1表示经过了QC
        nCount = SDatalist[[i]]@meta.data$nCount_RNA,
        nFeature = SDatalist[[i]]@meta.data$nFeature_RNA,
        percent.MT = SDatalist[[i]]@meta.data$percent.MT,
        percent.HB = SDatalist[[i]]@meta.data$percent.HB
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
    p014 <- ggplot(data = ggplotdatalist1[[i]], aes(x = QC, y = percent.HB, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("0" = "Raw Data")) +
        guides(fill = FALSE) +
        theme_light()
    p01 <- p011 + p012 + p013 + p014 +
        patchwork::plot_layout(ncol = 4) +
        patchwork::plot_annotation(
            title = paste("Raw Data Metrics Distribution(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/Raw Data Metric Distribution (orig.ident:", i, ").pdf", sep = ""),
        plot = p01, width = 18
    )

    # 散点图可视化(ggplot编号为02）
    p021 <- ggplot(data = ggplotdatalist1[[i]], aes(x = nCount, y = percent.MT)) +
        geom_point(size = 0.05) +
        labs(title = "nCount ~ percent.MT") +
        theme_light()
    p022 <- ggplot(data = ggplotdatalist1[[i]], aes(x = nCount, y = percent.HB)) +
        geom_point(size = 0.05) +
        labs(title = "nCount ~ percent.HB") +
        theme_light()
    p023 <- ggplot(data = ggplotdatalist1[[i]], aes(x = nCount, y = nFeature)) +
        geom_point(size = 0.05) +
        labs(title = "nCount ~ nFeature") +
        theme_light()
    p02 <- p021 + p022 + p023 +
        patchwork::plot_annotation(
            title = paste("Raw Data Metrics Correlation(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/Raw Data Metric Correlation (orig.ident:", i, ").pdf", sep = ""),
        plot = p02, width = 21
    )
}

saveRDS(SDatalist, file = paste(data_path, "SeuratDataList.RDS", sep = "/"))