# 2025/8/14  QC1 (限制条件比较松)

# Input: SeuratDataList.RDS 
# Output: 
#   1.QC1SeuratDataList.RDS
#   2.QC1 Data Metric Distribution (orig.ident:i).pdf  第一次QC后的各样本metric分布小提琴图
    



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


#载入数据   

RawDataList <- readRDS(paste(data_path, "/SeuratDataList.RDS", sep = ""))  # about 2GB
#RawDataList <- list()   #测试数据
#RawDataList[[1]] <- RawData[[1]]
#RawDataList[[2]] <- RawData[[2]]



#第一次QC（主要针对percent.HB和percent.MT）
#设置QC阈值
maxMT <- c(15, 15, 15, 15, 15, 15, 15, 20, 25, 15, 10, 15, 40)
maxHB <- c(3, 3, 3, 3, 3, 25, 3, 3, 3, 25, 3, 3, 15)
QC1DataList <- list()
QC1ggplotDataList <- list()
#QC及可视化
for (i in 1 : length(RawDataList)){
    QC1DataList[[i]] <- subset(RawDataList[[i]], 
                               subset = percent.MT < maxMT[i] & percent.HB < maxHB[i])
    #准备QC后ggplot的data
    QC1ggplotDataList[[i]] <- data.frame(
      QC = rep(1, times = ncol(QC1DataList[[i]])),
      nCount = QC1DataList[[i]]@meta.data$nCount_RNA,
      nFeature = QC1DataList[[i]]@meta.data$nFeature_RNA,
      percent.MT = QC1DataList[[i]]@meta.data$percent.MT,
      percent.HB = QC1DataList[[i]]@meta.data$percent.HB
    )
    #画图  ggplot编号为11、12
    p011 <- ggplot(data = QC1ggplotDataList[[i]], aes(x = QC, y = nCount, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC1 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p012 <- ggplot(data = QC1ggplotDataList[[i]], aes(x = QC, y = nFeature, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC1 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p013 <- ggplot(data = QC1ggplotDataList[[i]], aes(x = QC, y = percent.MT, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC1 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p014 <- ggplot(data = QC1ggplotDataList[[i]], aes(x = QC, y = percent.HB, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("1" = "QC1 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p01 <- p011 + p012 + p013 + p014 +
        patchwork::plot_layout(ncol = 4) +
        patchwork::plot_annotation(
            title = paste("QC1 Data Metrics Distribution(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/QC1 Data Metric Distribution (orig.ident:", i, ").pdf", sep = ""),
        plot = p01, width = 18
    )
}

saveRDS(QC1DataList, file = paste(data_path, "QC1SeuratDataList.RDS", sep = "/"))


#orig,ident 13 样本percent.MT值高：10~30   疑似死细胞过多
#orig.ident 6 样本有少数几个细胞percent.HB在15~25
#orig.ident 10 样本有半数细胞percent.HB在5~20   混入了部分红细胞，聚类后观察是否高percent.HB的这些细胞聚集在一起。




