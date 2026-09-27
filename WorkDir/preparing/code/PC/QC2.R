#2025/8/15      QC2      主要针对count数过多的细胞（multiblets）

# Input: QC1SeuratDataList.RDS 
# Output: 
#   1.QC2SeuratDataList.RDS
#   2.QC2 Data Metric Distribution (orig.ident:i).pdf  第二次QC后的各样本metric分布小提琴图
 


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
QC1DataList <- readRDS(paste(data_path, "QC1SeuratDataList.RDS", sep = "/"))


#第二次QC（主要针对nCount值）
maxcount <- c(13000, 11000, 30000, 15000, 15000, 13000, 45000, 35000, 15000, 8000, 13000, 5000, 6000)
QC2DataList <- list()
QC2ggplotDataList <- list()
#QC及可视化
for (i in 1 : length(QC1DataList)){
    QC2DataList[[i]] <- subset(QC1DataList[[i]], 
                               subset = nCount_RNA < maxcount[i])
    #准备QC后ggplot的data
    QC2ggplotDataList[[i]] <- data.frame(
      QC = rep(2, times = ncol(QC2DataList[[i]])),
      nCount = QC2DataList[[i]]@meta.data$nCount_RNA,
      nFeature = QC2DataList[[i]]@meta.data$nFeature_RNA,
      percent.MT = QC2DataList[[i]]@meta.data$percent.MT,
      percent.HB = QC2DataList[[i]]@meta.data$percent.HB
    )
    #画图  ggplot编号为21、22
    p211 <- ggplot(data = QC2ggplotDataList[[i]], aes(x = QC, y = nCount, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("2" = "QC2 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p212 <- ggplot(data = QC2ggplotDataList[[i]], aes(x = QC, y = nFeature, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("2" = "QC2 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p213 <- ggplot(data = QC2ggplotDataList[[i]], aes(x = QC, y = percent.MT, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("2" = "QC2 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p214 <- ggplot(data = QC2ggplotDataList[[i]], aes(x = QC, y = percent.HB, fill = QC)) +
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_x_discrete(labels = c("2" = "QC2 Data")) +
        guides(fill = FALSE) +
        theme_light()
    p21 <- p211 + p212 + p213 + p214 +
        patchwork::plot_layout(ncol = 4) +
        patchwork::plot_annotation(
            title = paste("QC2 Data Metrics Distribution(orig.ident:", i, ")", sep = ""),
            theme = theme(plot.title = element_text(hjust = 0.5))
        )
    ggsave(
        filename = paste(plot_path, "/QC2 Data Metric Distribution (orig.ident:", i, ").pdf", sep = ""),
        plot = p21, width = 18
    )
}

saveRDS(QC2DataList, file = paste(data_path, "QC2SeuratDataList.RDS", sep = "/"))
