#2025/10/21     进行大类注释（共9类）,其中Macrophage所在大类有72148个细胞



# Input: 
#   1.SData_Marker_Processed2.RDS        Marker_Processed2.R的output

# Output: 
#   1.SData_Annotation1.RDS
#   2.Annotation1.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "SData_Marker_Processed2.RDS"))



#大类注释
SData@meta.data$Lineage <- "NA"
SData$Lineage[SData$seurat_clusters == "0"] <- "Epithelial"
SData$Lineage[SData$seurat_clusters == "1"] <- "T/NK"
SData$Lineage[SData$seurat_clusters == "2"] <- "Myeloid"
SData$Lineage[SData$seurat_clusters == "4"] <- "B"
SData$Lineage[SData$seurat_clusters == "5"] <- "Profilerating"
SData$Lineage[SData$seurat_clusters == "7"] <- "Fibroblast"
SData$Lineage[SData$seurat_clusters == "8"] <- "Endothelial"
SData$Lineage[SData$seurat_clusters == "9"] <- "Plasma"
SData$Lineage[SData$seurat_clusters == "10"] <- "Mast"
table(SData$Lineage)
#            B   Endothelial    Epithelial    Fibroblast    Macrophage 
#        30804         12193        109606         20098         72148 
#         Mast        Plasma Profilerating          T/NK 
#         5658         10535         30538        109970 



#saveRDS(SData, file = file.path(data_path, "SData_Annotation1.RDS"))



#再画一张UMAP图
umapdata <- Embeddings(SData, reduction = "umap")
clusterdata <- SData$Lineage
ggplotDataframe <- data.frame(umapdata, clusterdata)
label_df <- ggplotDataframe %>%         # 计算cluster中心坐标
    group_by(clusterdata) %>%
    summarise(
        umap_1 = mean(umap_1, na.rm = TRUE),
        umap_2 = mean(umap_2, na.rm = TRUE)
    )
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
    geom_point(shape = 1, aes(color = clusterdata), size = 0.1, stroke = 0.1) +     #shape=1指绘图元素为空心圆圈， stroke指圆圈边厚度
    theme_light() +
    # 添加 cluster 标签
    geom_text(
        data = label_df,
        aes(label = clusterdata),
        size = 3,
        color = "black",
        vjust = -0.5,       # 微调位置，避免覆盖点
        fontface = "bold"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) # 调整图例点大小
#p1 <- p1 +
#    patchwork::plot_annotation(
#        title = "Lineage Annotation Result ",
#        theme = theme(plot.title = element_text(hjust = 0.5))
#    )
ggsave(
    filename = paste(plot_path, "/Annotation1.pdf", sep = ""),
    plot = p1,
    width = 8,
    height = 6
)
ggsave(
    filename = paste(plot_path, "/Annotation1.svg", sep = ""),
    plot = p1,
    width = 8,
    height = 6
)