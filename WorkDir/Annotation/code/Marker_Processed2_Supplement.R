#本脚本仅仅用与美化umap图



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




#Umap图
umapdata <- Embeddings(SData, reduction = "umap")
tsnedata <- Embeddings(SData, reduction = "tsne")
clusterdata <- SData@active.ident
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata)
label_df <- ggplotDataframe %>% # 计算cluster中心坐标
    group_by(clusterdata) %>%
    summarise(
        umap_1 = mean(umap_1, na.rm = TRUE),
        umap_2 = mean(umap_2, na.rm = TRUE)
    )
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
    geom_point(shape = 1, aes(color = clusterdata), size = 0.1, stroke = 0.1) +
    theme_light() +
    # 添加 cluster 标签
    geom_text(
        data = label_df,
        aes(label = clusterdata),
        size = 3,
        color = "black",
        vjust = -0.5, # 微调位置，避免覆盖点
        fontface = "bold"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) # 调整图例点大小
ggsave(
    filename = paste(plot_path, "/Marker_UMAP_Processed2.pdf", sep = ""),
    plot = p1,
    width = 7,
    height = 7
)
