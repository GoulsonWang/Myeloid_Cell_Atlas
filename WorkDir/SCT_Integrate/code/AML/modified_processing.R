#2025/9/16      修正processing.R的bug


# Input: SData_Integrated.RDS   (AML数据)
# Output: 
#   1.SData_Integrated_Processed.RDS
#   2.Integrate_Data_Cluster_Plot.pdf



suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/AML", sep = "")
data_path <- paste(file_path, "/data/AML", sep = "")
plot_path <- paste(file_path, "/plot/AML", sep = "")

#
IntegratedSData <- readRDS(paste(data_path, "SData_Integrated.RDS", sep = "/"))




#processing again
IntegratedSData_Processed <- RunPCA(IntegratedSData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:50, n.neighbors = 40, verbose = F) %>%  #####################n.neighbors参数按需调整
    RunTSNE(dims = 1:50, verbose = F)
#print("输出str(IntegratedSData_Processed)")
#str(IntegratedSData_Processed)     #结果保存为str(IntegratedSData_Processed).txt
saveRDS(IntegratedSData_Processed, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))

#画图 umap+tsne
umapdata <- Embeddings(IntegratedSData_Processed, reduction = "umap")
tsnedata <- Embeddings(IntegratedSData_Processed, reduction = "tsne")
clusterdata <- IntegratedSData_Processed@active.ident
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata)

p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe, aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "AML Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot(Nneighbors40+npcs50).pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
