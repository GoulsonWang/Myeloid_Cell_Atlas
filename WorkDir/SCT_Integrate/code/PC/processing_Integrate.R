#2025/9/15   对整合后的数据再次进行预处理

# Input: SData_Integrated.RDS   #整合后的数据
# Output: 
#   1.Integrate_Data_Cluster_Plot.pdf  包含umap和Tsne两部分
#   2.SData_Integrated_Processed.RDS  
    

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")

#载入数据
SData <- readRDS(paste(data_path, "SData_Integrated.RDS", sep = "/"))

#预处理  active.assay: chr "integrated"
#SData <- ScaleData(SData, assay = "integrated", verbose = F)        #注：该行仅仅用做Demo
SData <- RunPCA(SData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    RunTSNE(dims = 1:30, verbose = F)
str(SData)     #结果保存为str(SData_Integrated_Processed).txt


#画图 umap+tsne
umapdata <- Embeddings(SData, reduction = "umap")
tsnedata <- Embeddings(SData, reduction = "tsne")
clusterdata <- SData@active.ident
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
        title = "PC Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)


saveRDS(SData, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))

#注：该脚本运行约14分钟