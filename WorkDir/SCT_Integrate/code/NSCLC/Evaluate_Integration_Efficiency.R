#2025/9/22     画图对比整合前后效果

# Input: SData_Integrated_Processed.RDS   (NSCLC数据)
# Output: 
#   1.Evaluate_Integration_Efficiency.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/NSCLC", sep = "")
data_path <- paste(file_path, "/data/NSCLC", sep = "")
plot_path <- paste(file_path, "/plot/NSCLC", sep = "")

#载入数据
IntegratedSData_Processed <- readRDS(paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))    #500MB
#str(IntegratedSData_Processed)
#names(IntegratedSData_Processed@meta.data)


#umap图：整合前后按样本染色
#整合前
IntegratedSData_Processed_before <- JoinLayers(IntegratedSData_Processed, assay = "RNA")
#names(IntegratedSData_Processed_before@assays$RNA@layers)      #"counts"
#length(IntegratedSData_Processed_before@assays$RNA@layers)     # 1

DefaultAssay(IntegratedSData_Processed_before) <- "RNA"
options(future.globals.maxSize = 100 * 1024^3)
IntegratedSData_Processed_before <- SCTransform(IntegratedSData_Processed_before)
#str(IntegratedSData_Processed_before)
DefaultAssay(IntegratedSData_Processed_before) <- "SCT"
IntegratedSData_Processed_before <- RunPCA(IntegratedSData_Processed_before, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) 
umapdata_before <- Embeddings(IntegratedSData_Processed_before, reduction = "umap")

#整合后
umapdata_after <- Embeddings(IntegratedSData_Processed, reduction = "umap")
origident <- IntegratedSData_Processed@meta.data$orig.ident
ggplotDataframe <- data.frame(umapdata_after, umapdata_before, origident)   
#colnames(ggplotDataframe)   #"umap_1"    "umap_2"    "umap_1.1"  "umap_2.1"  "origident"
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = origident), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 1))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe, aes(x = umap_1.1, y = umap_2.1)) +
        geom_point(aes(color = origident), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 1))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "NSCLC Integrated After vs Before ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Evaluate_Integration_Efficiency.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
