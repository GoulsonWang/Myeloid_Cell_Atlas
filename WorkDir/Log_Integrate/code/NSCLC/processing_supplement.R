#2025/9/25      对processing.R的补充（重新绘制整合前后对比umap图）



# Input: 
# Output: 
#   1.Evaluate_Integration_Plot.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(future))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")
#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 95 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/NSCLC", sep = "")
data_path <- paste(file_path, "/data/NSCLC", sep = "")
plot_path <- paste(file_path, "/plot/NSCLC", sep = "")


#载入数据
IntegratedSData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/NSCLC/SData_Integrated.RDS")
IntegratedSData_Processed <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/data/NSCLC/SData_Integrated_Processed.RDS")



#整合前后对比
DefaultAssay(IntegratedSData) <- "RNA"
UnintegratedSData <- GetAssayData(IntegratedSData, assay = "RNA", layer = "counts")
metadata <- IntegratedSData@meta.data
UnintegratedSData <- CreateSeuratObject(counts = UnintegratedSData, 
    meta.data = metadata)
UnintegratedSData <- NormalizeData(UnintegratedSData, verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)
umapdata2 <- Embeddings(UnintegratedSData, reduction = "umap")
groupdata2 <- UnintegratedSData@meta.data$orig.ident
umapdata <- Embeddings(IntegratedSData_Processed, reduction = "umap")
ggplotDataframe2 <- data.frame(umapdata, umapdata2, groupdata2)
ggplotDataframe2$groupdata2 <- as.factor(ggplotDataframe2$groupdata2)
p1 <- ggplot(ggplotDataframe2, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = groupdata2), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe2, aes(x = umap_1.1, y = umap_2.1)) +
        geom_point(aes(color = groupdata2), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "NSCLC Integration after vs before ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Evaluate_Integration_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
