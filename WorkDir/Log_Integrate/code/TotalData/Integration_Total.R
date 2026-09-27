#2025/9/29      将所有数据整合起来（rpca）



# Input: 
#   1.Merge_Data.RDS        Merge合并后的Seurat对象

# Output: 
#   1.Total_SData.RDS     整合四个数据集后的Seurat对象。
#   2.SData_Integrated_Processed.RDS       #合并后对所有数据再次进行处理
#   3.Integrate_Data_CancerType_Plot.pdf       #整合后对所有数据处理得到UMAP和TSNE图
#   4.Evaluate_Integration_Plot.pdf         #整合前后对比UMAP图结果


#提交命令时候需要修改的地方：
#1.plan()
#2.SData <- readRDS()
#3.options()





suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")

#载入数据
SData_Merged <- readRDS(file.path(data_path, "Merge_Data_Supplement.RDS"))
#注意：PC数据集中orig.ident=4的样本细胞只有177个，故筛除。剩余68个样本
SData_Merged <- subset(SData_Merged, orig.ident != 4) 



#preparing
SData_Merged <- NormalizeData(SData_Merged, verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
    RunPCA(verbose = F)

#str(CancerData)



#Integration
SData_Integrated <- IntegrateLayers(SData_Merged, 
    method = RPCAIntegration, 
    normalization.method = "LogNormalize",
    scale.layer = "scale.data",
    verbose = F)

DefaultAssay(SData_Integrated) <- "RNA"
SData_Integrated <- JoinLayers(SData_Integrated, assay = "RNA")
cat("\n", "输出str(SData_Integrated)")
str(SData_Integrated)        #结果保存为str(Total_SData).txt
saveRDS(SData_Integrated, paste(data_path, "SData_Integrated.RDS", sep = "/"))



#precessing again
SData_Integrated_Processed <- FindNeighbors(SData_Integrated, reduction = "integrated.dr", verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F) %>%
    RunUMAP(dims = 1:30, reduction = "integrated.dr", verbose = F) %>%
    RunTSNE(dims = 1:30, reduction = "integrated.dr", verbose = F)
#注意：此处所有reduction参数均要改为integrated.dr，否则默认为pca
cat("\n", "输出str(SData_Integrated_Processed)")
str(SData_Integrated_Processed)     #结果保存为str(SData_Integrated_Processed).txt
saveRDS(SData_Integrated_Processed, file = paste(data_path, "SData_Integrated_Processed", sep = "/"))



#画图 umap+tsne(整合后)
umapdata <- Embeddings(SData_Integrated_Processed, reduction = "umap")
tsnedata <- Embeddings(SData_Integrated_Processed, reduction = "tsne")
ClusterData <- SData_Integrated_Processed@active.ident
ggplotDataframe <- data.frame(umapdata, tsnedata, ClusterData)

p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = ClusterData), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe, aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(color = ClusterData), size = 0.02) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "Total Integrated Data ClusterData Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_ClusterData_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)



#整合前后对比
UnintegratedSData <- LayerData(SData_Integrated_Processed, assay = "RNA", layer = "counts")
metadata <- SData_Integrated_Processed@meta.data
UnintegratedSData <- CreateSeuratObject(
    counts = UnintegratedSData, 
    meta.data = metadata
    )
UnintegratedSData <- NormalizeData(UnintegratedSData, verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(vars.to.regress = "percent.MT", verbose = F) %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)
umapdata2 <- Embeddings(UnintegratedSData, reduction = "umap")
groupdata2 <- UnintegratedSData@meta.data$CancerType
ggplotDataframe2 <- data.frame(umapdata, umapdata2, groupdata2)
p1 <- ggplot(ggplotDataframe2, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = groupdata2), size = 0.02) +
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
        title = "Total Integration after vs before ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Evaluate_Integration_Plot.pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)
