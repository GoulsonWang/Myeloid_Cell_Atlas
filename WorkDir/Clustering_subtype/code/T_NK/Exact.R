#2026/1/15
#本脚本用于提取T/NK细胞



#input：
#   1.Description1.RDS

#output：
#   1.T_NK_SData.RDS    T/NK提取后的细胞
#   2.SData_Processed.RDS   对每个样本的数据计算MT、HSP，计算SCT
#   3.deviation_contribution.pdf
#   4.SData_Integrated_rpca.RDS     rpca整合后的数据
#   5.SData_Integrated_Processed.RDS    整合后、再处理后的数据
#   6.Integration_RPCA_Result.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/T_NK")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/T_NK", sep = "")
data_path <- paste(file_path, "/data/T_NK", sep = "")
plot_path <- paste(file_path, "/plot/T_NK", sep = "")



#提交时使用
plan(multicore, workers = 18)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
SData_Overall <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/data/Description1.RDS")
countsData <- LayerData(SData_Overall, layer = "counts")
metaData <- SData_Overall@meta.data
New_SData <- CreateSeuratObject(countsData, meta.data = metaData)
str(New_SData)
SData <- subset(New_SData, subset = Lineage %in% c("T/NK"))
saveRDS(SData, file = file.path(data_path, "T_NK_SData.RDS"))
cat("已保存T_NK_SData.RDS")



#预处理
cat("检查SData Lineage")
table(SData$Lineage)
dim(SData)
cat("开始处理数据")
SData[["RNA"]] <- split(SData[["RNA"]], f = SData$orig.ident)
SData[["percent.HSP"]] <- PercentageFeatureSet(SData, pattern = "^HSP-")
SData_Processed <- SCTransform(SData, vars.to.regress = c("percent.MT", "percent.HSP"), verbose = F) %>%
    RunPCA(verbose = F)
p1 <- ElbowPlot(SData_Processed, ndims = 50)
ggsave(file.path(plot_path, "deviation_contribution.pdf"), plot = p1)
saveRDS(SData_Processed, file = file.path(data_path, "SData_Processed.RDS"))



#整合
cat("开始矫正批次效应")
start_time <- Sys.time()
SData_Integrated <- IntegrateLayers(
    SData_Processed, 
    method = RPCAIntegration, 
    normalization.method = "SCT",
    dims = 1:40 ,           #累计方差贡献率为87.5%
    scale.layer = "scale.data",
    k.weight = 39,          #当选择43（orig所含细胞数最小的值）时候报错信息为降低k.weight小于40， k.weight越小，矫正越严格。
    verbose = F             #注意：k.weight默认值为100  此处存在校正过度的风险
)
end_time <- Sys.time()
integration_time <- end_time - start_time
cat("批次效应矫正完成，耗时:", format(integration_time), "\n") 
saveRDS(SData_Integrated, file = file.path(data_path, "SData_Integrated_rpca.RDS"))



#precessing again
cat("开始聚类")
SData_Integrated_Processed <- FindNeighbors(
    SData_Integrated, 
    dims = 1:40, 
    reduction = "integrated.dr", 
    verbose = F
    ) %>%
    FindClusters(resolution = 0.3, verbose = F) %>%
    RunUMAP(dims = 1:40, reduction = "integrated.dr", verbose = F) %>%
    RunTSNE(dims = 1:40, reduction = "integrated.dr", verbose = F)
#注意：此处所有reduction参数均要改为integrated.dr，否则默认为pca
end_time2 <- Sys.time()
integration_time2 <- end_time2 - end_time
cat("批次效应矫正完成，耗时:", format(integration_time2), "\n") 
saveRDS(SData_Integrated_Processed, file = file.path(data_path, "SData_Integrated_Processed.RDS"))



#绘制整合后的umap图
umapdata <- Embeddings(SData_Integrated_Processed, reduction = "umap")
ClusterData <- SData_Integrated_Processed$SCT_snn_res.0.3
Orig <- SData_Integrated_Processed$orig.ident
PrimaryData <- SData_Integrated_Processed$PrimaryOrMet
table(SData_Integrated_Processed$PrimaryOrMet)      #转移灶：原发灶≈1：6， 与细胞量1：2不同
#          brain      chest wall           liver      lymph node peritoneal node 
#           1710            3036            4597            2338              61 
#        Primary 
#          66064 
ggplotDataframe <- data.frame(umapdata, ClusterData, Orig, PrimaryData)
library(RColorBrewer)
num_orig_levels <- length(unique(Orig))
color_palette <- if(num_orig_levels <= 8) brewer.pal(num_orig_levels, "Set2") else colorRampPalette(brewer.pal(9, "Set1"))(num_orig_levels)
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = ClusterData), size = 0.2, stroke = 0.2) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke  = 1.5))) 
p2 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
    geom_point(shape = 1, aes(color = Orig), size = 0.2, stroke = 0.2) +
    scale_color_manual(values = color_palette) + # 应用新的颜色调色板
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) 
p3 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = PrimaryData), size = 0.2, stroke = 0.2) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke  = 1.5))) 
p <- p1 + p2 + p3 +
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(
        title = "RPCA Integration Result", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integration_RPCA_Result.pdf", sep = ""),
    plot = p, 
    width = 23, 
    height = 7
)
