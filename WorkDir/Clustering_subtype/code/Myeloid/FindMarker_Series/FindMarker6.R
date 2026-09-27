##2025/11/18    本脚本用于重新SCTransform，主要是为了排除掉HSP系列基因对基因方差的影响
#代码主要对照Integration1.R及FindMarker2.R



#input：
#   1.Myeloid_SData.RDS

#output：
#   1.SData_Processed.RDS               SCTransform后的数据
#   2.deviation_contribution.pdf        主成分标准差贡献图
#   3.SData_Integrated_rpca.RDS
#   4.SData_Integrated_Processed.RDS
#   5.Integration_RPCA_Result.pdf       整合后的umap图
#   6.Mast_Markers_Distribution.pdf
#   7.DCs_Markers_Distribution.pdf
#   8.Macro_Mono_Markers_Distribution.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid/FindMarker_Series")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid/FindMarker_Series", sep = "")
data_path <- paste(file_path, "/data/Myeloid/HSP_Exclude", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid/HSP_Exclude", sep = "")



#提交时使用
plan(multicore, workers = 30)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 10 * 1024^3)   



#载入数据
#SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/Myeloid_SData.RDS")
#SData[["RNA"]] <- split(SData[["RNA"]], f = SData$orig.ident)



#预处理
#SData[["percent.HSP"]] <- PercentageFeatureSet(SData, pattern = "^HSP-")
#SData_Processed <- SCTransform(SData, vars.to.regress = c("percent.MT", "percent.HSP"), verbose = F) %>%
#    RunPCA(verbose = F)
#p1 <- ElbowPlot(SData_Processed, ndims = 50)
#ggsave(file.path(plot_path, "deviation_contribution.pdf"), plot = p1)
#saveRDS(SData_Processed, file = file.path(data_path, "SData_Processed.RDS"))



#整合
SData_Processed <- readRDS(file.path(data_path, "SData_Processed.RDS"))
#计算累计方差贡献率
#pct <- SData_Integrated [["pca"]]@stdev / sum( SData_Integrated [["pca"]]@stdev) * 100
#cumu <- cumsum(pct)
#cumu
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



#利用Marker区分Major_Subtype
Mast_Markers <- c("TPSAB1","CPA3","KIT","TPSB2")
DCs_Markers <- c("CD1C", "FLT3")
Macro_Mono_Markers <- c("CD68", "CD163", "APOE", "CD14")

#绘制各Marker的分布图
p1 <- FeaturePlot(SData_Integrated_Processed,
    features = Mast_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 4,
    min.cutoff = "q10",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Mast_Markers_Distribution.pdf"),
    plot = p1,
    width = 12, 
    height = 12
)

p2 <- FeaturePlot(SData_Integrated_Processed,
    features = DCs_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 2,     #注意：此处cutoff值为绝对值啊
    raster = F
)
ggsave(
    filename = file.path(plot_path, "DCs_Markers_Distribution.pdf"),
    plot = p2,
    width = 12, 
    height = 6
)

p3 <- FeaturePlot(SData_Integrated_Processed,
    features = Macro_Mono_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 4,
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Macro_Mono_Markers_Distribution.pdf"),
    plot = p3,
    width = 12, 
    height = 12
)