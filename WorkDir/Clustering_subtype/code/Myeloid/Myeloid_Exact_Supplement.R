#2025/11/1      本脚本用于对聚类的resolution参数改为0.3， Myeloid_Exact.R的UMAP图进行美化


#input：
#   1.SData_Merged_Processed.RDS

#output：
#   1.SData_Merged_Processed.RDS        将分类的resolution参数改为0.3（原来是1）
#   2.Batch_Identification(Supplement).pdf



suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "SData_Merged_Processed.RDS"))    #6GB
#str(SData)
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值
#SData <- FindClusters(SData, resolution = 0.3, verbose = F)        #改resolution = 0.3      聚类个数是12类
#saveRDS(SData, file = file.path(data_path, "SData_Merged_Processed.RDS"))


#绘图
umapdata <- Embeddings(SData, reduction = "umap")
ClusterData <- SData$SCT_snn_res.0.3
CancerType <- SData$CancerType
Orig <- SData$orig.ident
Patient <- SData$Patient.ID
ggplotDataframe <- data.frame(umapdata, ClusterData, CancerType, Orig, Patient)

# 为'Orig'创建一个色彩丰富的调色板
num_orig_levels <- length(unique(Orig))
color_palette <- if(num_orig_levels <= 8) brewer.pal(num_orig_levels, "Set2") else colorRampPalette(brewer.pal(9, "Set1"))(num_orig_levels)

#计算每个 orig.ident 的 UMAP 中心坐标
label_data <- ggplotDataframe %>%
  group_by(Orig) %>%
  summarise(
    umap_1 = mean(umap_1, na.rm = TRUE),
    umap_2 = mean(umap_2, na.rm = TRUE),
    .groups = 'drop'
  )
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = ClusterData), size = 0.1, stroke = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) 
p2 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = CancerType), size = 0.1, stroke = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) 
p3 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
    geom_point(shape = 1, aes(color = Orig), size = 0.1, stroke = 0.1) +
    scale_color_manual(values = color_palette) + # 应用新的颜色调色板
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) 
p4 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = Patient), size = 0.1, stroke = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) 
p <- p1 + p2 + p3 + p4 + 
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "Identify the origin of batch effect(Supplement) ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Batch_Identification(Supplement).pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 14
)
p_Orig_CancerType <- p2 + p3 + 
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "Identify the origin of batch effect(Supplement) ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Batch_Identification(Supplement)_Orig_CancerType.pdf", sep = ""),
    plot = p_Orig_CancerType, 
    width = 16, 
    height = 7
)


#思考：批次效应应该是以orig.ident为单位而存在的。
#但此时orig.ident之间的细胞量差距已经非常大，45~5155。Seurat-RPCA方法失败
#Patient.ID之间的细胞量差距也是：43~5155