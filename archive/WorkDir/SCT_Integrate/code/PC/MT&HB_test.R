#2025/9/14      针对UMAP图表现不佳的8，9，10，13号样本，观察其碎细胞簇是否是MT&HB基因高表达


# Input: SDataList_Processed.RDS  (PC数据)
# Output: 
#   1.HB_test-1.pdf
#   2.MT_test-1.pdf
    

suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
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
SDataList <- readRDS(paste(data_path, "/SDataList_Processed.RDS", sep = ""))
SDataList_TargetSample <- list()
SDataList_TargetSample[[1]] <- SDataList[[8]]
SDataList_TargetSample[[2]] <- SDataList[[9]]
SDataList_TargetSample[[3]] <- SDataList[[10]]
SDataList_TargetSample[[4]] <- SDataList[[13]]


#检查HBgene表达情况
for(i in 1:length(SDataList_TargetSample)){
    #选择HBgene
    hb_genes <- rownames(SDataList_TargetSample[[i]]@assays$SCT@scale.data) # 假设您想在第一个样本上查看
    hb_genes <- hb_genes[grep("^HB", hb_genes)]
    #计算平均表达量
    scale_data <- GetAssayData(SDataList_TargetSample[[i]], assay = "SCT", layer = "scale.data")
    hb_avg_expression <- colMeans(scale_data[hb_genes, ]) # 注意：矩阵是 genes x cells
    SDataList_TargetSample[[i]]@meta.data$HB_Avg_Expression <- hb_avg_expression
    p_avg_feature <- FeaturePlot(SDataList_TargetSample[[i]], features = "HB_Avg_Expression", 
                               reduction = "umap", 
                               # 可选参数
                               cols = c("lightgrey", "blue"), # 从低到高: 灰色 -> 蓝色
                               pt.size = 0.8,
                               order = TRUE) # 高表达的点在上层
    ggsave(paste(plot_path, "/HB_test-", i, ".pdf", sep = ""), plot = p_avg_feature)
}


#检查MTgene表达情况
for(i in 1:length(SDataList_TargetSample)){
    #选择MTgene
    MT_genes <- rownames(SDataList_TargetSample[[i]]@assays$SCT@scale.data) # 假设您想在第一个样本上查看
    MT_genes <- MT_genes[grep("^MT", MT_genes)]
    #计算平均表达量
    scale_data <- GetAssayData(SDataList_TargetSample[[i]], assay = "SCT", layer = "scale.data")
    MT_avg_expression <- colMeans(scale_data[MT_genes, ]) # 注意：矩阵是 genes x cells
    SDataList_TargetSample[[i]]@meta.data$MT_Avg_Expression <- MT_avg_expression
    p_avg_feature <- FeaturePlot(SDataList_TargetSample[[i]], features = "MT_Avg_Expression", 
                               reduction = "umap", 
                               # 可选参数
                               cols = c("lightgrey", "blue"), # 从低到高: 灰色 -> 蓝色
                               pt.size = 0.8,
                               order = TRUE) # 高表达的点在上层
    ggsave(paste(plot_path, "/MT_test-", i, ".pdf", sep = ""), plot = p_avg_feature)
}


