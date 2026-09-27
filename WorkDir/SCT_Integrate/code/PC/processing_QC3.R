# 2025/9/14      processing

# Input: QC3SDataList.RDS   (PC数据)
# Output:
#   1.SDataList_Processed(QC3).RDS  -----把QC3后的数据进行预处理
#   2.QC3DimPlot_Origident-i.pdf  -----预处理结果画图UMAP


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


# 载入数据
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/PC/QC3SDataList.RDS")
#输入形式为RDS，文件格式为list，其中每个子文件为一个Seurat对象
# str(SDataList)     结果保存至/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/data/PC/str(QC2SeuratDataList).txt



# SCTransform
SDataListTransformed <- list()
nlist <- length(SDataList)

options(future.globals.maxSize = 10 * 1024^3) # 提高内存保护机制阈值为10G
for (i in 1:nlist) {
    SDataListTransformed[[i]] <- SCTransform(SDataList[[i]], vars.to.regress = "percent.MT", verbose = FALSE)
}




# 聚类画图
ggplotDataList <- list()
for (i in 1:nlist) {
    SDataListTransformed[[i]] <- RunPCA(SDataListTransformed[[i]], verbose = F) %>%
        RunUMAP(dims = 1:30, verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(verbose = F)
    # str(SDataListTransformed[[1]])   #结果保存为str(QC3SDataListTransformed).txt
    # 准备画图的dataframe
    embeddings <- Seurat::Embeddings(SDataListTransformed[[i]], reduction = "umap") # 返回一个矩阵，行是细胞，列是维度
    group_info <- SDataListTransformed[[i]]@active.ident
    names(group_info) <- colnames(SDataListTransformed[[i]]) # 确保 names 对应细胞名
    ggplotDataList[[i]] <- data.frame(embeddings, stringsAsFactors = FALSE)
    ggplotDataList[[i]]$ggplotgroup <- group_info[rownames(ggplotDataList[[i]])] # 按细胞名匹配
    # 画图
    p <- ggplot(ggplotDataList[[i]], aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = ggplotgroup), size = 0.1) +
        theme_light() +
        labs(title = paste0("QC3DimPlot", "origident-", i)) +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
    ggsave(paste(plot_path, "/QC3DimPlot_Origident-", i, ".pdf", sep = ""), plot = p)
}


# 看图，8，9，10，13四个样本效果较差，具体表现为UMAP图上细碎小细胞簇太多。
# 猜测可能原因：QC时候percent.mt和percent.HB阈值过于宽松
# 处理方式：进行QC3————收紧percent.MT&HB阈值为15和3


saveRDS(SDataListTransformed, file = paste(data_path, "/SDataList_Processed(QC3).RDS", sep = ""))
