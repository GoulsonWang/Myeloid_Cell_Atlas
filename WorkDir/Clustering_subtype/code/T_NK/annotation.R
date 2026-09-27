#2026/1/15  本脚本用于注释T/NK细胞



#input：
#   1.SData_Integrated_Processed.RDS

#output：
#   1.NK_Marker_Distribution.pdf
#   2.SData_Annotated_By_Starcat_T.RDS



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
SData <- readRDS(file.path(data_path, "SData_Integrated_Processed.RDS"))    #6GB
colnames(SData@meta.data)
dim(SData)      #23609 109970
table(SData$Lineage)



#从中剔除NK细胞——marker:KLRF1,KLRD1,NCAM1
NK_Markers <- c("FGFBP2", "FCG3RA", "CX3CR1", "GNLY", "NKG7", "TYOBP", "PRF1")
p1 <- FeaturePlot(SData,
    features = NK_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "NK_Marker_Distribution.pdf"),
    plot = p1,
    width = 18, 
    height = 12
)
#判定cluster2为NK细胞，并删除
SData <- subset(SData, !SCT_snn_res.0.3 == 2)
table(SData$SCT_snn_res.0.3)
DimPlot(SData, reduction = "umap", label = TRUE)



#利用starcat注释
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(Matrix)
    library(Seurat)
    library(R.utils)
})
SData@active.assay <- "SCT"

counts = GetAssayData(SData, assay = "SCT", layer = "counts")
writeMM(counts, paste0(data_path, '/starcat_T/matrix.mtx'))
gzip(paste0(data_path, '/starcat_T/matrix.mtx'))

barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes), paste0(data_path, '/starcat_T/barcodes.tsv'),
           col_names = FALSE)
gzip(paste0(data_path, '/starcat_T/barcodes.tsv'))

gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
write_delim(as.data.frame(features), delim = "\t", paste0(data_path, '/starcat_T/features.tsv'),
           col_names = FALSE)
gzip(paste0(data_path, '/starcat_T/features.tsv'))



#提交starcat命令
output_name = 'starcat_result'
counts_fn = paste0(data_path, '/starcat_T/matrix.mtx.gz')

cmd = paste0('starcat', 
             ' --reference ', '"TCAT.V1"',
             ' --counts ', '"', counts_fn, '"', 
             ' --output-dir ', '"', data_path, '/"', 
             ' --name ', '"', output_name, '"' 
           )
cmd

system(cmd)



#载入结果
usage = read.table(paste0(data_path, "/", output_name, '.rf_usage_normalized.txt'))
scores = read.table(paste0(data_path, "/", output_name, '.scores.txt'))
usage %>% head(2)
scores %>% head()

cell_types <- scores[, "Multinomial_Label", drop = FALSE]  # 保持数据框结构
SData <- AddMetaData(SData,metadata = as.data.frame(cell_types))
SData@meta.data %>% colnames
table(SData$Multinomial_Label)
saveRDS(SData, file.path(data_path, "SData_Annotated_By_Starcat_T.RDS"))



#可视化
SData <- readRDS(file.path(data_path, "SData_Annotated_By_Starcat_T.RDS"))
p1 <- DimPlot(
    SData, 
    reduction = "umap", 
    label = TRUE, 
    group.by = "Multinomial_Label"
    )
ggsave(
    filename = file.path(plot_path, "T_Starcat_Annotation.pdf"),
    plot = p1,
    width = 7, 
    height = 7
)
