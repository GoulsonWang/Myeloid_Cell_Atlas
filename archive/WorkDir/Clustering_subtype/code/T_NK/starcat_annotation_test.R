#2026/1/15
#本脚本用于测试starcat用作注释的效果



suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(Matrix)
    library(Seurat)
    library(R.utils)
})
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/T_NK")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/T_NK", sep = "")
data_path <- paste(file_path, "/data/T_NK", sep = "")
plot_path <- paste(file_path, "/plot/T_NK", sep = "")



#下载示例数据
#library(curl)
#curl_download('https://zenodo.org/records/13368041/files/COMBAT-CITESeq-DATA.Raw.T.ADTfixed20230831FiltForcNMF.Downsampled.rds?download=1',
#              file.path(data_path, 'example_data.rds'), 
#              handle = new_handle(timeout = 1000))
seu_object = readRDS(file.path(data_path, 'example_data.rds'))
seu_object
#An object of class Seurat 
#20957 features across 13800 samples within 1 assay 
#Active assay: RNA (20957 features, 0 variable features)
# 2 layers present: counts, data
seu_object@meta.data %>% colnames
#[1] "orig.ident"   "nCount_RNA"   "nFeature_RNA"
counts = seu_object@assays$RNA@counts
table(seu_object$orig.ident)
str(seu_object)

# Output counts matrix
writeMM(counts, paste0(data_path, '/starcat_test/matrix.mtx'))
gzip(paste0(data_path, '/starcat_test/matrix.mtx'))

# Output cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes), paste0(data_path, '/starcat_test/barcodes.tsv'),
           col_names = FALSE)
gzip(paste0(data_path, '/starcat_test/barcodes.tsv'))

# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
write_delim(as.data.frame(features), delim = "\t", paste0(data_path, '/starcat_test/features.tsv'),
           col_names = FALSE)
gzip(paste0(data_path, '/starcat_test/features.tsv'))



#提交starcat命令
#示例：starcat --reference "TCAT.V1" --counts "counts_fn" --output-dir "output_dir" --name "outuput_name
output_name = 'example_data'
counts_fn = paste0(data_path, '/starcat_test/matrix.mtx.gz')
     

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

seu_object <- AddMetaData(seu_object,metadata = as.data.frame(cell_types))
seu_object@meta.data %>% colnames() 
