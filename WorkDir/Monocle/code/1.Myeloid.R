#2026/1/3
#本脚本用于对单核细胞-巨噬细胞拟时序分析



#input
#   1.SData_Minor_Cell_Type.RDS

#output
#   1.umap_preparing_By_monocle.pdf     monocle预处理后得到的umap图
#   2.umap_preparing_By_Seurat.pdf      将Seurat得到的umap结果导入到cds文件中
#   3.Cluster_Partition_By_monocle.pdf"     利用monocle对细胞进行分类，对各类分别进行拟时序分析
#   4.trajectory1.pdf   拟时序分析轨迹图
#   5. trajectory_setPresudoTime0_1.pdf     单起点拟时序分析结果
#   6.trajectory_setPresudoTime0_2.pdf      双起点结果
#   7.gene_fits.RDS     用pseudotime来作为广义线性回归的解释变量，被解释变量为基因的表达量
#   8.



suppressMessages(library(monocle3))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(ggplot2))



# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/code/")

#提交时使用
plan(multicore, workers = 6)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

# 载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS")    #6GB



#从数据中剔除DC、Mast、Neutro、以及组织驻留巨噬细胞(Macro_5)
table(SData$Minor_Cell_Type)
cells_to_remove <- colnames(SData)[SData@meta.data$Minor_Cell_Type %in% c("DC_1", "DC_2",  "DC_3", "Mast", "Neutro_1", "Macro_5")]
cat("将删除", length(cells_to_remove), "个细胞")

# 然后从SData中剔除这些细胞
SData_clean <- SData[, !colnames(SData) %in% cells_to_remove]



#构建monocle3对象cds
SData_clean
data <- LayerData(SData_clean, assay = "RNA", layer =  "counts")
metadata <-  SData_clean@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data, cell_metadata = metadata, gene_metadata = gene_annotation)

#预处理cds
cat(" 开始预处理...\n")
startTime <- Sys.time()
cds <- preprocess_cds(cds, method = "PCA") 
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", cores = 2)
endTime <- Sys.time()
cat("预处理耗时:", endTime - startTime, "\n")
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "Minor_Cell_Type", show_trajectory_graph = T)
ggsave(file.path(plot_path, "umap_preparing_By_monocle.pdf"), width = 10, height = 10)
cat("已输出Monocle预处理得到的umap图")

#将Seurat的UMAP图导入
Seurat_umap <- Embeddings(SData_clean, reduction = "umap")
Seurat_umap <- Seurat_umap[rownames(cds@int_colData$reducedDims$UMAP), ]    #排序
cds@int_colData$reducedDims$UMAP <- Seurat_umap
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "Minor_Cell_Type", show_trajectory_graph = T)
ggsave(file.path(plot_path, "umap_preparing_By_Seurat.pdf"), width = 10, height = 10)
cat("已输出Seurat得到的umap坐标导入到cds后的结果图，请用做对比检查导入是否正确")
#检查正确



#构建细胞轨迹
cat("开始构建细胞轨迹...\n")
cds <- cluster_cells(cds)   #聚类，结果只有1类，符合预期
p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = F)
ggsave(file.path(plot_path, "Cluster_Partition_By_monocle.pdf"), width = 10, height = 10)
cds <- learn_graph(
    cds, 
    close_loop = F,     #结果图中尽量避免闭环存在
    learn_graph_control = list(minimal_branch_len = 12)     #小于12的分支被修剪
    )
p <- plot_cells(
    cds, 
    color_cells_by = "partition", 
    label_principal_points = TRUE, 
    label_roots = T, 
    label_leaves = T,
    label_branch_points = T,
    graph_label_size = 1
    )
ggsave(file.path(plot_path, "trajectory1.pdf"), plot = p, width = 6, height = 5)
cat( "已构建细胞轨迹")

#拟时序排列细胞
cds <- order_cells(cds, root_pr_nodes = c("Y_312"))        #单起点图
p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5
           )
ggsave(file.path(plot_path, "trajectory_setPresudoTime0_2.pdf"), plot = p, width = 6.5, height = 5)

cds <- order_cells(cds, root_pr_nodes = c("Y_312", "Y_149"))        #双起点图
p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5
           )
ggsave(file.path(plot_path, "trajectory_setPresudoTime0_2.pdf"), plot = p, width = 5, height = 6)
#比较后发现双起点更符合预期结果



#查看随时间发生显著变化的基因
Variable_GeneSet <- read.csv("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SCT_HVG_3000_genes.csv")
Variable_GeneSet_vector <- as.vector(unlist(Variable_GeneSet))
cds_subset <- cds[rowData(cds)$gene_short_name %in% Variable_GeneSet_vector,]
gene_fits <- fit_models(
    cds_subset, 
    model_formula_str = "~pseudotime", 
    cores = 5
    )
saveRDS(gene_fits, file = file.path(data_path, "gene_fits.RDS"))
#3000个基因中2200+个基因都显著。

save_monocle_objects(
    cds_subset, 
    directory_path = file.path(data_path, "CDS_VariableGene_Only"), 
    comment='CDS_Only_VariableGene_3000'
    )
save_monocle_objects(
    cds, 
    directory_path = file.path(data_path, "CDS_AllGene"), 
    comment='CDS_AllGene'
    )