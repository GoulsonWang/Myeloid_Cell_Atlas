# 2025/11/25     本脚本用与对R组样本进行NMF



# Input:
#   1.SData_Cluster030.RDS      FindMarker6_Supplement4.R的output

# Output:
#   1.R_Myeloid_List(For_NMF).RDS       以样本为单位保存为list，仅仅包含Macro/Mono的细胞
#   2.Macro_NMF_Result/Genes_nmf_w_basis_[sample_id].RData      NMF结果



suppressMessages(library(ggplot2))
suppressMessages(library(NMF))
suppressMessages(library(Seurat))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/NMF/code")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/NMF"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#以下注释掉的代码在本地运行，未提交至作业



## 载入数据
#SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/SData_Cluster030.RDS") # 6GB
## 修改其中的错误：orig = 79的样本，改Response信息"NA"为"No"
#idx_to_modify <- which(SData@meta.data$orig.ident == "79" & is.na(SData@meta.data$Response.Comprehensive))
#cat("找到", length(idx_to_modify), "个orig.ident为79且Response.Comprehensive为NA的样本\n")
#SData@meta.data[idx_to_modify, "Response.Comprehensive"] <- "No"
#print(table(SData@meta.data[SData@meta.data$orig.ident == "79", "Response.Comprehensive"], useNA = "always"))
## 添加Major_Subtype信息
#major_annotations <- rep("Macro/Mono", ncol(SData))
#current_clusters <- Idents(SData)
#major_annotations[which(current_clusters == "4")] <- "Mast"
#major_annotations[which(current_clusters %in% c("5", "15", "17"))] <- "DCs"
#names(major_annotations) <- colnames(SData)
#SData <- AddMetaData(SData, metadata = major_annotations, col.name = "Major_Cell_Type")
#table(SData$Major_Cell_Type)
##       DCs     Macro/Mono       Mast 
##        7031      64756         6019 
#
#
#
##筛选R组样本
#SData_R <- subset(SData, Response.Comprehensive == "Yes" & Major_Cell_Type == "Macro/Mono")
##dim(SData_R)    #25279 27147
#table(SData_R$orig.ident)
#
## 将SData_R按样本分成若干Seurat对象，并集成为一个list
#cat("\n将SData_R按样本分割成多个Seurat对象...\n")
#orig_idents <- unique(SData_R$orig.ident)
#seurat_list <- list()
#
#for(orig in orig_idents) {
#    seurat_list[[as.character(orig)]] <- subset(SData_R, orig.ident == orig)
#    cat("样本", orig, "已添加到列表，基因数:", nrow(seurat_list[[as.character(orig)]]), 
#        "，细胞数:", ncol(seurat_list[[as.character(orig)]]), "\n")
#}
#
## 对每个样本处理：删除在所有细胞中表达量都为0的基因
#cat("\n开始处理每个样本，删除在所有细胞中表达量都为0的基因...\n")
#for(i in 1:length(seurat_list)) {
#    sample_name <- names(seurat_list)[i]
#    cat("处理样本:", sample_name, "\n")
#    
#    # 查看处理前的维度
#    dim_before <- dim(seurat_list[[i]])
#    cat("  处理前 - 基因数:", dim_before[1], ", 细胞数:", dim_before[2], "\n")
#    
#    # 删除在所有细胞中表达量都为0的基因
#    # 获取表达矩阵
#    expr_matrix <- GetAssayData(seurat_list[[i]], slot = "counts")
#    
#    # 找到至少在一个细胞中表达的基因（表达量>0）
#    expressed_genes <- rowSums(expr_matrix) > 0
#    
#    # 子集化只保留在至少一个细胞中表达的基因
#    seurat_list[[i]] <- subset(seurat_list[[i]], features = names(expressed_genes)[expressed_genes])
#    
#    # 查看处理后的维度
#    dim_after <- dim(seurat_list[[i]])
#    cat("  处理后 - 基因数:", dim_after[1], ", 细胞数:", dim_after[2], "\n")
#    cat("  共删除", (dim_before[1] - dim_after[1]), "个在所有细胞中都不表达的基因\n")
#}
#
## 输出处理后各样本的信息
#cat("\n处理完成后各样本信息:\n")
#for(i in 1:length(seurat_list)) {
#    name <- names(seurat_list)[i]
#    obj <- seurat_list[[i]]
#    cat("样本", name, "- 基因数:", nrow(obj), ", 细胞数:", ncol(obj), "\n")
#}
#saveRDS(seurat_list, file = file.path(data_path, "R_Myeloid_List(For_NMF).RDS"))



#以下代码提交至作业
seurat_list <- readRDS(file.path(data_path, "R_Myeloid_List(For_NMF).RDS"))
# 定义NMF程序函数
nmf_programs <- function(cpm, is.log = F, rank, method = "snmf/r", seed = 1) {
    if (is.log == F) CP100K_log <- log2((cpm / 10) + 1) else CP100K_log <- cpm
    CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > 3.5)) > ncol(CP100K_log) * 0.02), ]
    CP100K_log <- CP100K_log - rowMeans(CP100K_log)
    CP100K_log[CP100K_log < 0] <- 0
    nmf_programs <- nmf(as.matrix(CP100K_log),
        rank = rank, method = method, seed = seed, nrun = 10,
        .opt = list(parallel = TRUE, verbose = TRUE)
    )
    return(nmf_programs)
}

# 对每个样本进行NMF分析
for (i in 1:length(seurat_list)) {
    sample_id <- names(seurat_list)[i]
    print(paste0("############# Start ", sample_id, " ", Sys.time()))
    Genes_nmf_w_basis <- list()
    # 获取表达数据
    expr <- seurat_list[[i]]@assays$RNA$counts
    # 标准化到CPM
    expr <- t(t(expr) / colSums(expr)) * 1000000
    # 设置并行核心数
    NMF::nmf.options(cores = 30)
    # 运行NMF
    NMFs <- nmf_programs(expr, is.log = T, rank = 4:9, method = "snmf/r", seed = 1)
    # 提取不同rank下的basis矩阵
    target_ranks <- 5:9
    temp_rank4_9_nruns10 <- NMF::basis(NMFs$fit[[as.character(4)]])
    colnames(temp_rank4_9_nruns10) <- paste0(sample_id, "_rank4_9_nruns10", ".RDS.4.", 1:4)
    for (r in target_ranks) {
        w <- NMF::basis(NMFs$fit[[as.character(r)]])
        colnames(w) <- paste0(sample_id, "_rank4_9_nruns10", ".RDS.", r, ".", 1:r)
        temp_rank4_9_nruns10 <- cbind(temp_rank4_9_nruns10, w)
    }
    Genes_nmf_w_basis[[sample_id]] <- temp_rank4_9_nruns10
    print(paste0("############# Stop ", sample_id, Sys.time()))

    # 保存结果
    file_path <- paste0(data_path, "/R_Macro_NMF_Result/Genes_nmf_w_basis_", sample_id, ".RData")
    save(Genes_nmf_w_basis, file = file_path)
}