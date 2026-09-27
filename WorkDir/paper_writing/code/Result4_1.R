#本脚本用于画不同通路随轨迹的变化曲线：找出目标通路及其相关基因



#input
#   1.FCList_Branch.RDS

#output
#   1.target_pathway_geneset.RDS



suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
FC_result <- readRDS(file.path("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/data/FCList_Branch.RDS"))

#分别对四条branch进行富集分析
KEGG_branch_gesa_result <- list()
for(i in names(FC_result)){
    GeneID_kegg <- bitr(rownames(FC_result[[i]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    geneList <- setNames(FC_result[[i]][GeneID_kegg$SYMBOL, "avg_log2FC"], GeneID_kegg$ENTREZID)
    geneList <- geneList[!is.na(names(geneList))]
    KEGG_branch_gesa_result[[i]] <- gseKEGG(
        gene = geneList,
        organism = "hsa",
        keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,          
        minGSSize = 10,
        maxGSSize = 500,
        verbose = F
    )
    if(nrow(KEGG_branch_gesa_result[[i]]) > 0){
        p1 <- dotplot(KEGG_branch_gesa_result[[i]], showCategory=30) + ggplot2::ggtitle(paste0("dotplot for GSEA", ":", i))
        ggplot2::ggsave(plot = p1, filename = paste(plot_path, "/", "KEGG_GSEA", i, ".pdf", sep = ""), width = 10, height = 10)
        cat("已保存", i, "的GSEA结果图\n")
    } else {
        cat("注意：", i, "未富集出任何通路\n")
    }
}


KEGG_branch_gesa_result[["MS4A6A_1"]]$Description
KEGG_branch_gesa_result[["MS4A6A_2"]]$Description
KEGG_branch_gesa_result[["APOE_1"]]$Description
KEGG_branch_gesa_result[["APOE_2"]]$Description
#提取补充通路的基因名称
pathway_name <- c(
    "MAPK signaling pathway", "Complement and coagulation cascades", 
    "Antigen processing and presentation", "Metabolic pathways", 
    "Carbon metabolism", "Cholesterol metabolism"
)

# 提取指定通路的core_enrichment并将其转换为SYMBOL格式
extract_pathway_genes <- function(gsea_results, pathway_names) {
    extracted_genes <- list()
    for(branch_name in names(gsea_results)) {
        branch_result <- gsea_results[[branch_name]]
        if(nrow(branch_result) == 0) next
        target_pathways <- list()
        for(pathway in pathway_names) {
            matching_rows <- which(grepl(pathway, branch_result$Description, fixed = TRUE))
            if(length(matching_rows) > 0) {
                best_match <- matching_rows[which.min(branch_result$p.adjust[matching_rows])]
                core_enrichment <- branch_result$core_enrichment[best_match]
                entrez_ids <- strsplit(core_enrichment, "/")[[1]]
                gene_symbols <- bitr(entrez_ids, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
                target_pathways[[pathway]] <- gene_symbols$SYMBOL
                cat("在", branch_name, "中找到通路", pathway, "，包含", length(gene_symbols$SYMBOL), "个基因\n")
            } else {
                cat("在", branch_name, "中未找到通路", pathway, "\n")
            }
        }
        
        if(length(target_pathways) > 0) {
            extracted_genes[[branch_name]] <- target_pathways
        }
    }
    return(extracted_genes)
}
target_pathway_genes <- extract_pathway_genes(KEGG_branch_gesa_result, pathway_name)
#MAPK signaling pathway ，只找到 14 个基因
#Complement and coagulation cascades ，找到 9 个基因
#Antigen processing and presentation ，找到 18 个基因
#Metabolic pathways ，包含 22 个基因
#Metabolic pathways ，包含 47 个基因
#Carbon metabolism ，包含 9 个基因
#Cholesterol metabolism ，包含 9 个基因



# 汇总不同branch中相同通路的基因，生成data.frame
# 行名为通路名称，行元素为基因名称，不同branch中相同通路合并时，基因名称取并集
summarize_pathway_genes <- function(pathway_genes_list) {
    all_pathways <- unique(unlist(lapply(pathway_genes_list, names)))
    combined_pathway_genes <- list()
    for(pathway in all_pathways) {
        genes_for_pathway <- c()
        for(branch_name in names(pathway_genes_list)) {
            if(pathway %in% names(pathway_genes_list[[branch_name]])) {
                genes_for_pathway <- union(genes_for_pathway, pathway_genes_list[[branch_name]][[pathway]])
            }
        }
        combined_pathway_genes[[pathway]] <- genes_for_pathway
    }
    max_genes <- max(sapply(combined_pathway_genes, length))
    pathway_df <- data.frame(matrix(ncol = max_genes, nrow = length(all_pathways)))
    rownames(pathway_df) <- all_pathways
    for(pathway in all_pathways) {
        genes <- combined_pathway_genes[[pathway]]
        pathway_df[pathway, 1:length(genes)] <- genes
    }
    colnames(pathway_df) <- paste0("Gene_", 1:max_genes)
    return(list(
        pathway_gene_matrix = pathway_df,
        pathway_gene_lists = combined_pathway_genes
    ))
}
summary_result <- summarize_pathway_genes(target_pathway_genes)
# 输出汇总信息
cat("\n=== 汇总结果信息 ===\n")
cat("总共找到", length(summary_result$pathway_gene_lists), "个通路\n")
for(pathway in names(summary_result$pathway_gene_lists)) {
    gene_count <- length(summary_result$pathway_gene_lists[[pathway]])
    cat("'", pathway, "' 包含 ", gene_count, " 个唯一基因\n", sep="")
}
saveRDS(summary_result, file = file.path(data_path, "target_pathway_geneset.RDS"))
