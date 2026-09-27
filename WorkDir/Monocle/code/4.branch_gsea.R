#2026/1/12  对五条branch的结果进行富集分析，从中挑出比较重要的通路及其相关基因集



#input
#   1.FCList_Branch.RDS

#output
#   1.KEGG_branch_gesa_result.RDS   #富集结果
#   2.target_pathway_geneset.RDS
#   3.KEGG_GSEAbranch1_1.pdf
#   4.KEGG_GSEAbranch2_1.pdf
#   5.KEGG_GSEAbranch1_2.pdf
#   6.KEGG_GSEAbranch2_2.pdf



suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))



# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/code/")

#提交时使用
#plan(multicore, workers = 6)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值



#载入数据
FC_result <- readRDS(file.path(data_path, "FCList_Branch.RDS"))     #已经经过排序了



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

#需要关注的通路：
#   最重要：
#       Lysosome， TNF signaling pathway， NOD-like receptor signaling pathway， PPAR signaling pathway， IL-17 signaling pathway
#   次重要：
#       NF-kappa B signaling pathway， RIG-I-like receptor signaling pathway， Lipid and atherosclerosis， Cell adhesion molecule (CAM) interaction



#提取关键通路的基因名称
pathway_name <- c("Lysosome", "TNF signaling pathway", "NOD-like receptor signaling pathway", 
    "PPAR signaling pathway", "IL-17 signaling pathway", "NF-kappa B signaling pathway", 
    "RIG-I-like receptor signaling pathway", "Lipid and atherosclerosis", "Cell adhesion molecule (CAM) interaction")
#important_macrophage_pathways_in_tme <- c(
#  "TNF signaling pathway",                    # TNF-α信号，调节炎症、细胞存活/凋亡、血管生成
#  "PPAR signaling pathway",                   # PPAR信号，调节代谢和M2极化，促免疫抑制
#  "NOD-like receptor signaling pathway",      # NOD样受体信号，识别DAMPs，激活炎症小体
#  "IL-17 signaling pathway",                  # IL-17信号，响应IL-17，促中性粒细胞招募/血管化
#  "MAPK signaling pathway",                   # MAPK信号，调控增殖、分化、存活、激活
#  "NF-kappa B signaling pathway",             # NF-κB信号，核心转录因子，调控炎症/免疫
#  "HIF-1 signaling pathway",                  # HIF-1信号，感应缺氧，促M2极化/血管生成
#  "Lysosome",                                 # 溶酶体，主要吞噬消化功能所需
#  "Cell cycle",                               # 细胞周期，指示增殖能力
#  "Cell adhesion molecule (CAM) interaction", # 细胞粘附分子互作，调节粘附/迁移
#  "Antigen processing and presentation",      # 抗原加工呈递，核心APC功能
#  "Metabolic pathways",                       # 代谢通路，与极化状态紧密相关
#  "Carbon metabolism",                        # 碳代谢，提供能量和生物合成原料
#  "Cholesterol metabolism",                   # 胆固醇代谢，影响膜结构和信号
#  "Glycerolipid metabolism",                  # 甘油酯代谢，与脂质存储/信号相关
#  "Glutathione metabolism",                   # 谷胱甘肽代谢，抗氧化防御
#  "Lipid and atherosclerosis",                # 脂质与动脉粥样硬化，与脂质代谢/PPAR相关
#  "Cytosolic DNA-sensing pathway",            # 胞质DNA感知通路，激活IFN等反应
#  "Necroptosis"                               # 坏死性凋亡，细胞死亡程序，可引发炎症
#)
    
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
