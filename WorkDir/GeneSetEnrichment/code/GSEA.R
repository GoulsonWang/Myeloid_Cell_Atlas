#2025/12/12
#本脚本用与GSEA



# Input: 
#   1.FCList_Macro.rds
#   2.FCList_DC.rds
#   3.FCList_Mono.rds
#   4.FCList_MastNeut.rds
#   5.SCT_HVG_3000_genes.csv

# Output: 
#   1.KEGG_GSEA", subclusterID[i], ".pdf"
#   2.KEGG_GSEA_Result.rds



suppressMessages((library(clusterProfiler)))
suppressMessages((library(ggplot2)))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据——差异基因
Gene_Macro <- readRDS(file.path(data_path, "FCList_Macro.rds"))
Gene_Mono <- readRDS(file.path(data_path, "FCList_Mono.rds"))
Gene_DC <- readRDS(file.path(data_path, "FCList_DC.rds"))
Gene_MastNeut <- readRDS(file.path(data_path, "FCList_MastNeut.rds"))
#str(Gene_Macro)     #list，每个子list为一个亚类，内包含data.frame，分别为Gene,pct.1,pct.2,avg_log2FC



#基因的筛选
Diff_Gene <- c(Gene_Macro, Gene_Mono, Gene_DC, Gene_MastNeut)
subclusterID <- names(Diff_Gene)        
#"Macro_1"  "Macro_2"  "Macro_3"  "Macro_4"  "Macro_5"  
#"Mono_1"   "Mono_2"   "DC_1"     "DC_2"     "DC_3"     "Mast"    
#"Neutro_1"



#GSEA富集
KEGG_GSEA_Result <- list()
for(i in 1:length(subclusterID)){
    # 对筛选后的基因按avg_log2FC从大到小排序
    sorted_genes <- Diff_Gene[[i]][order(Diff_Gene[[i]]$avg_log2FC, decreasing = TRUE), ]

    # 转换基因符号为ENTREZID
    GeneID_kegg <- bitr(rownames(sorted_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    if(nrow(GeneID_kegg) > 0) {
      # 创建基因列表，名称为ENTREZID，值为avg_log2FC
      geneList <- setNames(sorted_genes[GeneID_kegg$SYMBOL, "avg_log2FC"], GeneID_kegg$ENTREZID)
      
      # 确保基因列表按照FC值从大到小排序
      geneList <- sort(geneList, decreasing = TRUE)
      
      # 移除任何NA值
      geneList_clean <- geneList[!is.na(names(geneList))]
      
      if(length(geneList_clean) > 0) {
        KEGG_GSEA_Result[[i]] <- gseKEGG(
            gene = geneList_clean,
            organism = "hsa",
            keyType = "kegg",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.2,          
            minGSSize = 10,
            maxGSSize = 500,
            verbose = F
        )
        
        # 检查是否富集出通路
        if(nrow(KEGG_GSEA_Result[[i]]) > 0) {
          p1 <- dotplot(KEGG_GSEA_Result[[i]], showCategory=30) + ggtitle(paste0("dotplot for GSEA", ":", subclusterID[i]))
          ggplot2::ggsave(plot = p1, filename = paste(plot_path, "/", "KEGG_GSEA", subclusterID[i], ".pdf", sep = ""), width = 10, height = 10)
          cat("已保存", subclusterID[i], "的GSEA结果图\n")
        } else {
          cat("注意：", subclusterID[i], "未富集出任何通路\n")
        }
      } else {
        cat("Warning: No valid genes for", subclusterID[i], "\n")
        KEGG_GSEA_Result[[i]] <- NULL
      }
    } else {
      cat("Warning: No genes could be converted to ENTREZID for", subclusterID[i], "\n")
      KEGG_GSEA_Result[[i]] <- NULL
    }
}

saveRDS(KEGG_GSEA_Result, file = paste(data_path, "/KEGG_GSEA_Result.rds", sep = ""))



