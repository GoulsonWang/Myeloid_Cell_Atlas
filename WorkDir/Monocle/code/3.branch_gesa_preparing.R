#2026/1/11  对五条branch的结果进行富集分析（以FC值为基因集排序）



#input
#   1.branch_analysis_result.RDS
#   2.SData_Minor_Cell_Type.RDS

#output
#   1.FCList_Branch.RDS     四条branch中随时间显著变化的基因，相较于Trunk的FC值排序结果
#   2.SData_branch_assignment.RDS     添加了branch信息的SData



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(future))



# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/code/")

#提交时使用
plan(multicore, workers = 6)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值



branch_analysis_result <- readRDS(file.path(data_path, "branch_analysis_result.RDS"))   
#包含五条branch，其中每条branch下包括三列数据：branch相关细胞名称，3000个基因是否随时间显著变化的检验结果，和branch图保存地址
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS")



# 为SData创建一个新的metadata列，标记每个细胞属于哪个branch
branch_assignment <- rep(NA, ncol(SData))
cell_names <- colnames(SData)
for (branch_name in names(branch_analysis_result)) {
  # 获取当前分支的细胞名称
  branch_cells <- branch_analysis_result[[branch_name]]$branch_cells
  # 找到这些细胞在SData中的位置
  matching_cells <- intersect(branch_cells, cell_names)
  # 将这些细胞标记为当前分支的名称
  branch_assignment[cell_names %in% matching_cells] <- branch_name
}
SData$branch_group <- branch_assignment
table(SData$branch_group, useNA = "always")
saveRDS(SData, file = file.path(data_path, "SData_branch_assignment.RDS"))



#筛选其中q_value小于0.05的基因名称
extract_significant_genes <- function(branch_analysis_result, threshold = 0.05) {
  significant_results <- list()
  for (branch_name in names(branch_analysis_result)) {
    cat("Processing", branch_name, "\n")
    current_branch <- branch_analysis_result[[branch_name]]
    test_result <- current_branch$test_result
    sig_indices <- which(test_result$q_value < threshold)
    if (length(sig_indices) > 0) {
      significant_data <- test_result[sig_indices, ]
      significant_results[[branch_name]] <- list(
        genes = significant_data$gene_short_name,
        q_values = significant_data$q_value,
        p_values = significant_data$p_value,
        cell_count = length(current_branch$branch_cells),
        significant_gene_count = nrow(significant_data)
      )
      
      cat("  Found", nrow(significant_data), "significant genes in", branch_name, "\n")
    } else {
      cat("  No significant genes found in", branch_name, "\n")
      significant_results[[branch_name]] <- list(
        genes = character(0),
        q_values = numeric(0),
        p_values = numeric(0),
        morans_I = numeric(0),
        morans_test_statistic = numeric(0),
        cell_count = length(current_branch$branch_cells),
        significant_gene_count = 0
      )
    }
  }
  
  return(significant_results)
}

# 调用函数提取显著基因
sig_genes_result <- extract_significant_genes(branch_analysis_result, threshold = 0.05)
#Processing Trunk 
#  Found 1853 significant genes in Trunk 
#Processing MS4A6A_1 
#  Found 1323 significant genes in MS4A6A_1 
#Processing MS4A6A_2 
#  Found 1172 significant genes in MS4A6A_2 
#Processing APOE_1 
#  Found 683 significant genes in APOE_1 
#Processing APOE_2 
#  Found 2282 significant genes in APOE_2 



#计算四个分支中这些基因相较于Trunk的FC值
SData$branch_group <- as.factor(SData$branch_group)
SData@active.ident <- SData$branch_group
branch_names <- c("MS4A6A_1", "MS4A6A_2", "APOE_1", "APOE_2")
FCList <- list()
for(i in branch_names){
  significant_genes <- sig_genes_result[[i]]$genes
  FCList[[i]] <- FoldChange(
    SData, 
    ident.1 = i, 
    ident.2 = "Trunk", 
    features = significant_genes
  )
  cat(i, "的", nrow(FCList[[i]]), "个基因的FC值计算完毕\n")
  FCList[[i]] <- FCList[[i]][order(FCList[[i]]$avg_log2FC, decreasing = TRUE),]
}
saveRDS(FCList, file = file.path(data_path, "FCList_Branch.RDS"))



# 将显著基因及其FC值保存为CSV文件
for(i in branch_names){
  significant_genes <- sig_genes_result[[i]]$genes
  fc_data <- FCList[[i]]
  
  # 创建包含显著基因及其FC值的数据框
  if(nrow(fc_data) > 0 && length(significant_genes) > 0){
    # 筛选出显著基因的FC数据
    sig_fc_data <- fc_data[rownames(fc_data) %in% significant_genes, ]
    
    # 添加分支名称列
    sig_fc_data$branch <- i
    
    # 保存为CSV文件
    csv_filename <- file.path(data_path, paste0("Significant_Genes_FC_", i, ".csv"))
    write.csv(sig_fc_data, file = csv_filename, row.names = TRUE)
    cat("已保存", nrow(sig_fc_data), "个显著基因及其FC值到", csv_filename, "\n")
  }else{
    cat("分支", i, "没有显著基因数据或FC数据\n")
  }
}

# 另外，合并所有分支的数据到一个总文件中
all_branches_data <- data.frame()
for(i in branch_names){
  significant_genes <- sig_genes_result[[i]]$genes
  fc_data <- FCList[[i]]
  
  if(nrow(fc_data) > 0 && length(significant_genes) > 0){
    # 筛选出显著基因的FC数据
    sig_fc_data <- fc_data[rownames(fc_data) %in% significant_genes, ]
    
    # 添加分支名称列
    sig_fc_data$branch <- i
    
    # 添加到总数据框
    all_branches_data <- rbind(all_branches_data, sig_fc_data)
  }
}

if(nrow(all_branches_data) > 0){
  # 按照分支和平均log2FC值排序
  all_branches_data <- all_branches_data[order(all_branches_data$branch, all_branches_data$avg_log2FC, decreasing = TRUE),]
  
  # 保存总文件
  combined_csv_filename <- file.path(data_path, "All_Significant_Genes_FC_Combined.csv")
  write.csv(all_branches_data, file = combined_csv_filename, row.names = TRUE)
  cat("已保存所有分支的显著基因及其FC值到", combined_csv_filename, "\n")
  cat("总共包含", nrow(all_branches_data), "条记录\n")
}else{
  cat("没有找到任何显著基因数据用于保存\n")
}