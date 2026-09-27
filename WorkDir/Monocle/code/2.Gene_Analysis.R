#2025/1/7  

#对不同的branch细胞进行检验，发现其中随时间显著变化的基因



#input
#   1.gene_fits.RDS     广义线性模型的回归结果
#   2.CDS_VariableGene_Only   仅仅包含3000个高变基因的cds对象

#output
#   1.branch_analysis_result.RDS
#   2.trajectory_BranchName.pdf



suppressMessages(library(monocle3))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(dplyr))
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
cds <- load_monocle_objects(directory_path = file.path(data_path, "CDS_VariableGene_Only"))



#Branch Analysis
#设置函数，输入branch起始和终止节点，输入branch_name，输出对应branch的umap图和graph_test()结果
select_trajectory_and_test <- function(cds, 
                                   starting_pr_node, 
                                   ending_pr_nodes,  
                                   branch_name,
                                   reduction_method = "UMAP",
                                   cores = 5) {
  # 选择轨迹上的细胞
  branch_result <- choose_graph_segments(
    cds, 
    reduction_method = reduction_method,
    starting_pr_node = starting_pr_node,
    ending_pr_nodes = ending_pr_nodes, 
    clear_cds = FALSE, 
    return_list = F
  )
  
  selected_cells <- colnames(branch_result)
  # 创建一个临时列来标记选择的细胞
  temp_col_name <- paste0("temp_", branch_name)
  colData(cds)[[temp_col_name]] <- "Unchosen"
  colData(cds)[[temp_col_name]][colnames(cds) %in% selected_cells] <- branch_name
  # 绘制轨迹图
  p <- plot_cells(cds,
                  color_cells_by = temp_col_name,
                  label_groups_by_cluster = FALSE,
                  label_cell_groups = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 5,
                  trajectory_graph_segment_size = 0.5) +
    scale_color_manual(values = c("purple", "grey"))
  # 保存图片
  plot_file_path <- file.path(plot_path, paste0("trajectory_", branch_name, ".pdf"))
  ggsave(plot_file_path, p, height = 5, width = 6)
  # 执行图测试
  test_result <- graph_test(
    branch_result,
    neighbor_graph = "knn",  # 使用knn而不是principal_graph以避免错误
    cores = cores
  )
  # 返回测试结果
  return(list(
    branch_cells = colnames(branch_result),
    test_result = test_result,
    plot_path = plot_file_path
  ))
}
branch_analysis <- list()
branch_analysis[["Trunk"]] <- select_trajectory_and_test(cds, "Y_149", c("Y_9"), "Trunk")
branch_analysis[["MS4A6A_1"]] <- select_trajectory_and_test(cds, "Y_250", c("Y_112", "Y_196", "Y_38"), "MS4A6A_1")
branch_analysis[["MS4A6A_2"]] <- select_trajectory_and_test(cds, "Y_9", c("Y_96", "Y_108", "Y_298", "Y_189", "Y_275"), "MS4A6A_2")
branch_analysis[["APOE_1"]] <- select_trajectory_and_test(cds, "Y_9", c("Y_152"), "APOE_1")
branch_analysis[["APOE_2"]] <- select_trajectory_and_test(cds, "Y_125", c("Y_220", "Y_217", "Y_230", "Y_129", "Y_276"), "APOE_2")
saveRDS(branch_analysis, file = file.path(data_path, "branch_analysis_result.RDS"))

#另外：将五条branch的通路画在同一张图里‘
branch_name <- c("Trunk", "MS4A6A_1", "MS4A6A_2", "APOE_1", "APOE_2")
colData(cds)[["branch"]] <- "Unchosen"
for(i in branch_name){
  if(i %in% names(branch_analysis) && !is.null(branch_analysis[[i]]$branch_cells)) {
    colData(cds)[["branch"]][colnames(cds) %in% branch_analysis[[i]]$branch_cells] <- i # 赋值为当前循环的分支名 i
    print(paste("Assigned", length(branch_analysis[[i]]$branch_cells), "cells to", i)) # 可选：打印信息确认
  } else {
    print(paste("Warning: branch_analysis does not contain element", i, "or branch_cells is NULL.")) # 警告信息
  }
}
# 绘制轨迹图
p <- plot_cells(
  cds,
  color_cells_by = "branch",
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 5,
  trajectory_graph_segment_size = 0.5) +
  scale_color_manual(values = c("#F8766D", "#F97600", "#7CAE00", "#00BFC4", "#C77CFF", "grey"))
  # 保存图片
plot_file_path <- file.path(plot_path, paste0("trajectory_branch_all", ".pdf"))
ggsave(plot_file_path, p, height = 5, width = 6)



#查看branch_analysis_result.RDS
branch_analysis <- readRDS(file.path(data_path, "branch_analysis_result.RDS"))
