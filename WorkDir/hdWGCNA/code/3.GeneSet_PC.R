#2026/5/30
#本脚本用于绘制hub_gene图、不同通路的基因集在各module的富集情况
#预期结果：按照4.FindDME.R的结果：
#       PC中Module6，Module7与Mast高度相关，Module7治疗后下调，Module6上调
#       PC中Module,Module2与MS4A6A高度相关，Module2治疗后下调，module1上调



# Input: 
#   1.hMEs_PC.rds

# Output: 
#   1.1.hub_gene_combined_PC.pdf



suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidygraph))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/hdWGCNA")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/hdWGCNA"
code_path <- file.path(file_path, "code")
data_path <- file.path(file_path, "data")
plot_path <- file.path(file_path, "plot")



#常规设置
theme_set(theme_cowplot())      ## using the cowplot theme for ggplot
set.seed(12345)         # set random seed for reproducibility
enableWGCNAThreads(nThreads = 8)        #启用8个核
options(future.globals.maxSize = 3 * 1024^3)



#载入数据
pc_data <- readRDS(file.path(data_path, "hMEs_PC.rds"))



# 处理 PC 数据

# 1.可视化 PC 数据的所有模块的 hub gene
pdf(file = file.path(plot_path, "1.hub_gene_combined_PC.pdf"), width = 7, height = 7)
p <- HubGeneNetworkPlot(
  pc_data,
  n_hubs = 10, n_other = 20, 
  edge_prop = 0.8, 
  mods = c('PC_Module6', "PC_Module7"), 
  hub.vertex.size = 5
)
dev.off()

hub_df <- GetHubGenes(pc_data, n_hubs = 10)
hub_df%>% filter(module %in% c("PC_Module6", "PC_Module7"))
hub_df$gene_name
