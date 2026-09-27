#本脚本用于学习hdWGCNA的流程


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/hdWGCNA/code")
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/hdWGCNA"
code_path <- file.path(file_path, "code")
data_path <- file.path(file_path, "data")
plot_path <- file.path(file_path, "plot")


# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading  启用8个核
enableWGCNAThreads(nThreads = 8)

# load the Zhou et al snRNA-seq dataset
#seurat_obj <- readRDS('/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/hdWGCNA/data/Zhou_2020_control.rds')
#p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
#   umap_theme() + ggtitle('Zhou et al Control Cortex') + NoLegend()
#ggsave(p, filename = file.path(plot_path, "UMAP_Example.pdf"))



#为WGCNA准备Seurat变量

#seurat_obj <- SetupForWGCNA(
#  seurat_obj,
#  gene_select = "fraction",     # 选择计算时所用的基因。
  #fraction：选择在整个数据集中只在大部分细胞中表达的基因
  #variable：选择VariableFeatures高变基因
  #custom：使用自定义基因集
#  fraction = 0.05,              # 至少0.05比例的细胞表达了该基因
#  wgcna_name = "tutorial"       # the name of the hdWGCNA experiment
#)
#hdWGCNA的计算结果储存在Seurat@misc中


#构建metacells:识别表达谱相似的细胞（聚类），将一类细胞的表达谱聚合加起来或求均值

#这一步的主要目的是为了处理稀疏性
#处理时必须保证一个metacells中的细胞数据都来自同一个样本
#metacells所包含细胞的数量推荐设置为20-75
#参数 max_shared：用来控制元细胞之间允许“共享”多少个重复的单细胞（即允许重叠的程度），用以调节平滑度
#min_cells 参数，设置一个硬性门槛（比如若少于 50 个细胞就直接忽略该组），如果这个值设得太低，后续分析非常容易报错卡死。
#seurat_obj <- MetacellsByGroups(
#  seurat_obj = seurat_obj,
#  group.by = c("cell_type", "Sample"), # specify the columns in seurat_obj@meta.data to group by
#  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
#  k = 25, # nearest-neighbors parameter
#  max_shared = 10, # maximum number of shared cells between two metacells
#  ident.group = 'cell_type' # set the Idents of the metacell seurat object
#)
#saveRDS(seurat_obj, file.path(data_path, "Construct_metacells.rds"))



#构建好metacell后，可以单独把metacell提取出来重新进行预处理，也可以保存在原Seurat对象中重新预处理
#主要目的是根据umap图检查metacell聚合的好不好
#提取metacell
#metacell_obj <- GetMetacellObject(seurat_obj)
#cat("开始处理metacell")
#seurat_obj <- NormalizeMetacells(seurat_obj)
#seurat_obj <- Seurat::FindVariableFeatures(
#  seurat_obj, 
#  selection.method = "vst", 
#  nfeatures = 2000  # 或者根据你的需要调整高变基因数量，通常 2000-3000
#)
#seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
#seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
#seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='Sample')
#seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)


#p1 <- DimPlotMetacells(seurat_obj, group.by='cell_type') + umap_theme() + ggtitle("Cell Type")
#p2 <- DimPlotMetacells(seurat_obj, group.by='Sample') + umap_theme() + ggtitle("Sample")

#p <- p1 + p2
#ggsave(p, filename = file.path(plot_path, "UMAP_MetaCell.pdf"))
#saveRDS(seurat_obj, file.path(data_path, "Preparing_metacells.rds"))
#cat("处理metacell完成，中间文件已经保存到：", file.path(plot_path, "Preparing_metacells.rds"))
seurat_obj <- readRDS(file.path(data_path, "Preparing_metacells.rds"))


#共表达网络分析

#处理表达阵
seurat_obj <- SetDatExpr(   #默认使用metacell
  seurat_obj,
  # 1. 目标群体：我这次只想研究 "INH"（抑制性神经元）这一类细胞
  group_name = "INH", 
  # 2. 分类依据：去 metadata 的 'cell_type' 这一列里找上面说的 "INH"
  group.by = 'cell_type', 
  # 3. 哪个测序层：使用最基础的 'RNA' 测序数据（而不是其他模态）
  assay = 'RNA', 
  # 4. 数据状态：使用 'data' 层（代表已经经过 LogNormalize 归一化后的干净数据）
  # 提示：如果你之前用的是 SCTransform 标准化，这里可以改成 layer = 'data' 并指定 assay = 'SCT'
  layer = 'data' 
)
#挑选多个组时
#seurat_obj <- SetDatExpr(
#  seurat_obj,
#  group_name = c("INH", "EX"),
#  group.by='cell_type'
#)

#选择soft-power阈值，soft-power的基本原理是将相关系数做一个β的幂
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
p <- wrap_plots(plot_list, ncol=2)
ggsave(p, filename = file.path(plot_path, "soft_power_test.pdf"))
#补充：真实的生物学体内基因网络，绝对不是像随机网络那样人人平等的。它符合“二八定律”，属于无尺度网络：
#网络中极少数的核心基因（称为 Hub 基因，就像机场的交通枢纽）拥有成百上千个连接。
#绝大多数的普通基因只拥有寥寥数个连接。所以最优阈值就是使网络最接近无尺度网络
power_table <- GetPowerTable(seurat_obj)
head(power_table)

#构建共表达网络
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'INH' ,# 给TOM矩阵大文件起个名字,并且被保存下来。用GetTOM()函数可以重新读取进来：TOM <- GetTOM(seurat_obj)
  overwrite_tom = TRUE  #是否覆盖已经保存下来的TOM矩阵
)
#这一步依据soft power去计算TOP（Topological Overlap Matrix，拓扑重叠矩阵）
#TOM 不仅看基因 A 和基因 B 是不是直接相关，还会看它们是不是拥有共同的“基因朋友”。
#如果两个基因的“朋友圈”高度重合，它们就会被判定为有极强的拓扑重叠性，从而被紧紧绑在一起。

#可视化Dendrogram基因层级聚类树   不是ggplot格式，不能用ggsave来保存
pdf(file = file.path(plot_path, "Dendrogram.pdf"), width = 10, height = 7)
PlotDendrogram(seurat_obj, main = 'INH hdWGCNA Dendrogram')
dev.off()



#计算模块特征基因与同质化
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="Sample"
)
#我们不可能拿这 500 个基因一个个去跟患者的 pCR 算相关性，这会导致极其严重的统计学多重假设检验偏倚。
#我们需要把这 500 个基因的整体表达趋势浓缩成一个代表性的数值。这个数值就是 Module Eigengene (ME)。
#底层数学原理：函数会对属于该模块的所有基因的表达矩阵进行主成分分析（PCA）。提取出来的第一主成分（PC1）就是这个模块的 ME。
#生物学含义：ME 是这个模块所有基因表达模式的“综合代言人”。如果某个细胞或样本的“蓝色模块 ME”值很高，就意味着这个模块里绝大多数基因在这个细胞里都处于高表达状态。
#hdWGNCA 引入了 ModuleEigengenes 函数，它的一大杀器就是直接在特征基因级别调用 Harmony 算法进行批次矫正，从而输出和谐后的模块特征基因（harmonized Module Eigengenes, hMEs）。
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)



# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'INH'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "INH-M"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)
ggsave(p, filename = file.path(plot_path, "kMEs.pdf"))



# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
saveRDS(seurat_obj, file = file.path(data_path, 'hdWGCNA_object.rds'))



# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)



# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
p <- wrap_plots(plot_list, ncol=6)
ggsave(p, filename = file.path(plot_path, "hMEs_per_module.pdf"))

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
p <- wrap_plots(plot_list, ncol=6)
ggsave(p, filename = file.path(plot_path, "hub_scores.pdf"))

seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$annotation), ' '))[,1]

p <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% subset(cell_type == 'INH') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)
ggsave(p, filename = file.path(plot_path, "hub_radar.pdf"))



# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)



# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
ggsave(plot=p, filename='dotplot.pdf', width=10, height=10)
