#2026/5/24
#构建共表达网络
#主要细胞群对象是Macro_MS4A6A



# Input: 
#   1.Preparing_metacells.rds

# Output: 
#   1.soft_power_test_ESCC.pdf    检查soft power
#   2.Dendrogram_ESCC.pdf       模块网络构建成功
#   3.hdWGCNA_ESCC.rds        共表达网络
#   4.kMEs_ESCC.pdf       在每个模块内，与该模块联系最紧密的几个gene的“紧密度”，kME
#   5.hMEs_ESCC.rds       hub gene 关键基因
#   6.hMEs_per_module_ESCC.pdf      每个模块的hub gene
#   7.dotplot_ESCC.pdf



suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
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



#载入数据
seurat_obj <- readRDS(file.path(data_path, "Preparing_metacells_ESCC.rds"))



#处理表达阵
Idents(seurat_obj) <- "Minor_Cell_Type"
seurat_obj <- SetDatExpr(   #默认使用metacell
    seurat_obj,
    assay = 'SCT',      # 3. 哪个测序层：使用最基础的 'RNA' 测序数据（而不是其他模态）
    layer = 'data' # 4. 数据状态：使用 'data' 层（代表已经经过 LogNormalize 归一化后的干净数据）
)

#选择soft-power阈值，soft-power的基本原理是将相关系数做一个β的幂
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'    # you can also use "unsigned" or "signed hybrid"
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
p <- wrap_plots(plot_list, ncol=2)
ggsave(p, filename = file.path(plot_path, "soft_power_test_ESCC.pdf"))
#补充：真实的生物学体内基因网络，绝对不是像随机网络那样人人平等的。它符合“二八定律”，属于无尺度网络：
#网络中极少数的核心基因（称为 Hub 基因，就像机场的交通枢纽）拥有成百上千个连接。
#绝大多数的普通基因只拥有寥寥数个连接。所以最优阈值就是使网络最接近无尺度网络
power_table <- GetPowerTable(seurat_obj)
head(power_table)



#构建共表达网络
cat("开始构建共表达网络", "\n")     #将表达模式相似的基因聚类
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'ESCC' ,# 给TOM矩阵大文件起个名字,并且被保存下来。用GetTOM()函数可以重新读取进来：TOM <- GetTOM(seurat_obj)
  overwrite_tom = T  #是否覆盖已经保存下来的TOM矩阵
)
#这一步依据soft power去计算TOM（Topological Overlap Matrix，拓扑重叠矩阵）
#TOM 不仅看基因 A 和基因 B 是不是直接相关，还会看它们是不是拥有共同的“基因朋友”。
#如果两个基因的“朋友圈”高度重合，它们就会被判定为有极强的拓扑重叠性，从而被紧紧绑在一起。

#可视化Dendrogram基因层级聚类树   不是ggplot格式，不能用ggsave来保存
pdf(file = file.path(plot_path, "Dendrogram_ESCC.pdf"), width = 10, height = 7)
PlotDendrogram(seurat_obj, main = 'ESCC hdWGCNA Dendrogram')
dev.off()
saveRDS(seurat_obj, file = file.path(data_path, 'hdWGCNA_ESCC.rds'))
cat("构建共表达网络完成，中间文件已经保存到：", file.path(data_path, "hdWGCNA_ESCC.rds"))



#计算模块特征基因与同质化
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

#将每个模块的所有基因加和，再求第一主成分（降维至一维）
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="orig.ident",
 verbose = F
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
  seurat_obj
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "ESCC_Module"
)
# plot genes ranked by kME for each module    #总共7个module，两行四列排列
p <- PlotKMEs(seurat_obj, ncol = 4)
ggsave(p, filename = file.path(plot_path, "kMEs_ESCC.pdf"), width = 10, height = 3)



# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
saveRDS(seurat_obj, file = file.path(data_path, 'hMEs_ESCC.rds'))



# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  alpha = 0.4,
  point_size = 0.1,
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
p <- wrap_plots(plot_list, ncol = 4)
ggsave(p, filename = file.path(plot_path, "hMEs_per_module_ESCC.pdf"), width = 16, height = 8)



# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)



# 创建更美观的点状图
p <- DotPlot(
  seurat_obj, 
  features = mods, 
  cols = c("lightblue", "darkblue"),  # 更专业的配色
  dot.scale = 6,  # 调整点的大小缩放
  cluster.idents = F  # 按细胞类型聚类排序
) +
  RotatedAxis() +
  theme_classic() +  # 使用经典主题，更简洁
  scale_color_gradient2(
    low = "blue",      # 低表达用蓝色
    mid = "white",     # 中等表达用白色
    high = "red",      # 高表达用红色
    midpoint = 0,      # 中点设为0
    name = "Module\nExpression"  # 更清晰的颜色图例标题
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # 旋转x轴标签
    axis.text.y = element_text(size = 12),  # y轴文本大小
    axis.title = element_blank(),  # 移除轴标题
    legend.position = "right",  # 图例位置
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题样式
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 主标题样式
    panel.grid.major = element_line(color = "grey90", size = 0.5),  # 添加轻微网格线
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),  # 分组标签背景
    strip.text = element_text(size = 11, face = "bold")
  )

# 保存高质量图片
ggsave(
  plot = p, 
  filename = file.path(plot_path, 'dotplot_ESCC.pdf'), 
  width = 7, 
  height = 8,
  dpi = 300,  # 高分辨率
  device = "pdf",
  bg = "white"  # 白色背景
)



#画雷达图，表现各module在不同类细胞中的分布
plot_list <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'Minor_Cell_Type',
  axis.label.size=4,
  grid.label.size=4, 
  grid.min = 0,     # 网格中心点设为 0
  grid.mid = 1.5,   # 网格中间线
  grid.max = 3.0,    # 将最大圈圈
  combine = F
)
p <- wrap_plots(plot_list, ncol = 4, widths = c(4,4,4,4))
print(p)
ggsave(p, filename = file.path(plot_path, "radar_ESCC.pdf"), width = 18, height = 7)
