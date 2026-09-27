#2026/5/24      
#本脚本用于做准备工作（为了构建基因调控网络）：metacell的构建与预处理



# Input: 
#   1./home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_Rename.RDS

# Output: 
#   1.UMAP_Minor_Cell_Type_ESCC.pdf      检查input数据
#   2.Construct_metacells_ESCC.rds       构建好的metacell对象
#   3.Preparing_metacells_ESCC.rds       处理后的metacell对象



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
#seurat_obj <- readRDS('/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_Rename.RDS')
#seurat_obj <- subset(seurat_obj, subset = CancerType == "ESCC")
#Idents(seurat_obj) <- "Minor_Cell_Type"
#p <- DimPlot(seurat_obj, label=TRUE) + umap_theme()
#ggsave(p, filename = file.path(plot_path, "UMAP_Minor_Cell_Type_ESCC.pdf"))



#为WGCNA准备Seurat变量
#seurat_obj <- SetupForWGCNA(
#    seurat_obj,
#    gene_select = "fraction",     # 选择计算时所用的基因。
        #fraction：选择在整个数据集中只在大部分细胞中表达的基因
        #variable：选择VariableFeatures高变基因
        #custom：使用自定义基因集
#    fraction = 0.05,              # 至少0.05比例的细胞表达了该基因
#    wgcna_name = "Myeloid"       # the name of the hdWGCNA experiment
#)
#hdWGCNA的计算结果储存在Seurat@misc中



#构建metacells:识别表达谱相似的细胞（聚类），将一类细胞的表达谱聚合加起来或求均值

#这一步的主要目的是为了处理稀疏性
#处理时必须保证一个metacells中的细胞数据都来自同一个样本
#metacells所包含细胞的数量推荐设置为20-75
#参数 max_shared：用来控制元细胞之间允许“共享”多少个重复的单细胞（即允许重叠的程度），用以调节平滑度
#min_cells 参数，设置一个硬性门槛（比如若少于 50 个细胞就直接忽略该组），如果这个值设得太低，后续分析非常容易报错卡死。
#cat("开始构建metacell", "\n")
#seurat_obj <- MetacellsByGroups(
#    seurat_obj = seurat_obj,
#    group.by = c("Minor_Cell_Type", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
#    reduction = 'pca',      # select the dimensionality reduction to perform KNN on
#    k = 25,                 # nearest-neighbors parameter
#    max_shared = 10,            # maximum number of shared cells between two metacells
#    ident.group = 'Minor_Cell_Type'       # set the Idents of the metacell seurat object
#)
#saveRDS(seurat_obj, file.path(data_path, "Construct_metacells_ESCC.rds"))
#cat("metacell构建完成，中间文件已经保存到：", file.path(data_path, "Construct_metacells_ESCC.rds"))



#构建好metacell后，可以单独把metacell提取出来重新进行预处理，也可以保存在原Seurat对象中重新预处理
#主要目的是根据umap图检查metacell聚合的好不好
#提取metacell
seurat_obj <- readRDS(file.path(data_path, "Construct_metacells_ESCC.rds"))
cat("开始处理metacell", "\n")
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)


p1 <- DimPlotMetacells(seurat_obj, group.by='Minor_Cell_Type', alpha = 0.5, size = 0.8) + 
  hdWGCNA::umap_theme() +  # 使用hdWGCNA包的umap_theme函数
  theme(
    # 坐标轴标题调整
    axis.title = element_text(size = 14, face = "bold"),
    # 图例设置
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    
    # 标题设置
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    # 整体布局
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")
print(p1)
ggsave(p1, filename = file.path(plot_path, "UMAP_MetaCell_ESCC.pdf"), width = 10, height = 8)
saveRDS(seurat_obj, file.path(data_path, "Preparing_metacells_ESCC.rds"))
cat("处理metacell完成，中间文件已经保存到：", file.path(plot_path, "Preparing_metacells_ESCC.rds"))
