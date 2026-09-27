#2025/11/6
#本脚本用与尝试利用map来注释和聚类Myeloid亚型



#input：
#   1.SData_Integrated_rpca.RDS     整合后未聚类的数据
#   2.TabulaTIME_integrated_expression.h5   reference图谱
#   3.TabulaTIME_metacell_meta_information.rds  reference图谱的meta信息

#output：
#   1.ref_SData_subset_umap.pdf     reference图谱中关于Macrophage的umap图



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
SData <- readRDS(file.path(data_path, "SData_Integrated_rpca(Demo).RDS"))
#DemoCellName <- sample(colnames(SData), size = 5000)  
#Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)       #2000 5000
#str(Demo)
#saveRDS(Demo, paste(data_path, "SData_Integrated_rpca(Demo).RDS", sep = "/"))
library(SeuratDisk)
ref <- Read10X_h5("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Reference/TabulaTIME_integrated_expression.h5")      #返回结果是一个dgCMatrix
ref_meta <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Reference/TabulaTIME_metacell_meta_information.rds")
head(ref_meta)
dim(ref_meta)   #140072     11
dim(ref)        #9031 140072
head(ref)[, 10:13]
colnames(ref)[10:13]     #"BCC_GSE123813_aPD1@bcc_su001_pre|0|1" "BCC_GSE123813_aPD1@bcc_su001_pre|0|0"  "BCC_GSE123813_aPD1@bcc_su001_pre|0|0|2"



#构建ref_SData
ref_SData <- CreateSeuratObject(ref, meta.data = ref_meta)      #4Gb
#colnames(ref_SData@meta.data)
#[1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
#[4] "UMAP1"                "UMAP2"                "curated_anno"        
#[7] "Tissue"               "Cancer_type"          "Dataset"             
#[10] "Source"               "Sample"               "curated_anno_subtype"
#table(ref_SData$Cancer_type)
#   BCC   BLCA   BRCA Breast   CESC   CHOL    CRC   ESCA   GIST     HB   HNSC 
#  3338   1006   2866   1271   2360   1552   6698   6311   1052   1517   6940 
# HNSCC   KICH Kidney  KIPAN   KIRC   LIHC   Live   LSCC   LUAD    MCC    NPC 
#  5324     87     82   1394   5583   2769    890    393    328    657   7074 
# NSCLC   Oral     OS   OSCC     OV  Ovary   PAAD   PBMC    PPB   PRAD    SCC 
# 39891    993   1421   1806   7506   1201   3148   3069    377   6157   2358 
#  SCLC   SKCM     SS   STAD   THCA   UCEC    UVM 
#    56   5296    288   1894    737    655   3727
table(ref_SData$Source)
#       Blood   Metastatic       Normal Precancerous        Tumor 
#       15332         5937        24454         2583        91766
table(ref_SData$curated_anno)
#             B       CD4Tconv           CD8T             DC    Endothelial 
#          5297          23336          35692           2387           4218 
#    Epithelial    Fibroblasts      Malignant           Mast     Mono/Macro 
#          4013           6635          22413            793          13411 
#Myofibroblasts             NK         Others         Plasma        Tprolif 
#          2206           7321             96           2631           1846 
#          Treg 
#          7777

#对ref取子集，只取Mast、DC、Mono/Macro
ref_subset_SData <- subset(ref_SData, subset = (curated_anno %in% c("DC", "Mono/Macro", "Mast") & Source == "Tumor"))
#table(ref_subset_SData$Source)
#Tumor 
# 9406 
#table(ref_subset_SData$curated_anno)
#        DC       Mast Mono/Macro 
#      1561        493       7352 
table(ref_subset_SData$curated_anno_subtype)
table(SData$PreOrPost)
table(SData$Response.RECIST)
#        Partial  Stable 
#   2054    1478    1468 
table(SData$Response.Pathologic)
#      MPR NMPR  pCR 
#1129  620 1986  935
table(SData$Response.Comprehensive)
#  No  Yes 
#2488 2182 


#Mapping
map.anchors <- FindTransferAnchors(
    reference = ref_subset_SData,  
    query = SData, 
    normalization.method = "SCT",
    dims = 1:30,
    reference.reduction = "pca")
predictions <- TransferData(
    anchorset = map.anchors, 
    refdata = ref_subset_SData$curated_anno, 
    dims = 1:30)
SData <- AddMetaData(SData, metadata = predictions, col.name = "Map1_Annotaiton")


umapdata1 <- ref_subset_SData$UMAP1
umapdata2 <- ref_subset_SData$UMAP2
celltype <- ref_subset_SData$curated_anno_subtype
ggplotdataframe <- data.frame(umapdata1, umapdata2, celltype)
p1 <- ggplot(ggplotdataframe, aes(x = umapdata1, y = umapdata2)) +
        geom_point(shape = 1, aes(color = celltype), size = 0.5, stroke = 0.3) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) # 调整图例点大小
ggsave(
    filename = paste(plot_path, "/ref_SData_subset_umap.pdf", sep = ""),
    plot = p1, 
    width = 12, 
    height = 10
)

#ref的umap图结果不好，遂放弃