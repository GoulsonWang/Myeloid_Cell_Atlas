#2025/10/20     搜集经典marker，查看其分布



# Input: 
#   1.SData_Marker_Processed2.RDS        Marker_Processed1.R的output

# Output: 
#   1.B_Plasma_Marker_Distribution.pdf      7种细胞大类的Marker分布图
#   2.Epithelial_Marker_Distribution.pdf
#   3.T_NK_Marker_Distribution.pdf
#   4.Fibroblast_Marker_Distribution.pdf
#   5.Macrophage_Marker_Distribution.pdf
#   6.Mast_Marker_Distribution.pdf
#   7.Endothelial_Marker_Distribution.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#Marker设置：

#来源于文章1：ESCC数据论文。共计11类，
#B:MS4A1,CD79A,CD79B
#Endothelial:CDH4,PECAM1,VWF,CLDN5
#Epithelial:KRT5,KRT14,EPCAM,SPRR1A,SCGB3A1
#Erythrocyte（红细胞）:HBB,HBA1,ALAS2,CA1
#Fibroblast:LUM,DCN,COL1A1,COL1A2
#Mononuclear phagocytes（单核吞噬细胞）:LYZ,CD1C,MRC1,CD68,CD14
#Mast（肥大细胞）:TPSAB1,TPSB2,CPA3
#Mural（壁细胞）:ACTA2,TAGLN,MYLK,MYH11,RGS5,MCAM
#Plasma（浆细胞/效应B细胞）:JCHAIN,MZB1,IGHG1,IGHA1
#Schwann（施万细胞）:S100B,CRYAB,MPZ,PLP1,SOX10
#T:CD2,CD3D,CD3E,CD3G,TRAC,TRBC1

#来源于文章2：PC数据论文，共计10类
#Epithelial:KRT19
#CD4/CD8:CD3E
#NK:GNLY(表达量高，但比例不高)
#Marcophage:APOE
#Myeloid-derived suppressor cells/Monocyte(巨噬细胞、树突状细胞、粒细胞的前提)：FCN1
#dendritic cell树突状细胞:CD1C(效果较差，建议和Marcophage或Monocyte合并)
#Mesenchyme(间质细胞)：DCN
#B/Plasma:CD79A(和Epithelial紧贴着)
#Endothelial:PECAM1(效果较差，和Epithelial紧贴着)

#来源于文章3：NSCLC数据论文，共计9类
#Epithelial:EPCAM
#Stromal cell:COL1A1,VWF(该类细胞很少)。注意:Stromal包括Fibroblast和Endothelial两大块
#T:CD3E
#B/Plasma:CD79A。其中B:MS4A1(CD20),Plasma:IGHG1
#Neutrophil:CSF3R(在Myeloid中也有少量表达)
#Myeloid：LYZ(在Neutrophil中也有少量表达)
#NK:FGFBP2(数量很少，和T细胞紧贴)
#Mast:KIT
#细胞增殖(T和Myeloid的中间态):MKI67
#pDC浆细胞样树突状细胞:LILRA4数量很少，和B紧贴着
#免疫细胞整体：CD45(PTPRC)

#来源于文章：Tumor-associated macrophage clusters linked to immunotherapy in a pan-cancer census）
#巨噬细胞Marker：CD14, CD16, CD163, CD68    

#来源于文章：A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells
#髓系细胞Marker：PTPRC和CST3

#来源于文章4：Fibroblast atlas。共计14类
#B：CD79A(在Plasma中也有一些表达)
#T/NK:CD3D,CD8A,FOXP3,KLRF1,
#Macrophage/DC/Monocyte：CD68、FCN1
#Endothelial cell:VWF,PECAME，PLVAP，FABP4（在巨噬细胞、上皮、成纤维细胞中也有部分少量表达）
#Epithelial:EPCAM
#Fibroblast:PDGFRA,COL1A2,COL3A1,PGS4
#Mast:TPSAB1
#Plasma:MZB1
#Platelet:PPBP

#来源于文章5：Pan-cancer brain metastases。共计7类
#B：CD79A，CD79B
#T/NK：CD3D，CD3E，KLRD1（表达很弱）
#Myeloid：LYZ、AIF1
#Fibroblast/Mural：TAGLN、DCN
#Endothelial：RAMP2,CLDN5
#转移灶肿瘤上皮细胞：EPCAM,KRT8,KRT19
#Glial cell神经胶状细胞：TF,PLP1

#来源于文章6：同济泛癌



#！！！在本文中所尝试的细胞类型及Marker
#B/Plasma:
#   1.(MS4A1,CD79A,CD79B,JCHAIN,MZB1,IGHG1,IGHA1)
#   2.(CD79A)
#   3.(CD79A。其中B:MS4A1(CD20),Plasma:IGHG1)
#   4.(B:CD79A,Plasma:MZB1)
#   5.(CD79A，CD79B)
#   6.CD79A,MS4A1,CD19,JCHAIN,IGKC
#测试：CD79A,CD79B,MZB1,CD79B,MS4A1(CD20),IGHG1,JCHAIN。都出现了至少两次

#Epithelial/Cancer:
#   1.KRT5,KRT14,EPCAM,SPRR1A,SCGB3A1
#   2.KRT19
#   3.EPCAM
#   4.EPCAM
#   5.EPCAM,KRT8,KRT19
#   6.KRT16,KRT19
#测试：EPCAM,KRT19

#T/NK:
#   1.CD2,CD3D,CD3E,CD3G,TRAC,TRBC1
#   2.CD3E,GNLY
#   3.CD3E,FGFBP2
#   4.CD3D,CD8A,FOXP3,KLRF1,
#   5.CD3D,CD3E，KLRD1（表达很弱）
#   6.PTPRC,CD3D,CD3E,CD8A,CD8B,CD4,CD28,FOXP3,IL2RA,KLRC1,KLRD1,MK267
#测试：CD3E,CD3D,CD8A,FOXP3,KLRD1

#Fibroblast/Mural:
#   1.Fibroblast:LUM,DCN,COL1A1,COL1A2.     Mural:ACTA2,TAGLN,MYLK,MYH11,RGS5,MCAM
#   2.缺失
#   3.缺失
#   4.Fibroblast:PDGFRA,COL1A2,COL3A1,PGS4
#   5.TAGLN、DCN
#   6.COL1A2,MYL9,MYH11,MYLK
#测试：DCN,COL1A2,TAGLN

#Macrophage/dendritic/Monocyte:
#   1.LYZ,CD1C,MRC1,CD68,CD14
#   2.APOE,FCN1,CD1C
#   3.Neutrophil:CSF3R   Myeloid：LYZ   dendritic:LILRA4        T和Myeloid的中间态:MKI67
#   4.CD68、FCN1
#   5.LYZ、AIF1
#   6.CD14,CD68,CSF1R,CD1C,CD86,FLT3,S100A9
#测试：LYZ,CD1C,CD68,FCN1,MKI67         CD14没测试

#Mast
#   1.TPSAB1,TPSB2,CPA3
#   2.缺失
#   3.KIT
#   4.TPSAB1
#   5.缺失
#   6.CPA3
#测试：TPSAB1,CPA3,KIT,TPSB2

#Endothelial
#   1.CDH4,PECAM1,VWF,CLDN5
#   2.PECAM1
#   3.COL1A1,VWF
#   4.VWF,PECAME，PLVAP，FABP4
#   5.RAMP2,CLDN5
#   6.VWF,PECAM1
#测试：PECAM1,VWF,CLDN5



#载入数据
SData <- readRDS(file.path(data_path, "SData_Marker_Processed2.RDS"))
#dim(SData)     #44296  5000



#Marker设置
B_Plasma_Markers <- c("CD79A","CD79B","MZB1","MS4A1","IGHG1","JCHAIN")
Epithelial_Markers <- c("EPCAM","KRT19")
T_NK_Markers <- c("CD3E","CD3D","CD8A","FOXP3","KLRD1")
Fibroblast_Markers <- c("DCN","COL1A2","TAGLN")
Macrophage_Markers <- c("LYZ","CD1C","CD68","FCN1","MKI67")
Mast_Markers <- c("TPSAB1","CPA3","KIT","TPSB2")
Endothelial_Markers <- c("PECAM1","VWF","CLDN5")



#绘制各Marker的分布图
p1 <- FeaturePlot(SData,
    features = B_Plasma_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "B_Plasma_Marker_Distribution.pdf"),
    plot = p1,
    width = 12, # 减小宽度
    height = 8, # 减小高度
    units = "in",
    dpi = 300
)

p2 <- FeaturePlot(SData,
    features = Epithelial_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 2, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Epithelial_Marker_Distribution.pdf"),
    plot = p2,
    width = 8, # 减小宽度
    height = 4, # 减小高度
    units = "in",
    dpi = 300
)

p3 <- FeaturePlot(SData,
    features = T_NK_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "T_NK_Marker_Distribution.pdf"),
    plot = p3,
    width = 12, # 减小宽度
    height = 8, # 减小高度
    units = "in",
    dpi = 300
)

p4 <- FeaturePlot(SData,
    features = Fibroblast_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Fibroblast_Marker_Distribution.pdf"),
    plot = p4,
    width = 12, # 减小宽度
    height = 4, # 减小高度
    units = "in",
    dpi = 300
)

p5 <- FeaturePlot(SData,
    features = Macrophage_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Macrophage_Marker_Distribution.pdf"),
    plot = p5,
    width = 12, # 减小宽度
    height = 8, # 减小高度
    units = "in",
    dpi = 300
)

p6 <- FeaturePlot(SData,
    features = Mast_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 2, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Mast_Marker_Distribution.pdf"),
    plot = p6,
    width = 8, # 减小宽度
    height = 8, # 减小高度
    units = "in",
    dpi = 300
)

p7 <- FeaturePlot(SData,
    features = Endothelial_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    label = TRUE, # 可选：标注基因名
    label.size = 3, # 标签字体大小，减小字体
    ncol = 3, 
    min.cutoff = "q9",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Endothelial_Marker_Distribution.pdf"),
    plot = p7,
    width = 12, # 减小宽度
    height = 4, # 减小高度
    units = "in",
    dpi = 300
)



#注释结果：
#T/NK:Cluster1
#Fibroblast:Cluster7
#Macrophage：Cluster2
#Mast:Cluster10
#Endothelial:Cluster8
#B:Cluster4
#Plasma:Cluster9
#Epithelial:Cluster0
#存疑的Cluster5：大部分表达Epithelial，下沿右部分表达MKI67（增值状态）



#结合Marker_Processed2.R的Marker在CellMarker数据库中的检索结果：
#CellMarker注释结果：
#Cluster0：最高得分仅0.14
#Cluster1：同时3个0.5分，分别是CD8+T cell,Naive cytotoxic Tcell,分别依靠PTPRC,CD3D,CD7
#Cluster2：满分，Pro-inflammatory macrophage，依靠CD14
#Cluster4：最高得分0.2，Activated effector cell，依靠HLA-DRA
#Cluster5：同时两个满分——Proliferative cell，Cycling cell，分别依靠MKI67；MKI67和TOP2A
#Cluster7：最高得分仅0.05
#Cluster8：最高得分仅0.12
#Cluster9：最高得分0.29，Plasma cell，依靠CD79A,IGHG1,IGHG3,XBP1,JCHAIN
#Cluster10：最高得分0.5，Mast cell，依靠KIT，TPSAB1



#最终判定Cluster5为Profilerating