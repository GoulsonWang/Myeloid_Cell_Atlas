#计算图谱的各分组占比


suppressMessages(library(Seurat))
library(dplyr)
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_and_Rename.RDS")
colnames(SData@meta.data)
# [1] "orig.ident"             "nCount_RNA"             "nFeature_RNA"          
# [4] "percent.MT"             "RNA_snn_res.0.1"        "seurat_clusters"       
# [7] "Doublet_Score"          "Is_Doublet"             "CancerType"            
#[10] "percent.HB"             "percent.RB"             "TorB"                  
#[13] "RNA_snn_res.0.2"        "Lineage"                "Patient.ID"            
#[16] "PreOrPost"              "PrimaryOrMet"           "Response.RECIST"       
#[19] "Response.Pathologic"    "Response.Comprehensive" "percent.HSP"           
#[22] "nCount_SCT"             "nFeature_SCT"           "SCT_snn_res.0.3"       
#[25] "Major_Cell_Type"        "Minor_Cell_Type"  

table(SData$orig.ident) 
table(SData$PreOrPost)
table(SData$Response.Comprehensive)
table(SData$CancerType)
> dim(SData)
#[1] 25279 74085
table(SData$Minor_Cell_Type)



GSEA_Result <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/KEGG_GSEA_Result.rds")
Markers_Result <- readRDS(file.path(data_path, "1FindAllMarkers_Result.RDS"))
Markers_Macro_MS4A6A <- Markers_Result %>%
    filter(cluster == "Macro_MS4A6A") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Macro_MS4A6A

Markers_Macro_APOE <- Markers_Result %>%
    filter(cluster == "Macro_APOE") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Macro_APOE
Markers_Macro_FOSB <- Markers_Result %>%
    filter(cluster == "Macro_FOSB") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Macro_FOSB
Markers_Macro_CCL <- Markers_Result %>%
    filter(cluster == "Macro_CCL") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Macro_CCL
Markers_Macro_MARCO <- Markers_Result %>%
    filter(cluster == "Macro_MARCO") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Macro_MARCO
Markers_Mono_FCN1 <- Markers_Result %>%
    filter(cluster == "Mono_FCN1") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Mono_FCN1
Markers_Mono_TIMP1 <- Markers_Result %>%
    filter(cluster == "Mono_TIMP1") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Mono_TIMP1
Markers_DC_HLA <- Markers_Result %>%
    filter(cluster == "DC_HLA") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_DC_HLA
Markers_DC_LAMP3 <- Markers_Result %>%
    filter(cluster == "DC_LAMP3") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_DC_LAMP3
Markers_DC_CPVL <- Markers_Result %>%
    filter(cluster == "DC_CPVL") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_DC_CPVL
Markers_Mast <- Markers_Result %>%
    filter(cluster == "Mast") %>%
    filter(abs(avg_log2FC) >= 1, abs(pct.1-pct.2) >= 0.15) %>%
    arrange(desc(avg_log2FC)) %>%
    select(avg_log2FC, pct.1, pct.2, gene) 
Markers_Mast
