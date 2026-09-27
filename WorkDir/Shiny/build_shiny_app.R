# ============================================================
# ShinyCell2 构建脚本 — Macrophage 项目
# 关键修复: Seurat v5 Assay5 格式需要先 JoinLayers()
# 精简: 只保留5个关键metadata列
# ============================================================
suppressMessages(library(Seurat))
suppressMessages(library(ShinyCell2))

# 读取 Seurat 对象（重命名版：Macro_MS4A6A, Mono_FCN1 等）
#seu <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/data/1PrepSCT_SData_Rename.RDS")
#seu <- NormalizeData(seu, assay = "RNA")
#saveRDS(seu, file = "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/seu_for_shiny.RDS")
seu <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/seu_for_shiny.RDS")

cat("原始 Seurat 对象:\n")
cat("  细胞数:", ncol(seu), "\n")
cat("  基因数:", nrow(seu), "\n")
cat("  Assays:", paste(names(seu@assays), collapse = ", "), "\n")

# ========== 精简 metadata ==========
keep_cols <- c("Minor_Cell_Type", "CancerType",
               "Response.RECIST", "Response.Pathologic", "Response.Comprehensive",
               "PreOrPost", "Patient.ID")
cat("保留的 metadata 列:", paste(keep_cols, collapse = ", "), "\n")

all_cols <- colnames(seu@meta.data)
drop_cols <- setdiff(all_cols, keep_cols)
for (col in drop_cols) {
  seu@meta.data[[col]] <- NULL
}
cat("精简后 metadata 列:", paste(colnames(seu@meta.data), collapse = ", "), "\n")

# ========== Step 1: 创建配置表 ==========
# Patient.ID 有 >50 个水平，需调高 maxLevels 防止被自动排除
scConf <- createConfig(seu, maxLevels = 100)

# ========== Step 2: 调整配置 ==========
scConf <- modMetaName(scConf,
  meta.to.mod = "Minor_Cell_Type",
  new.name = "Myeloid Subtype"
)
scConf <- modMetaName(scConf,
  meta.to.mod = "CancerType",
  new.name = "Cancer Type"
)

scConf <- modDefault(scConf, default1 = "CancerType", default2 = "Minor_Cell_Type")

cat("配置表:\n")
print(scConf[, .(ID, UI, grp, default)])

# ========== Step 3: 导出数据文件 ==========
shiny_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/shinyApp"

# RNA/counts (原始计数) — counts layer 的 dim() 始终正常，无需 JoinLayers
makeShinyFiles(
  seu,
  scConf,
  assay = "RNA",
  assay.slot = "data",
  shiny.prefix = "sc1",
  dimred.to.use = "umap",
  shiny.dir = shiny_dir,
  default.gene1 = "CD68",
  default.gene2 = "CD14",
  default.multigene = c("CD68","CD14","MS4A6A","FCGR3A","APOE","MARCO","CD163",
                         "FOSB","CCL2","TIMP1","FCN1","LAMP3","CPVL","HLA-DRA")
)

cat("数据文件导出完成。\n")

# ========== Step 4: 生成 Shiny 代码 ==========
makeShinyCodes(
  shiny.title = "髓系细胞图谱 — 四种实体瘤中 Macrophage/Monocyte/DC 亚型分析",
  shiny.footnotes = "",
  shiny.prefix = "sc1",
  shiny.dir = shiny_dir,
  defPtSiz = 1.0
)

cat("========================================\n")
cat("构建完成！产出文件在:", shiny_dir, "\n")
cat("========================================\n")
