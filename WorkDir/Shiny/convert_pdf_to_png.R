# ============================================================
# 将三种分析的 PDF 转为 PNG，优化 Shiny App 加载速度
# 在计算节点运行: Rscript convert_pdf_to_png.R
# ============================================================

# 安装 pdftools（如果尚未安装）
if (!require("pdftools")) {
  install.packages("pdftools", repos = "https://cloud.r-project.org")
}

library(pdftools)

base_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir"

# 三种分析的 plot 目录
dirs <- c(
  "Monocle"   = file.path(base_dir, "Monocle/plot"),
  "Cellchat"  = file.path(base_dir, "Cellchat/plot"),
  "hdWGCNA"   = file.path(base_dir, "hdWGCNA/plot")
)

for (analysis in names(dirs)) {
  plot_dir <- dirs[[analysis]]
  png_dir <- file.path(plot_dir, "png")
  dir.create(png_dir, showWarnings = FALSE)
  
  pdf_files <- list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE)
  cat("Converting", length(pdf_files), "PDFs in", analysis, "...\n")
  
  for (pdf_file in pdf_files) {
    png_file <- file.path(png_dir, gsub("\\.pdf$", ".png", basename(pdf_file)))
    
    # 跳过已有 PNG
    if (file.exists(png_file)) {
      cat("  SKIP:", basename(png_file), "\n")
      next
    }
    
    tryCatch({
      # 转换第一页，DPI=150（平衡清晰度和文件大小）
      pdf_convert(pdf_file, format = "png", dpi = 150,
                  filenames = png_file, pages = 1)
      cat("  DONE:", basename(png_file), "\n")
    }, error = function(e) {
      cat("  FAIL:", basename(pdf_file), ":", e$message, "\n")
    })
  }
}

cat("\n========================================\n")
cat("转换完成！PNG 文件位置:\n")
for (analysis in names(dirs)) {
  cat("  ", analysis, ":", file.path(dirs[[analysis]], "png/"), "\n")
}
cat("========================================\n")