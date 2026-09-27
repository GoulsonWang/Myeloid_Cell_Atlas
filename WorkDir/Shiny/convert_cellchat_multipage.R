# 将多页 CellChat PDF 按页转换为 PNG
library(pdftools)

plot_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Cellchat/plot"
png_dir <- file.path(plot_dir, "png")

files <- c("Aggregate_Network_Count.pdf", "Aggregate_Network_Count2.pdf")

for (f in files) {
  pdf_path <- file.path(plot_dir, f)
  cat("Processing:", f, "\n")
  npages <- pdf_info(pdf_path)$pages
  cat("  Pages:", npages, "\n")
  for (p in 1:npages) {
    out_name <- paste0(gsub("\\.pdf$", "", f), "_p", p, ".png")
    out_path <- file.path(png_dir, out_name)
    if (file.exists(out_path)) {
      cat("  SKIP:", out_name, "\n")
      next
    }
    pdf_convert(pdf_path, format = "png", dpi = 150, filenames = out_path, pages = p)
    cat("  DONE:", out_name, "\n")
  }
}
cat("Done.\n")