# 将 paper_writing/plot/ 下的 PDF 转为 PNG
library(pdftools)

plot_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/plot"
png_dir <- file.path(plot_dir, "png")
dir.create(png_dir, showWarnings = FALSE)

pdf_files <- list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE)
cat("Converting", length(pdf_files), "PDFs...\n")

for (pdf_file in pdf_files) {
  png_file <- file.path(png_dir, gsub("\\.pdf$", ".png", basename(pdf_file)))
  if (file.exists(png_file)) {
    cat("  SKIP:", basename(png_file), "\n")
    next
  }
  tryCatch({
    pdf_convert(pdf_file, format = "png", dpi = 150, filenames = png_file, pages = 1)
    cat("  DONE:", basename(png_file), "\n")
  }, error = function(e) {
    cat("  FAIL:", basename(pdf_file), ":", e$message, "\n")
  })
}
cat("Done.\n")