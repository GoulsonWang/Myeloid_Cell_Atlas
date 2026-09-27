#2026/6/5
#本脚本用于查看两个网络中那几个module是有显著差异的.



# Input: 
#   1.hMEs_ESCC.rds

# Output: 
#   1.DMEs_volcano_ESCC.pdf



suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))



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
options(future.globals.maxSize = 3 * 1024^3)



#载入数据
escc_data <- readRDS(file.path(data_path, "hMEs_ESCC.rds"))
escc_data_R <- subset(escc_data, Response.Comprehensive == "Yes")


#设定分组
group1 <- escc_data_R@meta.data %>% subset(PreOrPost == "Pre") %>% rownames
group2 <- escc_data_R@meta.data %>% subset(PreOrPost == "Post") %>% rownames



#检验
DMEs <- FindDMEs(
  escc_data_R,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

head(DMEs)
#                     p_val avg_log2FC pct.1 pct.2     p_val_adj       module
#ESCC-Module1  0.000000e+00  6.6008520 0.492 0.242  0.000000e+00 ESCC-Module1
#ESCC-Module4 6.364939e-251 -2.2366661 0.369 0.641 4.455457e-250 ESCC-Module4
#ESCC-Module7  6.613356e-85 -1.5718462 0.431 0.539  4.629349e-84 ESCC-Module7
#ESCC-Module2  9.063452e-46 -3.8542075 0.457 0.352  6.344416e-45 ESCC-Module2
#ESCC-Module3  8.680202e-15 -0.1222018 0.449 0.500  6.076142e-14 ESCC-Module3
#ESCC-Module6  2.584525e-02 -0.7922006 0.459 0.438  1.809167e-01 ESCC-Module6
#总结，相较于治疗前，治疗后Module1表现为下调（6.6008520），Module2主要表现为上调（-3.8542075）

#可视化
p <- PlotDMEsVolcano(
  escc_data_R,
  DMEs, 
  label_size = 3
)
print(p)
ggsave(p, filename = file.path(plot_path, "DMEs_volcano_ESCC_R.pdf"), width = 5, height = 3)



# 对NR组进行同样的分析
escc_data_NR <- subset(escc_data, Response.Comprehensive == "No")

group1 <- escc_data_NR@meta.data %>% subset(PreOrPost == "Pre") %>% rownames
group2 <- escc_data_NR@meta.data %>% subset(PreOrPost == "Post") %>% rownames

DMEs_NR <- FindDMEs(
  escc_data_NR,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

head(DMEs_NR)
#                     p_val avg_log2FC pct.1 pct.2     p_val_adj       module
#ESCC-Module4 4.818379e-224 -5.7138927 0.266 0.567 3.372865e-223 ESCC-Module4
#ESCC-Module1 1.004975e-132  8.2129861 0.579 0.364 7.034827e-132 ESCC-Module1
#ESCC-Module3  1.641463e-45  0.2220518 0.525 0.374  1.149024e-44 ESCC-Module3
#ESCC-Module2  6.171082e-31 -1.0969774 0.363 0.466  4.319757e-30 ESCC-Module2
#ESCC-Module6  3.624980e-05  0.5123573 0.475 0.434  2.537486e-04 ESCC-Module6
#ESCC-Module7  1.056861e-02 -1.1256589 0.470 0.471  7.398026e-02 ESCC-Module7
#总结，在治疗后，Module1表现下调（8.2129861），Module2表现上升（-1.0969774）


#可视化
p <- PlotDMEsVolcano(
  escc_data_NR,
  DMEs_NR, 
  label_size = 3
)
print(p)
ggsave(p, filename = file.path(plot_path, "DMEs_volcano_ESCC_NR.pdf"), width = 5, height = 3)


# 提取R组和NR组中Module1和Module2的avg_log2FC值
# R组
R_module1_fc <- DMEs[DMEs$module == "ESCC-Module1", "avg_log2FC"]
R_module2_fc <- DMEs[DMEs$module == "ESCC-Module2", "avg_log2FC"]

# NR组
NR_module1_fc <- DMEs_NR[DMEs_NR$module == "ESCC-Module1", "avg_log2FC"]
NR_module2_fc <- DMEs_NR[DMEs_NR$module == "ESCC-Module2", "avg_log2FC"]

# 创建合并数据框，直接处理NR组中Module2的值
combined_plot_data <- rbind(
  data.frame(Group = "R", DMEs[DMEs$module %in% c("ESCC-Module1", "ESCC-Module2"),]),
  data.frame(Group = "NR", DMEs_NR[DMEs_NR$module %in% c("ESCC-Module1", "ESCC-Module2"),])
)

# 为了更好地可视化，我们创建模块名称的简化版本
combined_plot_data$SimpleModule <- gsub("ESCC-", "", combined_plot_data$module)

# 设置因子顺序以确保正确的图例顺序
combined_plot_data$Group <- factor(combined_plot_data$Group, levels = c("R", "NR"))

p_combined <- ggplot(combined_plot_data, aes(x=SimpleModule, y=avg_log2FC, fill=Group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, color = "black", size = 0.3) +
  geom_text(aes(label = round(avg_log2FC, 3), group = Group), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "Differential Module Expression: Pre vs Post Treatment Comparison",
    x = "Module",
    y = "Average log2 Fold Change",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("R" = "#F8766D", "NR" = "#00BFC4"), 
                    labels = c("Responders (R)", "Non-responders (NR)"),
                    guide = guide_legend(order = 1)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylim(c(min(combined_plot_data$avg_log2FC) - 0.5, max(combined_plot_data$avg_log2FC) + 1.5))

# 保存专业版组合图
ggsave(p_combined, filename = file.path(plot_path, "DME_Comparison_ESCC.pdf"), 
       width = 6, height = 7, dpi = 300)

# 输出数值信息
cat("R组中Module1的avg_log2FC:", R_module1_fc, "\n")
cat("R组中Module2的avg_log2FC:", R_module2_fc, "\n")
cat("NR组中Module1的avg_log2FC:", NR_module1_fc, "\n")
cat("NR组中Module2的avg_log2FC:", NR_module2_fc, "\n")