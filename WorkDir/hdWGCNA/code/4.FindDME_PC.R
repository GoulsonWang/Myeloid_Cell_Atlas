#2026/6/5
#本脚本用于查看两个网络中那几个module是有显著差异的.



# Input: 
#   1.hMEs_PC.rds

# Output: 
#   1.DMEs_volcano_PC.pdf



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
pc_data <- readRDS(file.path(data_path, "hMEs_PC.rds"))
pc_data_R <- subset(pc_data, Response.Comprehensive == "Yes")


#设定分组
group1 <- pc_data_R@meta.data %>% subset(PreOrPost == "Pre") %>% rownames
group2 <- pc_data_R@meta.data %>% subset(PreOrPost == "Post") %>% rownames



#检验
DMEs <- FindDMEs(
  pc_data_R,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

head(DMEs)
#                   p_val  avg_log2FC pct.1 pct.2     p_val_adj     module
#PC-Module2 1.402996e-226  -8.0144158 0.137 0.569 9.820971e-226 PC-Module2
#PC-Module5 4.116997e-154  -8.3018149 0.102 0.411 2.881898e-153 PC-Module5
#PC-Module7 1.017447e-143 -20.7250269 0.509 0.213 7.122131e-143 PC-Module7
#PC-Module6 5.537948e-125   6.6290763 0.761 0.437 3.876564e-124 PC-Module6
#PC-Module1  3.042356e-80   3.1711541 0.537 0.312  2.129649e-79 PC-Module1
#PC-Module4  5.231486e-05   0.1988701 0.395 0.354  3.662040e-04 PC-Module4
#总结，相较于治疗前，治疗后Module6表现为下调（6.6290763），Module7主要表现为上调（20.7250269）

#可视化
p <- PlotDMEsVolcano(
  pc_data_R,
  DMEs, 
  label_size = 3, 
  xlim_range = c(-21, 8)
)
print(p)
ggsave(p, filename = file.path(plot_path, "DMEs_volcano_PC_R.pdf"), width = 5, height = 3)



# 对NR组进行同样的分析
pc_data_NR <- subset(pc_data, Response.Comprehensive == "No")

group1 <- pc_data_NR@meta.data %>% subset(PreOrPost == "Pre") %>% rownames
group2 <- pc_data_NR@meta.data %>% subset(PreOrPost == "Post") %>% rownames

DMEs_NR <- FindDMEs(
  pc_data_NR,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

DMEs_NR
#                  p_val avg_log2FC pct.1 pct.2    p_val_adj     module
#PC-Module6 3.936084e-21  2.2796451 0.555 0.360 2.755259e-20 PC-Module6
#PC-Module7 6.888698e-01  3.7976238 0.198 0.211 1.000000e+00 PC-Module7
#总结，Module6在治疗后下调2.2796451，Module7无明显变化


#可视化
p <- PlotDMEsVolcano(
  pc_data_NR,
  DMEs_NR, 
  label_size = 3
)
print(p)
ggsave(p, filename = file.path(plot_path, "DMEs_volcano_PC_NR.pdf"), width = 5, height = 3)


# 创建合并数据框，直接处理NR组中Module7的值
combined_plot_data <- rbind(
  data.frame(Group = "R", DMEs[DMEs$module %in% c("PC-Module6", "PC-Module7"),]),
  data.frame(Group = "NR", 
             transform(DMEs_NR[DMEs_NR$module %in% c("PC-Module6", "PC-Module7"),], 
                       avg_log2FC = ifelse(module == "PC-Module7", 0, avg_log2FC)))
)

# 为了更好地可视化，我们创建模块名称的简化版本
combined_plot_data$SimpleModule <- gsub("PC-", "", combined_plot_data$module)

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
ggsave(p_combined, filename = file.path(plot_path, "DME_Comparison_PC.pdf"), 
       width = 6, height = 7, dpi = 300)

# 输出数值信息
cat("R组中Module6的avg_log2FC:", R_module6_fc, "\n")
cat("R组中Module7的avg_log2FC:", R_module7_fc, "\n")
cat("NR组中Module6的avg_log2FC:", NR_module6_fc, "\n")
cat("NR组中Module7的avg_log2FC:", NR_module7_fc, "\n")
