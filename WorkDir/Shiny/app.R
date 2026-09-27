# ============================================================
# Macrophage 项目结果展示 — 简化版 Shiny App
# 展示内容: 拟时序分析 / 细胞通讯 / 共表达网络
library(shiny)
library(bslib)

# 将分析结果的 plot 目录暴露给 Shiny（通过HTTP 而不是 file:// 协议）
work_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir"
addResourcePath("plots", file.path(work_dir, "Monocle/plot"))
addResourcePath("plots_cc", file.path(work_dir, "Cellchat/plot"))
addResourcePath("plots_hd", file.path(work_dir, "hdWGCNA/plot"))

# 定义各分析的图片文件名和标签（只需要文件名，路径由 addResourcePath 处理）
monocle_plots <- list(
  "分化轨迹(全图)"     = "trajectory_setPresudoTime0.pdf",
  "分化轨迹(分支0)"    = "trajectory_branch0.pdf",
  "分化轨迹(分支1_1)"  = "trajectory_branch1_1.pdf",
  "分化轨迹(分支1_2)"  = "trajectory_branch1_2.pdf",
  "通路得分-分支2"     = "pathway_branch_score2.pdf",
  "通路得分-分支3"     = "pathway_branch_score3.pdf",
  "通路得分-分支4"     = "pathway_branch_score4.pdf",
  "通路得分-分支5"     = "pathway_branch_score5.pdf",
  "通路得分-分支6"     = "pathway_branch_score6.pdf",
  "通路得分-分支8"     = "pathway_branch_score8.pdf",
  "通路得分-分支9"     = "pathway_branch_score9.pdf",
  "KEGG-GSEA(分支2)"   = "KEGG_GSEAbranch2_1.pdf",
  "UMAP-Seurat聚类"    = "umap_preparing_By_Seurat.pdf",
  "UMAP-Monocle分区"   = "umap_preparing_By_monocle.pdf"
)

cellchat_plots <- list(
  "总网络图"           = "Aggregate_Network.pdf",
  "总网络图(计数)"     = "Aggregate_Network_Count.pdf",
  "总网络图(计数2)"    = "Aggregate_Network_Count2.pdf",
  "总网络图(计数3)"    = "Aggregate_Network_Count3.pdf",
  "全部通路网络"       = "Pathway_Network_all.pdf",
  "亚型权重网络"       = "Subtype_weight_Networks.pdf",
  "亚型权重网络2"      = "Subtype_weight_Networks2.pdf",
  "CellChatDB"         = "CellChatDB.pdf",
  "通路1网络"          = "pathway_1network.pdf",
  "通路2网络"          = "pathway_2network.pdf",
  "通路3网络"          = "pathway_3network.pdf",
  "通路4网络"          = "pathway_4network.pdf",
  "通路5网络"          = "pathway_5network.pdf",
  "通路6网络"          = "pathway_6network.pdf",
  "通路7网络"          = "pathway_7network.pdf",
  "通路8网络"          = "pathway_8network.pdf",
  "通路9网络"          = "pathway_9network.pdf",
  "通路10网络"         = "pathway_10network.pdf",
  "通路11网络"         = "pathway_11network.pdf"
)

hdwgcna_plots <- list(
  "模块树状图(ESCC)"    = "Dendrogram_ESCC.pdf",
  "UMAP-MetaCell"       = "UMAP_MetaCell.pdf",
  "UMAP-MetaCell(ESCC)" = "UMAP_MetaCell_ESCC.pdf",
  "UMAP-细胞类型"       = "UMAP_Minor_Cell_Type.pdf",
  "UMAP-细胞类型(ESCC)" = "UMAP_Minor_Cell_Type_ESCC.pdf",
  "模块表达(ESCC)"      = "hMEs_per_module_ESCC.pdf",
  "模块表达(PC)"        = "hMEs_per_module_PC.pdf",
  "Hub基因得分(ESCC)"   = "hub_scores_ESCC.pdf",
  "Hub基因得分(PC)"     = "hub_scores_PC.pdf",
  "Hub基因合并(PC)"     = "1.hub_gene_combined_PC.pdf",
  "雷达图(ESCC)"        = "radar_ESCC.pdf",
  "雷达图(PC)"          = "radar_PC.pdf",
  "Dotplot(ESCC)"       = "dotplot_ESCC.pdf",
  "Dotplot(PC)"         = "dotplot_PC.pdf",
  "DME火山图(PC)"       = "DMEs_volcano_PC_NR.pdf",
  "DME比较(PC)"         = "DME_Comparison_PC.pdf",
  "kMEs(ESCC)"          = "kMEs_ESCC.pdf",
  "T细胞激活通路(PC)"   = "2.PC_six_layouts_combined_positive_regulation_of_T_cell_activation.pdf",
  "T细胞激活通路(ESCC)" = "2.ESCC_six_layouts_combined_positive_regulation_of_T_cell_activation.pdf",
  "分组柱状图(PC)"      = "4.Grouped_Barplot_PC.pdf",
  "软阈值检验(ESCC)"    = "soft_power_test_ESCC.pdf",
  "软阈值检验(PC)"      = "soft_power_test_PC.pdf"
)

# ========================== UI ==========================
ui <- page_navbar(
  title = "髓系细胞单细胞分析结果展示",

  # Tab 1: 拟时序分析
  nav_panel(
    "拟时序分析",
    page_sidebar(
      sidebar = sidebar(
        selectInput("monocle_plot", "选择图片:",
          choices = names(monocle_plots),
          selected = names(monocle_plots)[1]
        ),
        helpText("Monocle3 拟时序分析结果。"),
        helpText("展示髓系细胞从 Monocyte 向 Macrophage/DC 分化的轨迹。")
      ),
      card(
        full_screen = TRUE,
        card_header("拟时序分析"),
        uiOutput("monocle_viewer")
      )
    )
  ),

  # Tab 2: 细胞通讯
  nav_panel(
    "细胞通讯",
    page_sidebar(
      sidebar = sidebar(
        selectInput("cellchat_plot", "选择图片:",
          choices = names(cellchat_plots),
          selected = names(cellchat_plots)[1]
        ),
        helpText("CellChat 配体-受体通讯分析。"),
        helpText("展示髓系细胞与T细胞之间的信号交流网络。")
      ),
      card(
        full_screen = TRUE,
        card_header("细胞通讯分析"),
        uiOutput("cellchat_viewer")
      )
    )
  ),

  # Tab 3: 共表达网络
  nav_panel(
    "共表达网络",
    page_sidebar(
      sidebar = sidebar(
        selectInput("hdwgcna_plot", "选择图片:",
          choices = names(hdwgcna_plots),
          selected = names(hdwgcna_plots)[1]
        ),
        helpText("hdWGCNA 基因共表达网络分析。"),
        helpText("按ESCC和PC分别构建，找出驱动肿瘤免疫的基因模块。")
      ),
      card(
        full_screen = TRUE,
        card_header("共表达网络分析"),
        uiOutput("hdwgcna_viewer")
      )
    )
  )
)

# ========================== Server ==========================
server <- function(input, output, session) {

  # 拟时序图片渲染
  output$monocle_viewer <- renderUI({
    tags$iframe(
      src = paste0("plots/", monocle_plots[[input$monocle_plot]]),
      width = "100%",
      height = "800px",
      style = "border: none;"
    )
  })

  # CellChat 图片渲染
  output$cellchat_viewer <- renderUI({
    tags$iframe(
      src = paste0("plots_cc/", cellchat_plots[[input$cellchat_plot]]),
      width = "100%",
      height = "800px",
      style = "border: none;"
    )
  })

  # hdWGCNA 图片渲染
  output$hdwgcna_viewer <- renderUI({
    tags$iframe(
      src = paste0("plots_hd/", hdwgcna_plots[[input$hdwgcna_plot]]),
      width = "100%",
      height = "800px",
      style = "border: none;"
    )
  })
}

# ======================== 启动 =========================
shinyApp(ui, server)