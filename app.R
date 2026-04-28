library(shiny)
library(DESeq2)
library(airway)
library(EnhancedVolcano)
library(pheatmap)
library(DT)

# --- Prepare data once at startup ---
data(airway)
dds <- DESeqDataSet(airway, design = ~ dex)
dds$dex <- relevel(dds$dex, ref = "untrt")
dds <- DESeq(dds)
res <- results(dds)
vsd <- vst(dds, blind = FALSE)

# Pre-save heatmap to a fixed path
heatmap_path <- file.path(getwd(), "heatmap.png")

# --- UI ---
ui <- fluidPage(
  titlePanel("DESeq2 Explorer - Airway Dataset"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("pval", "P-value cutoff:",
                  min=0.001, max=0.05, value=0.05, step=0.001),
      sliderInput("lfc", "Log2 Fold Change cutoff:",
                  min=0.5, max=3, value=1.5, step=0.1),
      sliderInput("topn", "Top N genes for heatmap:",
                  min=10, max=100, value=50, step=10)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot", plotOutput("volcano", height="500px")),
        tabPanel("Heatmap",     plotOutput("heatmap",  height="600px")),
        tabPanel("Results Table", DT::dataTableOutput("table"))
      )
    )
  )
)

# --- Server ---
server <- function(input, output) {
  
  output$volcano <- renderPlot({
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "pvalue",
                    pCutoff = input$pval,
                    FCcutoff = input$lfc,
                    title = "Dexamethasone vs Untreated",
                    pointSize = 2.0,
                    labSize = 3.0)
  })
  
  output$heatmap <- renderPlot({
    library(ComplexHeatmap)
    
    topN <- head(order(res$padj, na.last = TRUE), input$topn)
    mat  <- assay(vsd)[topN, ]
    
    # Scale rows for better visualization
    mat_scaled <- t(scale(t(mat)))
    
    # Color for treatment annotation
    col_anno <- HeatmapAnnotation(
      treatment = colData(vsd)$dex,
      col = list(treatment = c("untrt" = "#a29bfe", "trt" = "#f9ca24"))
    )
    
    Heatmap(mat_scaled,
            name = "Z-score",
            top_annotation = col_anno,
            show_row_names = TRUE,
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 6),
            column_title = paste("Top", input$topn, "DEGs"),
            clustering_distance_rows = "euclidean",
            clustering_distance_columns = "euclidean")
  })
  
  output$table <- DT::renderDataTable({
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::arrange(padj)
  })
}

# --- Run ---
shinyApp(ui, server)