library(shiny)
library(OlinkAnalyze)
library(dplyr)
library(arrow)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Olink NPX Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("npx_file", "Upload Olink NPX parquet file", accept = ".parquet"),
      fileInput("meta_file", "Upload metadata CSV (SexManifest)", accept = ".csv"),
      selectInput("ontology", "Select ontology:",
                  choices = c("KEGG", "GO", "Reactome", "All"),
                  selected = "GO"),
      actionButton("run", "Run Analysis"),
      br(),
      downloadButton("download_pathway", "Download Pathway Results (.parquet)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
        tabPanel("Pathway Heatmap", plotOutput("heatmapPlot")),
        tabPanel("Pathway Bar Chart", plotOutput("barPlot"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive NPX dataframe
  npx_data <- reactive({
    req(input$npx_file, input$meta_file)
    olink <- read_parquet(input$npx_file$datapath)
    meta <- read.csv(input$meta_file$datapath, header = TRUE, sep = ",")
    
    olink |>
      filter(!grepl("CONTROL", SampleType)) |>
      filter(!grepl("ctrl", AssayType)) |>
      filter(!is.na(NPX)) |>
      left_join(meta, by = c("SampleID"))
  })
  
  # Reactive T-test results
  ttest_results <- eventReactive(input$run, {
    olink_ttest(
      df = npx_data(),
      variable = "Gender",
      alternative = "two.sided"
    )
  })
  
  # Reactive pathway enrichment
  pathway_results <- eventReactive(input$run, {
    ont <- input$ontology
    if (ont == "All") {
      olink_pathway_enrichment(data = npx_data(), test_results = ttest_results())
    } else {
      olink_pathway_enrichment(data = npx_data(),
                               test_results = ttest_results(),
                               ontology = ont)
    }
  })
  
  # Volcano plot
  output$volcanoPlot <- renderPlot({
    req(ttest_results())
    as.data.frame(ttest_results()) |>
      ggplot(aes(x = estimate, y = -log10(Adjusted_pval), col = Threshold)) +
      geom_point() +
      set_plot_theme()
  })
  
  # Pathway heatmap
  output$heatmapPlot <- renderPlot({
    req(pathway_results())
    olink_pathway_heatmap(enrich_results = pathway_results(),
                          test_results = ttest_results())
  })
  
  # Pathway bar chart
  output$barPlot <- renderPlot({
    req(pathway_results())
    olink_pathway_visualization(enrich_results = pathway_results(),
                                method = "GSEA",
                                number_of_terms = 15)
  })
  
  # Download pathway results
  output$download_pathway <- downloadHandler(
    filename = function() {
      paste0("Pathway_results_", Sys.Date(), ".parquet")
    },
    content = function(file) {
      write_parquet(pathway_results(), file)
    }
  )
}

shinyApp(ui, server)
