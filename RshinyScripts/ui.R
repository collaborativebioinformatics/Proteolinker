#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
options(shiny.maxRequestSize = 2000*1024^2)  
library(shiny)
library(plotly)

fluidPage(
  titlePanel("Olink Data Analysis Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data Upload", style = "color: #2c3e50;"),
      fileInput("prot_file", "Proteomics Parquet File", accept = ".parquet"),
      fileInput("meta_file", "Metadata File (CSV/Parquet)", accept = c(".csv", ".parquet")),
      
      tags$hr(),
      h4("LOD Options", style = "color: #2c3e50;"),
      radioButtons("lod_method", "LOD Method:",
                   choices = c("FixedLOD", "NC_STDEV"),
                   selected = "NC_STDEV"),
      conditionalPanel(
        condition = "input.lod_method == 'FixedLOD'",
        fileInput("fixed_lod_file", "Fixed LOD File (CSV/Parquet)", 
                  accept = c(".csv", ".parquet"))
      ),
      
      tags$hr(),
      checkboxInput("header", "Header (for CSV files)", TRUE),
      tags$hr(),
      downloadButton("download_qc_pdf", "Download QC Report", 
                     class = "btn-primary", style = "width: 100%;"),
      selectInput("ontology", "Select ontology:",
                  choices = c("KEGG", "GO", "Reactome", "All"),
                  selected = "GO"),
      actionButton("run", "Run Analysis"),
      br(),
      downloadButton("download_pathway", "Download Pathway Results (.parquet)",
                     class = "btn-primary", style = "width: 100%;")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Summary", 
                 h4("Sample Counts by Gender and Sample Type"),
                 tableOutput("summary_table")),
        
        tabPanel("Data Completeness", 
                 h4("Olink QC – Data Completeness"),
                 verbatimTextOutput("completeness")),
        
        tabPanel("QC Tables",
                 h4("Sample FAIL"), tableOutput("tbl_sample_fail"),
                 h4("Sample WARN"), tableOutput("tbl_sample_warn"),
                 h4("Assay WARN"), tableOutput("tbl_assay_warn")),
        
        # LOD Analysis moved to the last position
        tabPanel("LOD Analysis",
                 fluidRow(
                   column(6, h4("Overall LOD Statistics")),
                   column(6, h4("Volcano Plot: NPX vs (NPX - LOD)"))
                 ),
                 fluidRow(
                   column(6, tableOutput("lod_stats")),
                   column(6, plotlyOutput("lod_plot", height = "500px"))
                 ),
                 br(),
                 h4("Detailed LOD Summary by Assay"),
                 dataTableOutput("lod_summary")),
        tabPanel("Differential Analysis",
                 h4("Volcano Plot: Gender T-test"),
                 plotOutput("volcanoPlot", height = "600px")),
        
        tabPanel("Pathway Analysis",
                 h4("Pathway Heatmap"),
                 plotOutput("heatmapPlot", height = "500px"),
                 br(),
                 h4("Pathway Bar Chart"),
                 plotOutput("barPlot", height = "500px"))
      )
    )
  )
)