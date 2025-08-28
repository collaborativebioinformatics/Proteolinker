#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

fluidPage(
  titlePanel("Upload Your Data"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("prot_file", "Choose Proteomics Parquet File", accept = ".parquet"),
      fileInput("meta_file", "Choose Metadata File (CSV/Parquet)", accept = c(".csv", ".parquet")),
      tags$hr(),
      checkboxInput("header", "Header (for CSV files)", TRUE),
      tags$hr(),
      downloadButton("download_qc_pdf", "Download QC PDF")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", 
                 h4("Summary: counts by Sex x SampleType"),
                 tableOutput("summary_table")),
        tabPanel("Data Completeness", 
                 h4("Olink QC – Data Completeness"),
                 verbatimTextOutput("completeness")),
        tabPanel("QC Tables",
                 h4("Sample FAIL"), tableOutput("tbl_sample_fail"),
                 h4("Sample WARN"), tableOutput("tbl_sample_warn"),
                 h4("Assay WARN"), tableOutput("tbl_assay_warn"))
      )
    )
  )
)
