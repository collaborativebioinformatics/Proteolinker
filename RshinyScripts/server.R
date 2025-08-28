#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(tidyr)
library(OlinkAnalyze)
library(gridExtra)    
library(grid)  

function(input, output, session) {
  
  # --------- Loaders ----------
  prot_data <- reactive({
    req(input$prot_file)
    read_parquet(input$prot_file$datapath)
  })
  
  meta_data <- reactive({
    req(input$meta_file)
    if (grepl("\\.csv$", input$meta_file$name, ignore.case = TRUE)) {
      read.csv(input$meta_file$datapath, header = input$header, check.names = FALSE)
    } else {
      read_parquet(input$meta_file$datapath)
    }
  })
  
  # --------- Merged data for summary (keep all proteomics samples) ----------
  merged_data <- reactive({
    df1 <- prot_data()
    df2 <- meta_data()
    validate(
      need("SampleID" %in% names(df1), "Proteomics file must have 'SampleID'"),
      need("SampleID" %in% names(df2), "Metadata file must have 'SampleID'")
    )
    left_join(df1, df2, by = "SampleID")
  })
  
  # --------- Summary: unique SampleIDs by Gender x SampleType ----------
  output$summary_table <- renderTable({
    df <- merged_data()
    validate(
      need("Gender" %in% names(df), "Merged data must contain 'Gender'"),
      need("SampleType" %in% names(df), "Merged data must contain 'SampleType'"),
      need("SampleID" %in% names(df), "Merged data must contain 'SampleID'")
    )
    
    all_types <- c("PLATE_CONTROL", "NEGATIVE_CONTROL", "SAMPLE_CONTROL", "SAMPLE")
    
    unique_df <- df %>%
      distinct(SampleID, Gender, SampleType) %>%
      mutate(
        Gender = ifelse(is.na(Gender) | Gender == "", "Unknown", as.character(Gender)),
        SampleType = factor(SampleType, levels = all_types)
      )
    
    counts <- unique_df %>%
      count(Gender, SampleType, name = "UniqueCount") %>%
      tidyr::complete(Gender, SampleType = all_types, fill = list(UniqueCount = 0))
    
    wide <- counts %>%
      tidyr::pivot_wider(names_from = SampleType, values_from = UniqueCount, values_fill = 0) %>%
      arrange(Gender)
    
    with_totals <- wide %>%
      mutate(across(-Gender, as.integer)) %>%
      mutate(Total = as.integer(rowSums(across(where(is.integer))))) %>%
      {
        totals <- summarise(., across(where(is.integer), sum)) %>%
          mutate(Gender = "Total") %>%
          select(Gender, everything())
        bind_rows(., totals)
      }
    
    with_totals
  })
  
  # --------- Olink QC: completeness + WARN/FAIL tables ----------
  qc_list <- reactive({
    df <- prot_data()
    # check_data_completeness is not always exported; keep ::: as in your script
    compl <- tryCatch(
      OlinkAnalyze:::check_data_completeness(df),
      error = function(e) paste("check_data_completeness() error:", e$message)
    )
    
    # Safe guards if columns are missing
    need_cols <- c("AssayQC","SampleQC","PlateID","Assay","Block","OlinkID","SampleID")
    missing <- setdiff(need_cols, names(df))
    validate(need(length(missing) == 0,
                  paste("Proteomics file missing columns:", paste(missing, collapse = ", "))))
    
    AssayWarn <- df %>% filter(grepl("WARN", AssayQC, ignore.case = TRUE))
    SampleWarn <- df %>% filter(grepl("WARN", SampleQC, ignore.case = TRUE))
    SampleFail <- df %>% filter(grepl("FAIL", SampleQC, ignore.case = TRUE))
    
    AssayWarnData <- AssayWarn %>%
      distinct(PlateID, Assay, Block, OlinkID) %>%
      group_by(Assay, Block, OlinkID) %>%
      summarise(across(everything(), ~toString(.)), .groups = "drop")
    
    SampleWarnData <- SampleWarn %>%
      distinct(SampleID, Block, PlateID) %>%
      group_by(SampleID, PlateID) %>%
      summarise(across(everything(), ~toString(.)), .groups = "drop")
    
    SampleFailData <- SampleFail %>%
      distinct(SampleID, Block, PlateID) %>%
      group_by(SampleID, PlateID) %>%
      summarise(across(everything(), ~toString(.)), .groups = "drop")
    
    list(
      completeness = compl,
      assay_warn = AssayWarnData,
      sample_warn = SampleWarnData,
      sample_fail = SampleFailData
    )
  })
  
  output$completeness <- renderPrint({
    qc <- qc_list()$completeness
    qc
  })
  
  output$tbl_sample_fail <- renderTable({
    dat <- qc_list()$sample_fail
    if (nrow(dat) == 0) return(data.frame(Message = "No FAILED samples"))
    head(dat, 50)
  })
  
  output$tbl_sample_warn <- renderTable({
    dat <- qc_list()$sample_warn
    if (nrow(dat) == 0) return(data.frame(Message = "No WARNED samples"))
    head(dat, 50)
  })
  
  output$tbl_assay_warn <- renderTable({
    dat <- qc_list()$assay_warn
    if (nrow(dat) == 0) return(data.frame(Message = "No WARNED assays"))
    head(dat, 50)
  })
  
  # --------- Download QC PDF ----------
  output$download_qc_pdf <- downloadHandler(
    filename = function() {
      paste0("Olink_QC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      qc <- qc_list()
      pdf(file = file, width = 11, height = 8, onefile = TRUE)
      
      # Page 1: Sample FAIL
      grid.arrange(
        textGrob("Sample Failed", gp = gpar(fontsize = 16, fontface = "bold")),
        tableGrob(if (nrow(qc$sample_fail) == 0) data.frame(Message = "No FAILED samples") else qc$sample_fail),
        ncol = 1, heights = c(0.12, 0.88)   # grid.arrange opens a new page automatically
      )
      
      # Page 2: Sample WARN
      grid.arrange(
        textGrob("Sample Warned", gp = gpar(fontsize = 16, fontface = "bold")),
        tableGrob(if (nrow(qc$sample_warn) == 0) data.frame(Message = "No WARNED samples") else qc$sample_warn),
        ncol = 1, heights = c(0.12, 0.88)
      )
      
      # Page 3: Assay WARN
      grid.arrange(
        textGrob("Assay Warned", gp = gpar(fontsize = 16, fontface = "bold")),
        tableGrob(if (nrow(qc$assay_warn) == 0) data.frame(Message = "No WARNED assays") else qc$assay_warn),
        ncol = 1, heights = c(0.12, 0.88)
      )
      
      dev.off()
    }
  )
  
}
