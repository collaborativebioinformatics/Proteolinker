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
library(purrr)
library(DT)
library(plotly)
library(ggplot2)
library(ggrepel)

function(input, output, session) {
  
   # --------- Loaders ----------
  prot_data <- reactive({
    req(input$prot_file)
    arrow::read_parquet(input$prot_file$datapath)
  })
  
  meta_data <- reactive({
    req(input$meta_file)
    if (grepl("\\.csv$", input$meta_file$name, ignore.case = TRUE)) {
      read.csv(input$meta_file$datapath, header = input$header, check.names = FALSE)
    } else {
      arrow::read_parquet(input$meta_file$datapath)
    }
  })
  
  # --------- Fixed LOD Data (from your local file) ----------
  fixed_lod_data <- reactive({
    # Read your existing fixed LOD file
    fixed_lod <- read.csv("fixed_lod_proper_format.csv", check.names = FALSE)
    
    # Ensure the file has the expected structure
    validate(
      need("Assay" %in% names(fixed_lod), "Fixed LOD file must contain 'Assay' column"),
      need(any(c("LOD", "LODNPX") %in% names(fixed_lod)), 
           "Fixed LOD file must contain 'LOD' or 'LODNPX' column")
    )
    
    # Rename LODNPX to LOD if needed
    if ("LODNPX" %in% names(fixed_lod) && !"LOD" %in% names(fixed_lod)) {
      fixed_lod <- fixed_lod %>% rename(LOD = LODNPX)
    }
    
    return(fixed_lod)
  })
  
  # --------- Display Fixed LOD Table ----------
  output$fixed_lod_table <- renderDataTable({
    fixed_lod <- fixed_lod_data()
    datatable(fixed_lod, options = list(
      pageLength = 10,
      lengthMenu = c(5, 10, 15),
      searching = TRUE
    ))
  })
  
  # --------- LOD Processing ----------
  lod_data <- reactive({
    req(prot_data(), input$lod_method)
    
    npx_data <- prot_data()
    
    # Remove AssayType filter if column doesn't exist - add check
    if ("AssayType" %in% names(npx_data)) {
      npx_data <- npx_data %>%
        filter(!grepl("ext_ctrl", AssayType, ignore.case = TRUE)) %>%
        filter(!grepl("inc_ctrl", AssayType, ignore.case = TRUE)) %>%
        filter(!grepl("amp_ctrl", AssayType, ignore.case = TRUE))
    }
    
    if (input$lod_method == "FixedLOD") {
      # Use the fixed LOD data from your file
      fixed_lod <- fixed_lod_data()
      
      # Apply fixed LOD - manual join
      result <- npx_data %>%
        left_join(fixed_lod %>% select(Assay, LOD), by = "Assay") %>%
        mutate(LOD = as.numeric(LOD))
      
    } else if (input$lod_method == "NC_STDEV") {
      # Calculate LOD based on negative controls
      # Check if SampleType column exists
      if ("SampleType" %in% names(npx_data)) {
        negative_controls <- npx_data %>% 
          filter(grepl('NEGATIVE_CONTROL', SampleType, ignore.case = TRUE))
      } else {
        # If no SampleType, try to identify negative controls by SampleID pattern
        negative_controls <- npx_data %>% 
          filter(grepl('NC|NEG', SampleID, ignore.case = TRUE))
      }
      
      if (nrow(negative_controls) == 0) {
        showNotification("No negative controls found for NC_STDEV method", type = "warning")
        return(NULL)
      }
      
      result <- negative_controls %>%
        group_by(PlateID, Assay) %>%
        summarise(
          LOD = mean(NPX, na.rm = TRUE) + 3 * sd(NPX, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        right_join(npx_data, by = c("PlateID", "Assay"))
    }
    
    # Add threshold classification - handle cases where LOD might be missing
    if ("LOD" %in% names(result)) {
      result <- result %>%
        mutate(
          Threshold = case_when(
            NPX < LOD ~ "Below_LOD",
            is.na(NPX) ~ "Missing",
            TRUE ~ "Above_LOD"
          )
        )
    } else {
      result <- result %>%
        mutate(
          Threshold = ifelse(is.na(NPX), "Missing", "Above_LOD")
        )
    }
    
    return(result)
  })
  
  # --------- LOD Summary Table with Pagination ----------
  output$lod_summary <- renderDataTable({
    df <- lod_data()
    
    # Ensure Threshold column exists
    if (!"Threshold" %in% names(df)) {
      return(data.frame(Message = "No threshold data available"))
    }
    
    # Create summary with all possible threshold categories
    all_categories <- c("Above_LOD", "Below_LOD", "Missing")
    
    summary_table <- df %>%
      group_by(Assay) %>%
      count(Threshold) %>%
      ungroup() %>%
      complete(Assay, Threshold = all_categories, fill = list(n = 0)) %>%
      pivot_wider(names_from = Threshold, values_from = n, values_fill = 0) %>%
      mutate(Total = rowSums(across(where(is.numeric)))) %>%
      select(Assay, any_of(all_categories), Total) %>%
      arrange(desc(Below_LOD))  # Sort by most Below_LOD first
    
    datatable(summary_table, options = list(
      pageLength = 25,
      lengthMenu = c(10, 25, 50, 100),
      searching = TRUE,
      ordering = TRUE
    ))
  })
  
  # --------- LOD Visualization with Volcano Plot ----------
  output$lod_plot <- renderPlotly({
    df <- lod_data()
    
    if (!"Threshold" %in% names(df) || !"LOD" %in% names(df)) {
      return(plotly_empty() %>% 
               layout(title = list(text = "No LOD data available for volcano plot", 
                                   x = 0.5, y = 0.5, xanchor = 'center', yanchor = 'center')))
    }
    
    # Filter for sample types you want - adjust the filter as needed
    # Using 'SC' as in your example, but you can modify this
    plot_data <- df %>%
      filter(grepl('SC', SampleID, ignore.case = TRUE))
    
    if (nrow(plot_data) == 0) {
      return(plotly_empty() %>% 
               layout(title = list(text = "No samples matching 'SC' filter criteria", 
                                   x = 0.5, y = 0.5, xanchor = 'center', yanchor = 'center')))
    }
    
    # Create volcano plot
    p <- plot_ly(plot_data, 
                 x = ~(NPX - LOD), 
                 y = ~NPX, 
                 color = ~Threshold,
                 colors = c("Above_LOD" = "#2E8B57", "Below_LOD" = "#DC143C", "Missing" = "#808080"),
                 type = 'scatter',
                 mode = 'markers',
                 hoverinfo = 'text',
                 hovertext = ~paste('SampleID:', SampleID,
                                    '<br>Assay:', Assay,
                                    '<br>NPX:', round(NPX, 3),
                                    '<br>LOD:', round(LOD, 3),
                                    '<br>NPX-LOD:', round(NPX - LOD, 3),
                                    '<br>Status:', Threshold),
                 marker = list(size = 6, opacity = 0.7)) %>%
      layout(
        title = list(text = "<b>Volcano Plot: NPX vs (NPX - LOD)</b>", 
                     x = 0.5, xanchor = 'center'),
        xaxis = list(title = "NPX - LOD", 
                     zeroline = TRUE,
                     zerolinecolor = 'black',
                     zerolinewidth = 2),
        yaxis = list(title = "NPX"),
        legend = list(orientation = 'h', x = 0.5, y = -0.2, xanchor = 'center'),
        margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
        hoverlabel = list(bgcolor = "white", font = list(size = 12))
      ) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE,
             modeBarButtonsToRemove = c('pan2d', 'lasso2d', 'select2d'))
    
    p
  })
  
  # --------- Overall LOD Statistics ----------
  output$lod_stats <- renderTable({
    df <- lod_data()
    
    if (!"Threshold" %in% names(df)) {
      return(data.frame(Message = "No threshold data available"))
    }
    
    # Overall statistics
    total_samples <- nrow(df)
    below_lod_count <- sum(df$Threshold == "Below_LOD", na.rm = TRUE)
    above_lod_count <- sum(df$Threshold == "Above_LOD", na.rm = TRUE)
    missing_count <- sum(df$Threshold == "Missing", na.rm = TRUE)
    
    data.frame(
      Metric = c("Total Samples", "Above LOD", "Below LOD", "Missing", "% Below LOD"),
      Value = c(total_samples, above_lod_count, below_lod_count, missing_count, 
                round(below_lod_count/total_samples * 100, 2))
    )
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
    
    all_types <- c("SAMPLE")
    
    unique_df <- df %>%
      distinct(SampleID, Gender, SampleType) %>%
      filter(Gender %in% c("Male", "Female")) %>%   # keep only Male/Female
      mutate(SampleType = factor(SampleType, levels = all_types))
    
    counts <- unique_df %>%
      count(Gender, SampleType, name = "UniqueCount") %>%
      tidyr::complete(Gender, SampleType = all_types, fill = list(UniqueCount = 0))
    
    wide <- counts %>%
      tidyr::pivot_wider(names_from = SampleType, values_from = UniqueCount, values_fill = 0) %>%
      arrange(Gender)
    
    wide
  })
  
  # --------- Olink QC: completeness + WARN/FAIL tables ----------
  qc_list <- reactive({
    df <- prot_data()
    # check_data-completeness is not always exported; keep ::: as in your script
    compl <- tryCatch(
      {
        OlinkAnalyze:::check_data_completeness(df)  # run QC
        "No issues with Input data"                   # if no error, return this text
      },
      error = function(e) {
        paste("check_data_completeness() error:", e$message)
      }
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
      lod_df <- lod_data()
      
      pdf(file = file, width = 11, height = 8, onefile = TRUE)
      
      # Page 1: Sample FAIL
      grid.arrange(
        textGrob("Sample Failed", gp = gpar(fontsize = 16, fontface = "bold")),
        tableGrob(if (nrow(qc$sample_fail) == 0) data.frame(Message = "No FAILED samples") else qc$sample_fail),
        ncol = 1, heights = c(0.12, 0.88)
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
      
      # Page 4: LOD Summary (only if threshold data exists)
      if ("Threshold" %in% names(lod_df)) {
        lod_summary <- lod_df %>%
          group_by(Assay) %>%
          count(Threshold) %>%
          ungroup() %>%
          complete(Assay, Threshold = c("Above_LOD", "Below_LOD", "Missing"), fill = list(n = 0)) %>%
          pivot_wider(names_from = Threshold, values_from = n, values_fill = 0)
        
        grid.arrange(
          textGrob("LOD Summary", gp = gpar(fontsize = 16, fontface = "bold")),
          tableGrob(lod_summary),
          ncol = 1, heights = c(0.12, 0.88)
        )
      }
      
      dev.off()
    }
  )
  # --------- T-test (reactive on button press) ----------
  ttest_results <- eventReactive(input$run, {
    olink_ttest(
      df = merged_data(),
      variable = "Gender",
      alternative = "two.sided"
    )
  })

df_to_label <- ttest_results%>%
  dplyr::slice_head(n = 10)
        
  # --------- Pathway Enrichment ----------
  pathway_results <- eventReactive(input$run, {
    ont <- input$ontology
    if (ont == "All") {
      olink_pathway_enrichment(data = merged_data(),
                               test_results = ttest_results())
    } else {
      olink_pathway_enrichment(data = merged_data(),
                               test_results = ttest_results(),
                               ontology = ont)
    }
  })
  library(ggrepel)
        
  # --------- Volcano plot ----------
  output$volcanoPlot <- renderPlot({
    req(ttest_results())
    as.data.frame(ttest_results()) %>%
      geom_text_repel(data = df_to_label, aes(label = Assay), 
                  box.padding = 0.5, point.padding = 0.5, 
                  max.overlaps = Inf, # Allows more labels to be shown
                  segment.color = 'grey50')+
  set_plot_theme()
  })
  
  # --------- Pathway heatmap ----------
  output$heatmapPlot <- renderPlot({
    req(pathway_results())
    olink_pathway_heatmap(enrich_results = pathway_results(),
                          test_results = ttest_results())+ theme(axis.text.x = element_blank())
  })
  
  # --------- Pathway bar chart ----------
  output$barPlot <- renderPlot({
    req(pathway_results())
    olink_pathway_visualization(enrich_results = pathway_results(),
                                method = "GSEA",
                                number_of_terms = 15)
  })
  
  # --------- Download pathway results ----------
  output$download_pathway <- downloadHandler(
    filename = function() {
      paste0("Pathway_results_", Sys.Date(), ".parquet")
    },
    content = function(file) {
      arrow::write_parquet(pathway_results(), file)
    }
  )
}
