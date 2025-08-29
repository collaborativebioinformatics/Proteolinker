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
library(org.Hs.eg.db) 
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
      read.csv(input$meta_file$datapath, header = TRUE, check.names = FALSE)
    } else {
      arrow::read_parquet(input$meta_file$datapath)
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
  
  # --------- Summary: unique SampleIDs by Category x SampleType ----------
  output$summary_table <- renderTable({
    df <- merged_data()
    validate(
      need("Category" %in% names(df), "Merged data must contain 'Category'"),
      need("SampleType" %in% names(df), "Merged data must contain 'SampleType'"),
      need("SampleID" %in% names(df), "Merged data must contain 'SampleID'")
    )
    
    all_types <- c("SAMPLE")
    
    unique_df <- df %>%
      distinct(SampleID, Category, SampleType) %>%
      filter(!is.na(Category)) %>%  
      mutate(SampleType = factor(SampleType, levels = all_types))
    
    counts <- unique_df %>%
      count(Category, SampleType, name = "UniqueCount") %>%
      tidyr::complete(Category, SampleType = all_types, fill = list(UniqueCount = 0))
    
    wide <- counts %>%
      tidyr::pivot_wider(names_from = SampleType, values_from = UniqueCount, values_fill = 0) %>%
      arrange(Category)
    
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
  
  # --------- Fixed LOD Data (from your local file) ----------
  fixed_lod_data <- reactive({
    # Read your existing fixed LOD file
    fixed_lod <- read.csv("../resources/fixed_lod_proper_format.csv", check.names = FALSE)
    
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
      dplyr::select(Assay, any_of(all_categories), Total) %>%
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
      Metric = c("Total Samples/Assays", "Above LOD", "Below LOD", "Missing", "% Below LOD"),
      Value = c(total_samples, above_lod_count, below_lod_count, missing_count, 
                round(below_lod_count/total_samples * 100, 2))
    )
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
  
  de_results <- eventReactive(input$run, {
    if (input$de_mode == "ttest") {
      # Run internal t-test
      olink_ttest(
        df = merged_data(),
        variable = "Category",
        alternative = "two.sided"
      ) 
    } else {
      # Load user-uploaded DE results
      req(input$de_file)
      ext <- tools::file_ext(input$de_file$name)
      df <- if (ext == "csv") {
        read.csv(input$de_file$datapath, stringsAsFactors = FALSE)
      } else if (ext == "parquet") {
        arrow::read_parquet(input$de_file$datapath)
      } else {
        validate("Unsupported file format. Please upload CSV or Parquet.")
      }
      # Validate required columns
      required_cols <- c("Protein", "estimate", "Adjusted_pval")
      missing <- setdiff(required_cols, colnames(df))
      validate(
        need(length(missing) == 0,
             paste("Missing required columns in DE results:", paste(missing, collapse = ", ")))
      )
      # Add Threshold if not present
      if (!"Threshold" %in% colnames(df)) {
        df$Threshold <- ifelse(df$Adjusted_pval < 0.05 & df$estimate > 0, "Up",
                               ifelse(df$Adjusted_pval < 0.05 & df$estimate < 0, "Down", "NS"))
      }
      df
    }
  })
  
  df_to_label <- reactive({
    res <- req(de_results())
    validate(need(!is.null(res), "No results from t-test"))
    
    res_df <- as.data.frame(res)
    dplyr::slice_head(res_df, n = 10)
  })

  
  # --------- Volcano plot ----------
  output$volcanoPlot <- renderPlot({
    df <- req(de_results())
    
    df <- as.data.frame(df)
    validate(need(all(c("estimate", "Adjusted_pval", "Threshold") %in% names(df)),
                  "Volcano plot needs columns: estimate, Adjusted_pval, Threshold"))
    
    labs <- df_to_label()   # <- call the reactive
    if (!is.null(labs)) labs <- as.data.frame(labs)
    
    ggplot(df, aes(x = estimate, y = -log10(Adjusted_pval), colour = Threshold)) +
      geom_point() +
      { if (!is.null(labs) && nrow(labs) > 0)
        geom_text_repel(
          data = labs,
          aes(label = Assay),
          box.padding = 0.5, point.padding = 0.5,
          max.overlaps = Inf,
          segment.color = "grey50"
        )
        else NULL } +
      
      set_plot_theme()
  })
  
  # --------- Pathway Enrichment ----------
  pathway_results <- eventReactive(input$run, {
    ont <- input$ontology
    if (ont == "All") {
      olink_pathway_enrichment(
        data = merged_data(),
        test_results = de_results()
      )
    } else {
      olink_pathway_enrichment(
        data = merged_data(),
        test_results = de_results(),
        ontology = ont
      )
    }
  })
  
  # --------- Pathway heatmap ----------
  output$heatmapPlot <- renderPlot({
    req(pathway_results())
    olink_pathway_heatmap(enrich_results = pathway_results(),
                          test_results = de_results())+ theme(axis.text.x = element_blank())
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