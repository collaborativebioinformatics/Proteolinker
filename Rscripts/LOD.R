library(tidyverse)      
library(palmerpenguins) 
library(OlinkAnalyze)   
library(arrow)          
library(magrittr)       
library(dplyr)          
library(umap)           
library(ggplot2)        
library(gtable)         
library(grid)           
library(gridExtra)      
library(readxl)         
library(reshape2)     
library(kableExtra)
library(purrr)

LOD_opt <- readline("Enter your LOD optiones (1 for FixedLOD, 2 for NC STDEV LOD): ")
npx_parquet<-Olink_datafile|>
  filter(!grepl("ext_ctrl", AssayType))|>
  filter(!grepl("inc_ctrl", AssayType))|>
  filter(!grepl("amp_ctrl", AssayType))

if (LOD_opt==1){
  fixedLOD_filepath <- "/Users/qiaoyanw/Desktop/ExploreHT_Fixed LOD_2024-12-19.csv"
  LOD_npx_parquet <-olink_lod(npx_parquet, lod_file_path = fixedLOD_filepath, lod_method = "FixedLOD")
  
  LOD_model<-LOD_npx_parquet|>
    dplyr::mutate(Threshold  = ifelse(NPX<LOD,"LOD","Above_LOD"))|>
    mutate(Threshold = replace_na(Threshold, "Above_LOD"))
}

if (LOD_opt==2){
  LOD_Assay <- npx_parquet %>% 
    filter(grepl('NEGATIVE_CONTROL',SampleType, ignore.case = TRUE)) %>%
    group_by(Assay) %>%
    summarise(LOD=median(NPX)+3*sd(NPX),PCNormalizedLOD=median(PCNormalizedNPX)+3*sd(PCNormalizedNPX))
  LOD_npx_parquet_NCStdev<-select(LOD_Assay, Assay, PCNormalizedLOD, LOD) %>%
    left_join(npx_parquet, by=c("Assay"))
  
  LOD_model<-LOD_npx_parquet_NCStdev|>
    dplyr::mutate(Threshold  = ifelse(NPX<LOD,"LOD","Above_LOD"))
  }
  
  # volcano plots, two plates combined
  pdf(file = "/Users/qiaoyanw/Desktop/LOD_VolcanoPlot.pdf",
      width = 11,
      height = 8)
  LOD_model |>
    filter(grepl('SC',SampleID, ignore.case = TRUE))|>
    ggplot(aes(x = (NPX-LOD), y = NPX, col = Threshold)) +
    geom_point() +
    set_plot_theme() +
    facet_wrap(.~SampleID)
  dev.off()
