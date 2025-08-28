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
  }
LOD_model<-LOD_npx_parquet|>
  dplyr::mutate(Threshold  = ifelse(NPX<LOD,"LOD","Non_LOD"))|>
  mutate(Threshold = replace_na(Threshold, "Non_LOD"))

if (LOD_opt==2){
  LOD_Assay <- npx_parquet %>% 
    filter(grepl('NEGATIVE_CONTROL',SampleType, ignore.case = TRUE)) %>%
    group_by(PlateID, Assay) %>%
    summarise(LOD=mean(NPX)+3*sd(NPX),PCNormalizedLOD=mean(PCNormalizedNPX)+3*sd(PCNormalizedNPX))
  LOD_npx_parquet_NCStdev<-select(LOD_Assay, PlateID, Assay, PCNormalizedLOD, LOD) %>%
    left_join(npx_parquet, by=c("Assay", "PlateID"))
}
LOD_model_NCStdev<-LOD_npx_parquet_NCStdev|>
  dplyr::mutate(Threshold  = ifelse(NPX<LOD,"LOD","Non_LOD"))
