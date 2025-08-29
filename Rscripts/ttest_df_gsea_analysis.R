library(OlinkAnalyze)
library(dplyr)
library(arrow)
setwd('Desktop/HACKATHON/')
# read in files
Olink_datafile <- read_parquet("Olink_datafile_1.parquet")
Gender_datafile <- read.csv("SexManifest_1.csv", header = TRUE, sep = ",")
npx_parquet<-Olink_datafile |>
  filter(!grepl("CONTROL", SampleType))|>
  filter(!grepl("ctrl", AssayType))|>
  filter(!is.na(NPX))|>
  left_join(Gender_datafile, by=c("SampleID"))
                
ttest_results <- olink_ttest(
    df = npx_parquet,
    variable = "Gender",
    alternative = "two.sided")

Results<-as.data.frame(ttest_results)

df_to_label <- Results%>%
  dplyr::slice_head(n = 10)

library(ggrepel)
# volcano plots
pdf(file = "/Users/qiaoyanw/Desktop/ttest_plot.pdf",
    width = 11,
    height = 8)
Results |>
  ggplot(aes(x = estimate, y = -log10(Adjusted_pval), col = Threshold)) +
  geom_point() +
  geom_text_repel(data = df_to_label, aes(label = Assay), 
                  box.padding = 0.5, point.padding = 0.5, 
                  max.overlaps = Inf, # Allows more labels to be shown
                  segment.color = 'grey50')+
  set_plot_theme()
dev.off()

# pathway analysis

ontology_opt <- readline("Enter your ontology options (1 for KEGG, 2 for GO, 3 for Reactome, 4 for all): ")
if (ontology_opt==1){
  gsea_results <- olink_pathway_enrichment(data = npx_parquet, test_results = ttest_results, ontology = "KEGG")
}
if (ontology_opt==2){
  gsea_results <- olink_pathway_enrichment(data = npx_parquet, test_results = ttest_results, ontology = "GO")
}
if (ontology_opt==3){
  gsea_results <- olink_pathway_enrichment(data = npx_parquet, test_results = ttest_results, ontology = "Reactome")
}

if (ontology_opt==4){
  gsea_results <- olink_pathway_enrichment(data = npx_parquet, test_results = ttest_results)
}

pdf(file = "/Users/qiaoyanw/Desktop/Pathway_Plot.pdf",
    width = 11,
    height = 8)
# GSEA Heatmap from t-test results
olink_pathway_heatmap(enrich_results = gsea_results, test_results = ttest_results)
# GSEA bar chart from t-test results 
olink_pathway_visualization(enrich_results = gsea_results, method = "GSEA", number_of_terms = 15)
dev.off()

#output the pathway results
write_parquet(gsea_results, '/Users/qiaoyanw/Desktop/Pathway_results')
