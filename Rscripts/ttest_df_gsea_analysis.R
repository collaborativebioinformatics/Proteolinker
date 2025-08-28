library(OlinkAnalyze)
library(dplyr)
library(arrow)
setwd('Desktop/HACKATHON/')
# read in files
npx_parquet <- read_parquet("Olink_datafile_1.parquet")
manifest_file <- read.csv("SexManifest_1.csv", header = TRUE, sep = ",")
# merge Gender info to npx dataframe
npx_with_manifest <- merge(npx_parquet, manifest_file[, c("SampleID", "Gender")], by = "SampleID", all.x = TRUE)
# remove controls
npx_df <- npx_with_manifest |>
     dplyr::filter(!grepl("control", SampleType, ignore.case = TRUE))
ttest_results <- olink_ttest(
     df = npx_df,
     variable = "Gender",
     alternative = "two.sided")
# pathway analysis
gsea_results <- olink_pathway_enrichment(data = npx_df, test_results = ttest_results)

