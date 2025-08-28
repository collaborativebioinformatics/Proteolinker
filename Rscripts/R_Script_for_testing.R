# dummy R script for testing
# may need to adjust code to reflect file paths
#usage: Rscript R_Script_for_testing.R

# load libraries
library(OlinkAnalyze)
library(dplyr)

# Read the parquet file
olink_data <- read_parquet("Olink_datafile_1.parquet")

# Get distinct sample IDs
unique_ids <- olink_data %>%
  distinct(SampleID)

#write out to a file so you can make the replace id file below
write.table(unique_ids, "test_output.txt", row.names = FALSE, quote = FALSE, sep = "\t")
