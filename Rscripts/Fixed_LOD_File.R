library(dplyr)
library(readr)

# Insert your input file path
file_path <- "insert_path_here"

# Read your original fixed LOD file (semicolon-separated)
fixed_lod_original <- read_delim(file_path, delim = ";")

# Convert to proper format
fixed_lod_proper <- fixed_lod_original %>%
  select(Assay, LOD = LODNPX) %>%
  distinct(Assay, .keep_all = TRUE)

# Insert your output file path
output_path <- "insert_output_path_here"

# Save the corrected file (comma-separated)
write_csv(fixed_lod_proper, output_path)

# View the first few rows to confirm
head(fixed_lod_proper)
