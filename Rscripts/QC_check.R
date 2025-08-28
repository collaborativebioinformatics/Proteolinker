Olink_datafile <- read_parquet("Desktop/Olink_datafile.parquet")

Olink_datafile%>%
  OlinkAnalyze:::check_data_completeness()

npx_parquet<-Olink_datafile

AssayWarn <- npx_parquet %>% 
  filter(grepl('WARN',AssayQC, ignore.case = TRUE))
SampleWarn <- npx_parquet %>% 
  filter(grepl('WARN',SampleQC, ignore.case = TRUE))
SampleFail <- npx_parquet %>% 
  filter(grepl('FAIL',SampleQC, ignore.case = TRUE))

AssayWarnData <-unique(select(AssayWarn,c('PlateID','Assay','Block','OlinkID'))) |>
  group_by(Assay,Block, OlinkID) |>
  summarise(across(everything(), ~toString(.)))

SampleWarnData <-unique(select(SampleWarn,c('SampleID','Block', 'PlateID'))) |>
  group_by(SampleID, PlateID) |>
  summarise(across(everything(), ~toString(.))) 

SampleFailData <-unique(select(SampleFail,c('SampleID','Block', 'PlateID'))) |>
  group_by(SampleID,PlateID) |>
  summarise(across(everything(), ~toString(.)))

# print out tables to the pdf for failed and warned items
pdf(file = "/Users/qiaoyanw/Desktop/QC.pdf",
    width = 11,
    height = 8)
par(mfrow = c(3,3))
table_grob <- tableGrob(SampleFailData)
title_grob <- textGrob("Sample Failed", gp = gpar(fontsize = 16, fontface = "bold"))
arranged_plot_1 <- grid.arrange(title_grob, table_grob, ncol = 1, heights = c(0.1, 0.9))
table_grob <- tableGrob(SampleWarnData)
title_grob <- textGrob("Sample Warned", gp = gpar(fontsize = 16, fontface = "bold"))
arranged_plot_2 <- grid.arrange(title_grob, table_grob, ncol = 1, heights = c(0.1, 0.9))
table_grob <- tableGrob(AssayWarnData)
title_grob <- textGrob("Assay Warned", gp = gpar(fontsize = 16, fontface = "bold"))
arranged_plot_3 <- grid.arrange(title_grob, table_grob, ncol = 1, heights = c(0.1, 0.9))
dev.off()
