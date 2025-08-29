library(msigdbr)

msig_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>%
  rbind(msigdbr::msigdbr(species = "Homo sapiens", category = "C5"))
msig_df <- msig_df %>%
  dplyr::select(gs_name, gene_symbol)
msig_db<-msig_df|>
  group_by(gs_name) %>%
  summarise(Total_gene = paste(gene_symbol, collapse = "/"))
colnames(msig_db) <-c("ID","Total_Genes")
write.csv(msig_db,'Desktop/MSig_DB.csv')
