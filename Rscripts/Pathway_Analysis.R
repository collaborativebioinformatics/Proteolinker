# Install if needed
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", 
#                        "enrichplot", "ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)     # Gene annotation
library(enrichplot)       # Visualization
library(ReactomePA)       # Reactome pathways
library(dplyr)

# Must contain at least: GeneSymbol, log2FC, pval

de_results <- read.csv("DE_results.csv")

# Step 1: Create a ranked list
de_results <- de_results %>%
  mutate(rank_stat = -log10(pval) * sign(log2FC))

gene_list <- de_results$rank_stat
names(gene_list) <- de_results$GeneSymbol

gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(names(gene_list))]


# Step 2: Function to run enrichment
run_enrichment <- function(method = c("GO", "KEGG", "Reactome", "All")) {
  method <- match.arg(method)
  
  if (method == "GO") {
    res <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "BP",   # Biological Process
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  } else if (method == "KEGG") {
    # KEGG needs ENTREZ IDs
    gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_list_entrez <- gene_list[gene_df$SYMBOL]
    names(gene_list_entrez) <- gene_df$ENTREZID
    
    res <- gseKEGG(geneList     = gene_list_entrez,
                   organism     = "hsa",   # human
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
  } else if (method == "Reactome") {
    gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_list_entrez <- gene_list[gene_df$SYMBOL]
    names(gene_list_entrez) <- gene_df$ENTREZID
    
    res <- gsePathway(geneList     = gene_list_entrez,
                      organism     = "human",
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
  } else if (method == "All") {
    res_go <- run_enrichment("GO")
    res_kegg <- run_enrichment("KEGG")
    res_reactome <- run_enrichment("Reactome")
    
    res <- list(GO = res_go,
                KEGG = res_kegg,
                Reactome = res_reactome)
  }
  
  return(res)
}


# Step 3: Run & visualize
# 1 "GO" (2)"KEGG" (3) "Reactome" (4)"All"
#gsea_res <- run_enrichment("GO") 

if (!is.list(gsea_res)) {
  dotplot(gsea_res, showCategory = 20)
}


# Step 4: Visualization
visualize_results <- function(res, method) {
  if (is.list(res)) {
    par(mfrow = c(1, 3))
    dotplot(res$GO, showCategory = 10, title = "GO BP")
    dotplot(res$KEGG, showCategory = 10, title = "KEGG")
    dotplot(res$Reactome, showCategory = 10, title = "Reactome")
    par(mfrow = c(1, 1))
  } else {
    print(dotplot(res, showCategory = 20, title = paste(method, "Enrichment")))
    print(heatplot(res, foldChange = gene_list, showCategory = 20))
  }
}

# Example
#visualize_results(gsea_res, "GO")



