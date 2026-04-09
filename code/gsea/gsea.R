library(ggplot2)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# Load DE results
de_results <- read.csv("DE_Macrophages_Infected_vs_Naive.csv")

cat("DE results loaded:", nrow(de_results), "genes\n")
cat("Upregulated:  ", sum(de_results$significance == "Upregulated"), "\n")
cat("Downregulated:", sum(de_results$significance == "Downregulated"), "\n")

# Convert gene symbols to Entrez IDs
gene_map <- bitr(
  de_results$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

cat("Genes successfully mapped:", nrow(gene_map), "/", nrow(de_results), "\n")

# Merge Entrez IDs back into DE results
de_results <- de_results %>%
  left_join(gene_map, by = c("gene" = "SYMBOL"))

# Build ranked gene list descending by avg_log2FC
# Higher rank = more upregulated in Infected vs Naive
ranked_df <- de_results %>%
  filter(!is.na(ENTREZID), !is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))

gene_list <- setNames(ranked_df$avg_log2FC, ranked_df$ENTREZID)
gene_list <- gene_list[!duplicated(names(gene_list))] 

cat("Genes in ranked list:", length(gene_list), "\n")
cat("log2FC range:", round(min(gene_list), 2),
    "to", round(max(gene_list), 2), "\n")

# --- GSEA on GO Biological Process ---
gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

cat("Significant GO terms (GSEA):",
    nrow(as.data.frame(gsea_go)), "\n")

# Plot
p_gsea_go <- dotplot(gsea_go,
                     showCategory = 15,
                     split        = ".sign") +
  facet_grid(~.sign) +
  ggtitle("GSEA - GO Biological Process: Infected vs Naive Macrophages") +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(size = 7),
    strip.text  = element_text(face = "bold", size = 9)
  )

ggsave("GSEA_GO_dotplot.png", 
       plot = p_gsea_go, width = 16, height = 9, bg = "white")
cat("GSEA GO dot plot saved.\n")

# --- GSEA on KEGG pathways ---
gsea_kegg <- gseKEGG(
  geneList      = gene_list,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

cat("Significant KEGG pathways (GSEA):",
    nrow(as.data.frame(gsea_kegg)), "\n")

# Plot
p_gsea_kegg <- dotplot(gsea_kegg,
                       showCategory = 15,
                       split        = ".sign") +
  facet_grid(~.sign) +
  ggtitle("GSEA - KEGG Pathways: Infected vs Naive Macrophages") +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(size = 7),
    strip.text  = element_text(face = "bold", size = 9)
  )

ggsave("GSEA_KEGG_dotplot.png", 
       plot = p_gsea_kegg, width = 16, height = 8, bg = "white")
cat("GSEA KEGG dot plot saved.\n")
