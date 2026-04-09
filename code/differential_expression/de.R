library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

# Load object
seurat_obj <- readRDS("annotated.rds")

# Confirm metadata looks correct
cat("Time points:\n")
print(table(seurat_obj$orig.ident))

cat("\nMacrophage cell counts:\n")
print(table(seurat_obj$cell_type)[
  grep("Macrophage", unique(seurat_obj$cell_type), value = TRUE)
])

# Merge data
DefaultAssay(seurat_obj) <- "RNA"
all_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
seurat_obj[["RNA"]] <- CreateAssayObject(counts = all_counts)
seurat_obj <- NormalizeData(seurat_obj)

# Create subset to Macrophages
cluster_of_interest <- "Macrophages"

seurat_macro <- subset(seurat_obj,
                       subset = cell_type == cluster_of_interest
)

# Add Naive vs Infected condition label
seurat_macro$condition <- ifelse(
  seurat_macro$orig.ident == "Naive", "Naive", "Infected"
)
print(table(seurat_macro$orig.ident, seurat_macro$condition))

# Check a gene that should be in Infected cells (e.g., an Interferon or Isg15)
cat("\nChecking Isg15 counts in Infected cells:\n")
isg_check <- sum(GetAssayData(seurat_macro, slot = "counts")["Isg15", seurat_macro$condition == "Infected"])
print(isg_check)

# Set condition as active identity
Idents(seurat_macro) <- "condition"

# Build Wilcoxon DE
de_results <- FindMarkers(seurat_macro,
  ident.1         = "Infected",
  ident.2         = "Naive",
  test.use        = "wilcox",
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

de_results$gene <- rownames(de_results)
de_results <- de_results %>% arrange(p_val_adj)

# Classify significance
de_results$significance <- "Not Significant"
de_results$significance[
  de_results$avg_log2FC > 0.5 & de_results$p_val_adj < 0.05
] <- "Upregulated"
de_results$significance[
  de_results$avg_log2FC < -0.5 & de_results$p_val_adj < 0.05
] <- "Downregulated"

cat("DE Summary: Macrophages - Infected vs Naive \n")
cat("Total genes tested:  ", nrow(de_results), "\n")
cat("Upregulated:         ", sum(de_results$significance == "Upregulated"), "\n")
cat("Downregulated:       ", sum(de_results$significance == "Downregulated"), "\n")
cat("Not significant:     ", sum(de_results$significance == "Not Significant"), "\n")

cat("\nTop 10 upregulated:\n")
print(de_results %>%
  filter(significance == "Upregulated") %>%
  head(10) %>%
  select(gene, avg_log2FC, pct.1, pct.2, p_val_adj))

cat("\nTop 10 downregulated:\n")
print(de_results %>%
  filter(significance == "Downregulated") %>%
  head(10) %>%
  select(gene, avg_log2FC, pct.1, pct.2, p_val_adj))

# Save full results
write.csv(de_results,
  "DE_Macrophages_Infected_vs_Naive.csv",
  row.names = FALSE
)
cat("Results saved.\n")

# Volcano Plot
top_up <- de_results %>%
  filter(significance == "Upregulated") %>%
  head(15)
top_down <- de_results %>%
  filter(significance == "Downregulated") %>%
  head(15)
top_label <- bind_rows(top_up, top_down)

p_volcano <- ggplot(
  de_results,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    color = significance
  )
) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#D62728",
      "Downregulated" = "#1F77B4",
      "Not Significant" = "grey70"
    ),
    drop = FALSE
  ) +
  geom_text_repel(
    data = top_label,
    aes(label = gene),
    size = 2.8,
    max.overlaps = 20,
    color = "black",
    segment.size = 0.3
  ) +
  labs(
    title = "Macrophages: Infected vs Naive (Wilcoxon)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "bottom"
  )

ggsave("volcano_Macrophages.png",
  plot = p_volcano, width = 9, height = 7, bg = "white"
)
cat("Volcano plot saved.\n")


# Feature plots of top DE genes split by condition
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj$condition <- ifelse(
  seurat_obj$orig.ident == "Naive", "Naive", "Infected"
)

top_feat <- c(
  de_results %>% filter(significance == "Upregulated") %>% head(3) %>% pull(gene),
  de_results %>% filter(significance == "Downregulated") %>% head(3) %>% pull(gene)
)
cat("\nFeature plot genes:", paste(top_feat, collapse = ", "), "\n")

p_feat <- FeaturePlot(seurat_obj,
  features = top_feat,
  split.by = "condition",
  ncol     = 2,
  pt.size  = 0.1,
  order    = TRUE
) &
  theme(plot.title = element_text(size = 9, face = "bold"))

ggsave("featureplot_DE_Macrophages.png",
  plot   = p_feat,
  width  = 12,
  height = 3 * length(top_feat),
  bg     = "white"
)
cat("Feature plot saved.\n")

# Violin plots
seurat_macro$orig.ident <- factor(
  seurat_macro$orig.ident,
  levels = c("Naive", "D02", "D05", "D08", "D14")
)
Idents(seurat_macro) <- "orig.ident"

top6_vln <- de_results %>%
  filter(significance != "Not Significant") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(6) %>%
  pull(gene)

cat("\nViolin plot genes:", paste(top6_vln, collapse = ", "), "\n")

p_vln <- VlnPlot(seurat_macro,
  features = top6_vln,
  group.by = "orig.ident",
  ncol     = 3,
  pt.size  = 0,
  raster   = FALSE
) &
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(face = "bold", size = 10),
    axis.title.x = element_blank()
  )

ggsave("violin_Macrophages_timepoints.png",
  plot = p_vln, width = 14, height = 8, bg = "white"
)
cat("Violin plot saved.\n")

# Save
saveRDS(seurat_obj, file = "de.rds")
