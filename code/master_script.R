# Yumna Khan
# April 1, 2026
# scRNA-seq Analysis: Influenza A Virus (IAV) Infection in Mouse Lung

# ---------- General Overview ----------
# This R master script provides a complete pipeline for scRNA-seq analysis of IAV-infected respiratory tissues, covering quality control, SCTransform normalization, unsupervised clustering, and manual annotation of 38 cell types. 
# The script specifically investigates the macrophage response through Wilcoxon-based differential expression and Gene Set Enrichment Analysis (GSEA) to map activated antiviral pathways. 
# To accommodate the high memory and processing demands of large-scale transcriptomic mapping and iterative clustering, the majority of this workflow was executed on a High-Performance Computing (HPC) cluster.


# ---------- Load libraries ----------
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# ----------- HPC SETTINGS ----------
options(future.globals.maxSize = 400 * 1024^3) # Increases limit to 400GB
library(future)
plan("multisession", workers = 16) # Parallel processing

# ---------- Data Exploration ----------
# Load and print seurat object
seurat_obj <- readRDS("seurat_ass4.rds")

print(seurat_obj)

# Read row and column names
head(rownames(seurat_obj))
head(colnames(seurat_obj))

cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")

# Inspect metadata — find the column holding sample / condition labels
head(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)

# Check infectious status
unique(seurat_obj$disease__ontology_label)

# Identify organ regions
table(seurat_obj$organ_custom)

# Identify time points
unique(seurat_obj$time)

unique(seurat_obj$orig.ident)

# Table for different time points samples
table(seurat_obj$organ_custom, seurat_obj$time)


# ---------- Quality Control ----------
# Calculate percentage of counts for mitochondria
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Summary statistics
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)

# Violin plots BEFORE filtering
vln_p <- VlnPlot(seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  ncol = 3,
  pt.size = 0,
  raster = FALSE
)

# Set up axes and title
vln_p <- vln_p &
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8)
  ) &
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Assign specific titles and axis labels to each panel
vln_p[[1]] <- vln_p[[1]] + labs(
  title = "Gene Count",
  y = "Number of Unique Genes"
)

vln_p[[2]] <- vln_p[[2]] + labs(
  title = "RNA Count",
  y = "Total RNA Molecules"
)

vln_p[[3]] <- vln_p[[3]] + labs(
  title = "Mitochondrial %",
  y = "Percentage of Reads (%)"
)

# Plot
vln_p

ggsave("QC_violin_prefilter.png", width = 14, height = 5)


# Number of cells BEFORE filtering
cells_before <- ncol(seurat_obj)

# Scatter plot for doublets and percent mt. (used for subsetting)
p1 <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident") +
  labs(
    x = "Total RNA Counts",
    y = "Number of Genes",
    title = "RNA Counts vs Gene Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

p2 <- FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt", group.by = "orig.ident") +
  labs(
    x = "Total RNA Counts",
    y = "Mitochondrial Percentage (%)",
    title = "RNA Counts vs Mitochondrial Content"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

p1 + p2

ggsave("QC_scatter_prefilter.png", width = 12, height = 5)

# Calculate the pearson correlation coefficient
cor(
  seurat_obj$nCount_RNA,
  seurat_obj$nFeature_RNA,
  method = "pearson"
)

cor(
  seurat_obj$nCount_RNA,
  seurat_obj$percent.mt,
  method = "pearson"
)

# Subset to new seurat object
seurat_obj_filtered <- subset(seurat_obj,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    nCount_RNA < 6500 &
    percent.mt < 10
)

# Number of cells AFTER filtering
cells_after <- ncol(seurat_obj_filtered)

# Print before and after cell count differences
cat("Cells before QC:", cells_before, "\n")
cat("Cells after QC: ", cells_after, "\n")
cat("Cells removed:  ", cells_before - cells_after, "\n")

# Violin plots AFTER filtering
vln_p <- VlnPlot(seurat_obj_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  ncol = 3,
  pt.size = 0,
  raster = FALSE
)

# Set up axes and title
vln_p <- vln_p &
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8)
  ) &
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Assign specific titles and axis labels to each panel
vln_p[[1]] <- vln_p[[1]] + labs(
  title = "Gene Count",
  y = "Number of Unique Genes"
)

vln_p[[2]] <- vln_p[[2]] + labs(
  title = "RNA Count",
  y = "Total RNA Molecules"
)

vln_p[[3]] <- vln_p[[3]] + labs(
  title = "Mitochondrial %",
  y = "Percentage of Reads (%)"
)

# Plot
vln_p

ggsave("QC_violin_postfilter.png", width = 14, height = 5)

# Save clean version of object for HPC
saveRDS(seurat_obj_filtered, file = "seurat_CLEAN.rds")

# ---------- Normalize & Scale ----------
# Load the clean object
seurat_obj <- readRDS("seurat_CLEAN.rds")

# SCTransform
seurat_obj <- SCTransform(seurat_obj,
  vars.to.regress     = "percent.mt",
  variable.features.n = 3000,
  method              = "glmGamPoi",
  verbose             = TRUE
)

# Save the Progress
saveRDS(seurat_obj, file = "seurat_processed_SCT.rds")
cat("Normalization complete. File saved as seurat_processed_SCT.rds\n")

print(seurat_obj)

# View the top variable genes
top20 <- head(VariableFeatures(seurat_obj), 20)
print(top20)


# ---------- Clustering ----------
# Load the normalized object
seurat_obj <- readRDS("seurat_processed_SCT.rds")

# Run PCA to produce principal components
seurat_obj <- RunPCA(seurat_obj,
  assay   = "SCT",
  npcs    = 50,
  verbose = FALSE
)

# Visualize the PC
message("Generating Elbow Plot...")
p_elbow <- ElbowPlot(seurat_obj, ndims = 50) +
  labs(
    x = "Number of PCs",
    y = "Standard Deviation",
    title = "Elbow Plot: Number of PCs vs. Standard Deviation"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave("elbow_plot.png",
  plot   = p_elbow,
  width  = 8,
  height = 5,
  bg     = "white"
)

# Check if samples separate in PCA for batch correction
p_pca <- DimPlot(seurat_obj,
  reduction = "pca",
  group.by  = "orig.ident"
) +
  ggtitle("PCA — Coloured by Time Point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCA_by_timepoint.png",
  plot   = p_pca,
  width  = 8,
  height = 6,
  bg     = "white"
)
print(p_pca)

# Save the Progress
saveRDS(seurat_obj, file = "seurat_PCA.rds")
cat("PCA and Elbow Plot complete. Files saved.\n")

# Load the object
seurat_obj <- readRDS("seurat_PCA.rds")

# Set dimensions after inspecting the elbow plot
dims_use <- 1:30

# Clustering
seurat_obj <- FindNeighbors(seurat_obj,
  reduction = "pca",
  dims = dims_use
)

# Test two resolutions
seurat_obj <- FindClusters(seurat_obj,
  resolution = c(0.6, 0.9)
)

# UMAP with clusters displayed
seurat_obj <- RunUMAP(seurat_obj,
  dims = dims_use,
  reduction = "pca"
)

# Plot Resolution 0.6
p1 <- DimPlot(seurat_obj,
  reduction = "umap",
  group.by = "SCT_snn_res.0.6",
  label = TRUE,
  raster = TRUE
) +
  ggtitle("IAV UMAP (Res 0.6)")

ggsave("UMAP_Res_0.6.png", plot = p1, width = 10, height = 8, bg = "white")

# Plot Resolution 0.9
p2 <- DimPlot(seurat_obj,
  reduction = "umap",
  group.by = "SCT_snn_res.0.9",
  label = TRUE,
  raster = TRUE
) +
  ggtitle("IAV UMAP (Res 0.9)")

ggsave("UMAP_Res_0.9.png", plot = p2, width = 10, height = 8, bg = "white")

# Save the Progress
saveRDS(seurat_obj, file = "seurat_UMAP.rds")
cat("Complete. File saved as seurat_UMAP.rds\n")

# ---------- Annotation ----------
# Load object
seurat_obj <- readRDS("seurat_UMAP.rds")

# Force the object to use the 0.6 resolution column
Idents(seurat_obj) <- "SCT_snn_res.0.6"

# Find top marker genes per cluster (manual annotation)
all_markers <- FindAllMarkers(seurat_obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

top5_genes <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print(top5_genes)

# Save
write.csv(top5_genes, "top5_cluster.csv", row.names = FALSE)
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

# Feature plots
p_feat <- FeaturePlot(seurat_obj,
  features = c(
    "Adcy3", "Gap43", "Neurod1",
    "Krt15", "Ascl3", "Muc2",
    "Vmo1", "Trpm5", "C1qa",
    "Ms4a1", "Col2a1", "Emcn"
  ),
  ncol = 4,
  pt.size = 0.2,
  order = TRUE
)

ggsave("featureplot_markers.png",
  plot = p_feat, width = 18, height = 9, bg = "white"
)

# Manual annotation
# Based on marker gene analysis of all_cluster_markers.csv
new_labels <- c(
  "0"  = "Olfactory Sensory Neurons", # Kirrel3, Calb2, Rims3
  "1"  = "Horizontal Basal Cells", # Krt15, Serpinb5, Defb1 (Nasal Stem Cells)
  "2"  = "Olfactory Sensory Neurons", # Cabp4, Nxnl2, Pcdhb1
  "3"  = "Mature Olfactory Neurons", # Adcy3, Umodl1, Cnga4 (Canonical OSN)
  "4"  = "Olfactory Sensory Neurons", # S100a5, Kirrel2, Rgs7
  "5"  = "Macrophages", # C1qa, C1qb, C1qc (Complement)
  "6"  = "B Cells", # Ms4a1, Ighd, Fcer2a
  "7"  = "Olfactory Sensory Neurons", # Pcp4l1, S100a5, Dlg2
  "8"  = "NK Cells", # Ncr1, Klrb1c, Prf1, Gzma
  "9"  = "Immature Neurons", # Gap43, Prph (Growth/Axon markers)
  "10" = "Endothelial Cells", # Emcn, Gpihbp1, Rbp7
  "11" = "Neutrophils", # Cxcr2, Clec4d, Acod1
  "12" = "Olfactory Ensheathing Cells", # Nrcam, Abca4, Car10
  "13" = "Tissue Resident Macrophages", # Adgre4, Treml4, Ear2
  "14" = "Fibroblasts", # Mfap4, Apod, Dpt
  "15" = "Dendritic Cells", # Cd209a, Flt3, Mgl2
  "16" = "M2-like Macrophages", # Mrc1, Pf4, P2ry12, Fcrls
  "17" = "Bowman Gland Cells", # Reg3g, Bpifa1, Tff2, Gp2
  "18" = "Olfactory Gland Cells", # Vmo1, Chil6, Bpifb4, Gpx6
  "19" = "Secretory Cells", # Tac1, Wfdc18, Barx2
  "20" = "Monocytes", # Ly6i, Ifi205, Ms4a4c
  "21" = "Ionocytes", # Ascl3, Clcnka, Kl
  "22" = "Neuronal Progenitors", # Neurod1, Neurog1 (Neurogenic TFs)
  "23" = "Smooth Muscle Cells", # Myh11, Pln, Higd1b
  "24" = "Immature B Cells", # Rag1, Vpreb3, Pou2af1
  "25" = "Activated Neutrophils", # Ly6g, Ngp, Lcn2, Mmp8
  "26" = "Goblet Cells", # Muc2, Sec14l3, Gldn
  "27" = "Serous Acinar Cells", # Bpifb9a/b (Nasal serous proteins)
  "28" = "Tuft Cells", # Trpm5, Il25, Avil
  "29" = "Serous Gland Cells", # Car6, Scgb2b27, Csn3
  "30" = "Osteoblasts", # Bglap/Bglap2, Ibsp (Bone)
  "31" = "Myeloid Progenitors", # Mpo, Ctsg, Elane
  "32" = "Microvillar Cells", # Dnase2b, Mal, Mfsd2a
  "33" = "Schwann Cells", # Mpz, Gpr37l1, Foxd3
  "34" = "Neutrophil Progenitors", # Ms4a3, Fcnb, Elane
  "35" = "Ciliated Cells", # Tmem212, Sntn (Motile cilia)
  "36" = "Chondrocytes", # Col2a1, Snorc, Ucma (Cartilage)
  "37" = "Squamous Epithelial Cells", # Csta1, Crct1, Cysrt1
  "38" = "Neuronal Stem Cells" # Fezf2, Meg3, Rian
)

seurat_obj <- RenameIdents(seurat_obj, new_labels)
seurat_obj$cell_type <- Idents(seurat_obj)

cat("\nCell type distribution:\n")
print(sort(table(seurat_obj$cell_type), decreasing = TRUE))


# UMAP plots
# Annotated UMAP with labels
p_annot <- DimPlot(seurat_obj,
  reduction  = "umap",
  group.by   = "cell_type",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3,
  raster     = FALSE,
  pt.size    = 0.3
) +
  ggtitle("Manual Cell Type Annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 3)
  ))

ggsave("UMAP_annotated.png",
  plot = p_annot, width = 14, height = 10, bg = "white"
)

# Unannotated cluster UMAP
p_clusters <- DimPlot(seurat_obj,
  reduction = "umap",
  group.by = "SCT_snn_res.0.6",
  label = TRUE,
  label.size = 3,
  raster = FALSE,
  pt.size = 0.3
) +
  ggtitle("Seurat Clusters (Res 0.6)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 3)
  ))

# Side by side comparison
ggsave("UMAP_clusters_vs_annotations.png",
  plot   = p_clusters + p_annot,
  width  = 24,
  height = 10,
  bg     = "white"
)

# UMAP coloured by time point to check if biology drives clustering
p_time <- DimPlot(seurat_obj,
  reduction = "umap",
  group.by  = "orig.ident",
  raster    = FALSE,
  pt.size   = 0.3
) +
  ggtitle("UMAP — Time Point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("UMAP_by_timepoint.png",
  plot = p_time, width = 10, height = 8, bg = "white"
)

# Save
saveRDS(seurat_obj, file = "annotated.rds")
cat("Complete. Saved as annotated.rds\n")

# ---------- DE ----------
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
  seurat_macro$orig.ident == "Naive", 
  "Naive", "Infected"
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

# ---------- GSEA ----------
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
