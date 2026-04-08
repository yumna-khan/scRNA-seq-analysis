library(Seurat)
library(ggplot2)

# Load the normalized object
seurat_obj <- readRDS("seurat_processed_SCT.rds")

# Run PCA to produce principal components
seurat_obj <- RunPCA(seurat_obj,
                     assay   = "SCT",
                     npcs    = 50,
                     verbose = FALSE)

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
       bg     = "white")

# Check if samples separate in PCA for batch correction 
p_pca <- DimPlot(seurat_obj,
                 reduction = "pca",
                 group.by  = "orig.ident") +
  ggtitle("PCA — Coloured by Time Point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCA_by_timepoint.png",
       plot   = p_pca,
       width  = 8,
       height = 6,
       bg     = "white")
print(p_pca)

cat("PCA and Elbow Plot complete. Files saved.\n")
