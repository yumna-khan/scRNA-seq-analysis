# Load libraries
library(Seurat)
library(ggplot2)

# Load object
seurat_obj <- readRDS("seurat_PCA.rds")

# Define dimensions
dims_use <- 1:30

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, 
                            reduction = "pca",
                            dims = dims_use)

# Test two resolutions
seurat_obj <- FindClusters(seurat_obj, 
                           resolution = c(0.6, 0.9))

# UMAP with clusters displayed
seurat_obj <- RunUMAP(seurat_obj, 
                      dims = dims_use,
                      reduction = "pca")

# Plot Resolution 0.6
p1 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "SCT_snn_res.0.6", 
              label = TRUE, 
              raster = TRUE) + 
  ggtitle("IAV UMAP (Res 0.6)")

ggsave("UMAP_Res_0.6.png", plot = p1, width = 10, height = 8, bg = "white")

# Plot Resolution 0.9
p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "SCT_snn_res.0.9", 
              label = TRUE, 
              raster = TRUE) + 
  ggtitle("IAV UMAP (Res 0.9)")

ggsave("UMAP_Res_0.9.png", plot = p2, width = 10, height = 8, bg = "white")

# Save the Progress
saveRDS(seurat_obj, file = "seurat_UMAP.rds")
cat("Complete. File saved as seurat_UMAP.rds\n")
