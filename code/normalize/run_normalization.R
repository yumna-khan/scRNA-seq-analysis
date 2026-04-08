# ----------- Load Libraries ----------
library(Seurat)
library(SeuratObject)
library(sctransform)
library(glmGamPoi)

# ----------- HPC SETTINGS ----------
options(future.globals.maxSize = 400 * 1024^3)
library(future)
plan("multisession", workers = 16) # Parallel processing

# ---------- Normalize & Scale ----------
# Load the clean object
seurat_obj <- readRDS("seurat_CLEAN.rds")

# SCTransform 
seurat_obj <- SCTransform(seurat_obj,
                          vars.to.regress     = "percent.mt",
                          variable.features.n = 3000,
                          method              = "glmGamPoi",
                          verbose             = TRUE)

# Save the Progress
saveRDS(seurat_obj, file = "seurat_processed_SCT.rds")
cat("Normalization complete. File saved as seurat_processed_SCT.rds\n")

print(seurat_obj)

# View the top variable genes
top20 <- head(VariableFeatures(seurat_obj), 20)
print(top20)
