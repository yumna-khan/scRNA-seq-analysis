library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load object
seurat_obj <- readRDS("seurat_UMAP.rds")

# Force the object to use the 0.6 resolution column
Idents(seurat_obj) <- "SCT_snn_res.0.6"

# Find top marker genes per cluster (manual annotation)
all_markers <- FindAllMarkers(seurat_obj,
                              only.pos        = TRUE,
                              min.pct         = 0.25,
                              logfc.threshold = 0.25)

top5_genes <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print(top5_genes)

# Save
write.csv(top5_genes, "top5_cluster.csv", row.names = FALSE)
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

# Feature plots
p_feat <- FeaturePlot(seurat_obj,
                      features = c("Adcy3", "Gap43", "Neurod1",
                                   "Krt15", "Ascl3", "Muc2",
                                   "Vmo1", "Trpm5", "C1qa",
                                   "Ms4a1", "Col2a1", "Emcn"),
                      ncol    = 4,
                      pt.size = 0.2,
                      order   = TRUE)

ggsave("featureplot_markers.png",
       plot = p_feat, width = 18, height = 9, bg = "white")

# Manual annotation
# Based on marker gene analysis of all_cluster_markers.csv
new_labels <- c(
  "0"  = "Olfactory Sensory Neurons",    # Kirrel3, Calb2, Rims3
  "1"  = "Horizontal Basal Cells",       # Krt15, Serpinb5, Defb1 (Nasal Stem Cells)
  "2"  = "Olfactory Sensory Neurons",    # Cabp4, Nxnl2, Pcdhb1
  "3"  = "Mature Olfactory Neurons",     # Adcy3, Umodl1, Cnga4 (Canonical OSN)
  "4"  = "Olfactory Sensory Neurons",    # S100a5, Kirrel2, Rgs7
  "5"  = "Macrophages",                  # C1qa, C1qb, C1qc (Complement)
  "6"  = "B Cells",                      # Ms4a1, Ighd, Fcer2a
  "7"  = "Olfactory Sensory Neurons",    # Pcp4l1, S100a5, Dlg2
  "8"  = "NK Cells",                     # Ncr1, Klrb1c, Prf1, Gzma
  "9"  = "Immature Neurons",             # Gap43, Prph (Growth/Axon markers)
  "10" = "Endothelial Cells",            # Emcn, Gpihbp1, Rbp7
  "11" = "Neutrophils",                  # Cxcr2, Clec4d, Acod1
  "12" = "Olfactory Ensheathing Cells",  # Nrcam, Abca4, Car10
  "13" = "Tissue Resident Macrophages", # Adgre4, Treml4, Ear2
  "14" = "Fibroblasts",                  # Mfap4, Apod, Dpt
  "15" = "Dendritic Cells",              # Cd209a, Flt3, Mgl2
  "16" = "M2-like Macrophages",         # Mrc1, Pf4, P2ry12, Fcrls
  "17" = "Bowman Gland Cells",           # Reg3g, Bpifa1, Tff2, Gp2
  "18" = "Olfactory Gland Cells",        # Vmo1, Chil6, Bpifb4, Gpx6
  "19" = "Secretory Cells",              # Tac1, Wfdc18, Barx2
  "20" = "Monocytes",                    # Ly6i, Ifi205, Ms4a4c
  "21" = "Ionocytes",                    # Ascl3, Clcnka, Kl
  "22" = "Neuronal Progenitors",         # Neurod1, Neurog1 (Neurogenic TFs)
  "23" = "Smooth Muscle Cells",          # Myh11, Pln, Higd1b
  "24" = "Immature B Cells",             # Rag1, Vpreb3, Pou2af1
  "25" = "Activated Neutrophils",        # Ly6g, Ngp, Lcn2, Mmp8
  "26" = "Goblet Cells",                 # Muc2, Sec14l3, Gldn
  "27" = "Serous Acinar Cells",          # Bpifb9a/b (Nasal serous proteins)
  "28" = "Tuft Cells",                   # Trpm5, Il25, Avil
  "29" = "Serous Gland Cells",           # Car6, Scgb2b27, Csn3
  "30" = "Osteoblasts",                  # Bglap/Bglap2, Ibsp (Bone)
  "31" = "Myeloid Progenitors",          # Mpo, Ctsg, Elane
  "32" = "Microvillar Cells",            # Dnase2b, Mal, Mfsd2a
  "33" = "Schwann Cells",                # Mpz, Gpr37l1, Foxd3
  "34" = "Neutrophil Progenitors",       # Ms4a3, Fcnb, Elane
  "35" = "Ciliated Cells",               # Tmem212, Sntn (Motile cilia)
  "36" = "Chondrocytes",                 # Col2a1, Snorc, Ucma (Cartilage)
  "37" = "Squamous Epithelial Cells",    # Csta1, Crct1, Cysrt1
  "38" = "Neuronal Stem Cells"           # Fezf2, Meg3, Rian
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
                   pt.size    = 0.3) +
  ggtitle("Manual Cell Type Annotations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 2,
                              override.aes = list(size = 3)))

ggsave("UMAP_annotated.png",
       plot = p_annot, width = 14, height = 10, bg = "white")

# Unannotated cluster UMAP
p_clusters <- DimPlot(seurat_obj,
                      reduction = "umap",
                      group.by  = "SCT_snn_res.0.6",
                      label     = TRUE,
                      label.size = 3,
                      raster    = FALSE,
                      pt.size   = 0.3) +
  ggtitle("Seurat Clusters (Res 0.6)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 2,
                              override.aes = list(size = 3)))

# Side by side comparison
ggsave("UMAP_clusters_vs_annotations.png",
       plot   = p_clusters + p_annot,
       width  = 24,
       height = 10,
       bg     = "white")

# UMAP coloured by time point to check if biology drives clustering
p_time <- DimPlot(seurat_obj,
                  reduction = "umap",
                  group.by  = "orig.ident",
                  raster    = FALSE,
                  pt.size   = 0.3) +
  ggtitle("UMAP — Time Point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("UMAP_by_timepoint.png",
       plot = p_time, width = 10, height = 8, bg = "white")

# Save
saveRDS(seurat_obj, file = "annotated.rds")
cat("Complete. Saved as annotated.rds\n")
