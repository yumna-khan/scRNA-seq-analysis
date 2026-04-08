# scRNA-seq

## General Overview


## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
  - [1. Data Description](#1-data-description)
  - [2. Quality Control](#2-quality-control)
  - [3. Normalization & Scaling](#3-normalization--scaling)
  - [4. Clustering](#4-clustering)
  - [5. Annotation](#5-annotation)
  - [6. Differential Expression](#6-differential-expression)
  - [7. GSEA](#7-gsea)
- [Results](#results)
- [Discussion](#discussion)
- [References](#references)

## Introduction
Influenza A virus (IAV) is a highly contagious respiratory pathogen that causes seasonal flu outbreaks, resulting in substantial global morbidity and mortality. IAV infections are responsible for hundreds of thousands of respiratory deaths annually, placing significant strain on healthcare systems during peak transmission seasons (Influenza (Seasonal), 2025; Mansell & Tate, 2017). Understanding cellular responses to IAV infection is therefore critical for elucidating mechanisms of viral pathogenesis and host defense. At the molecular level, these responses are reflected in changes in gene expression across different cell types.

Although most cells in the body contain the same genetic material, their transcriptomes differ, reflecting cell-specific gene activity that governs identity, function, and responses to environmental cues. Profiling gene expression at the single-cell level provides the most accurate insight into these cellular differences, enabling the study of cell states, lineage, and functional responses. Over the past two decades, single-cell transcriptomics has emerged as a transformative technology, allowing high-resolution analysis of gene expression in individual cells (Jovic et al., 2022).

Building on this, single-cell RNA sequencing (scRNA-seq) enables the dissection of transcriptional heterogeneity, identifying distinct cell types, and capturing their dynamic responses during infection. This technology is particularly valuable for understanding how infections alter tissue composition and immune memory, informing the development of mucosal therapeutics and vaccines against respiratory viral pathogens (Kazer et al., 2024). The nasal mucosa, which plays roles in filtration, air conditioning, and olfaction, is a critical site for viral entry and host defense.

For scRNA-seq data analysis, Seurat provides a robust framework encompassing preprocessing, quality control, dimensionality reduction, clustering, and visualization. Quality control typically involves filtering cells based on the number of detected genes, total RNA counts, and the proportion of mitochondrial gene expression (Seurat – Guided Clustering Tutorial, 2023). Dimensionality reduction using PCA is typically performed prior to UMAP, which is preferred over alternatives such as t-SNE due to its computational efficiency and ability to preserve both local and global structure across heterogeneous cell populations (Marx, 2024). Moreover, feature plots allow visualization of gene expression at the single-cell level, facilitating validation of cluster identity.

Furthermore, annotation is a critical step for interpreting biological meaning from clusters. Manual annotation using known marker genes enables accurate cell-type identification, while automated approaches such as Seurat’s anchor-based label transfer are robust for datasets with lower sequencing depth. For fine-grained classification, particularly among closely related immune subtypes, SingleR offers high sensitivity and complements Seurat’s annotations (Pasquini et al., 2021).

Once cell types have been annotated, differential expression analysis of specific clusters enables the identification of genes and pathways altered in response to infection. Tools such as edgeR and DESeq2 provide statistical frameworks for differential expression analysis, although their performance may be reduced at moderate sequencing depth (Nguyen et al., 2023). Functional interpretation of differentially expressed genes can be enhanced through Gene Set Enrichment Analysis (GSEA) or Over-Representation Analysis (ORA), providing insight into the biological processes underlying observed changes.

Overall, the objective of this study is to characterize the transcriptional response of nasal mucosa cell populations to IAV infection in mice, leveraging scRNA-seq to identify distinct cell types, annotate clusters, and investigate infection-induced changes at both the gene and pathway levels.


## Methods
### 1. Data Description
The dataset used in this study was obtained from a previously published study by Kazer et al., 2024. It consists of single-cell RNA sequencing data from mouse nasal mucosa following Influenza A virus (IAV) infection, where viral exposure was restricted to the nasal cavity.

The dataset includes cells collected from distinct anatomical regions, including the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LNG). Samples were obtained across five time points: naïve (uninfected) and 2, 5, 8, and 14 days post-infection (dpi).

In total, the dataset comprises 156,572 cells and 25,129 gene features. Cells are categorized into two infection states: naïve and infected, enabling comparative analysis of transcriptional responses over time.

### 2. Quality Control
The dataset was explored prior to quality control to examine metadata and assess the number of cells and genes. Quality control was performed by calculating the percentage of mitochondrial gene expression (`percent.mt`) as an indicator of cellular stress. Cells undergoing stress or apoptosis, which can occur during viral infections such as Influenza A virus (IAV), often exhibit elevated mitochondrial RNA content due to leakage of cytoplasmic RNA. Cells with mitochondrial percentages exceeding 10–20% were considered low quality and potentially non-viable (Osorio & Cai, 2020).

To visualize the distribution of quality control metrics prior to filtering, violin plots were generated for each sample (grouped by `orig.ident`), displaying the number of detected genes per cell (`nFeature_RNA`), total RNA counts per cell (`nCount_RNA`), and the percentage of mitochondrial gene expression (`percent.mt`).

Additionally, scatter plots were used to assess relationships between quality control metrics. A strong positive correlation was observed between `nCount_RNA` and `nFeature_RNA` (Pearson correlation = 0.83), indicating that cells with higher RNA counts tend to have more detected genes. In contrast, a weak correlation was observed between mitochondrial percentage and gene count (Pearson correlation = -0.06), suggesting that mitochondrial content is largely independent of gene complexity.

Based on these metrics, cells were filtered to retain high-quality cells with `percent.mt < 10, nFeature_RNA > 200, and nFeature_RNA < 4000`. A new Seurat object containing only filtered cells was then created for downstream analysis.

### 3. Normalization & Scaling
Prior to normalization, the filtered Seurat object was transferred to a high-performance computing (HPC) environment for efficient processing. Normalization was performed using the `SCTransform` method implemented in the Seurat package. SCTransform performs normalization, scaling, and variance stabilization in a single step using a regularized negative binomial regression model, which accounts for sequencing depth and technical variability (Hafemeister & Satija, 2019).

During normalization, unwanted sources of technical variation were regressed out using the `vars.to.regress` parameter. Highly variable genes (HVGs) were identified using `variable.features.n = 3000`, and the `method = "glmGamPoi"` option was used to improve computational efficiency during model fitting (Choudhary et al., 2023).

SCTransform normalization was executed on the HPC using an R script submitted via a shell job script.

### 4. Clustering


### 5. Annotation


### 6. Differential Expression


### 7. GSEA


## Results



## Discussion




## References

