# scRNA-seq

## General Overview


## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
  - [1. Data Description](#1-data-description)
  - [2. Quality Control](#2-quality-control)
  - [3. ](#3-)
  - [4. ](#4-)
  - [5. (#5-)
  - [6. ](#6-)
  - [7. ](#7)
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


### 2. Quality Control


### 3.


## Results



## Discussion




## References

