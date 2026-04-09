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
The dataset was explored prior to quality control to examine metadata and assess the number of cells and genes. Quality control was performed in R (v4.5.2) by calculating the percentage of mitochondrial gene expression (`percent.mt`) as an indicator of cellular stress. Cells undergoing stress or apoptosis, which can occur during viral infections such as Influenza A virus (IAV), often exhibit elevated mitochondrial RNA content due to leakage of cytoplasmic RNA. Cells with mitochondrial percentages exceeding 10–20% were considered low quality and potentially non-viable (Osorio & Cai, 2020).

To visualize the distribution of quality control metrics prior to filtering, violin plots were generated for each sample (grouped by `orig.ident`), displaying the number of detected genes per cell (`nFeature_RNA`), total RNA counts per cell (`nCount_RNA`), and the percentage of mitochondrial gene expression (`percent.mt`).

Additionally, scatter plots were used to assess relationships between quality control metrics. A strong positive correlation was observed between `nCount_RNA` and `nFeature_RNA` (Pearson correlation = 0.83), indicating that cells with higher RNA counts tend to have more detected genes. In contrast, a weak correlation was observed between mitochondrial percentage and gene count (Pearson correlation = -0.06), suggesting that mitochondrial content is largely independent of gene complexity.

Based on these metrics, cells were filtered to retain high-quality cells with `percent.mt < 10, nFeature_RNA > 200, and nFeature_RNA < 4000`. A new Seurat object containing only filtered cells was then created for downstream analysis.

Quality control analyses were performed in R(v4.5.2), and the corresponding [code](code/master_script.R) is provided in the project repository.

### 3. Normalization & Scaling
Prior to normalization, the filtered Seurat object was transferred to a high-performance computing (HPC) environment for efficient processing. Normalization was performed in R (4.5.0) using the `SCTransform` method implemented in the Seurat package. SCTransform performs normalization, scaling, and variance stabilization in a single step using a regularized negative binomial regression model, which accounts for sequencing depth and technical variability (Hafemeister & Satija, 2019).

During normalization, unwanted sources of technical variation were regressed out using the `vars.to.regress` parameter. Highly variable genes (HVGs) were identified using `variable.features.n = 3000`, and the `method = "glmGamPoi"` option was used to improve computational efficiency during model fitting (Choudhary et al., 2023).

SCTransform normalization was executed on the HPC using an [R Script](code/normalize/run_normalization.R) submitted via a [shell job script](code/normalize/normalize_job.sh).

### 4. Clustering
Following SCTransform normalization, the dataset was subjected to Principal Component Analysis (PCA) using the `RunPCA` function in the Seurat package. To determine the number of principal components (PCs) to retain, an elbow plot was generated for the first 50 PCs. Based on the point at which the standard deviation plateaued, the first 30 PCs were selected for downstream analysis.

PCA embeddings were also examined to assess potential batch effects across samples. As cells from different time points did not form distinct, segregated clusters, batch correction methods such as Harmony were not applied.

To identify distinct cell populations, a shared nearest neighbor (SNN) graph was constructed using the `FindNeighbors` function based on the selected 30 PCs. Graph-based clustering was then performed using the `Louvain` algorithm implemented in the `FindClusters` function. To evaluate clustering resolution, two resolutions (0.6 and 0.9) were tested, consistent with parameters reported in the original study (Kazer et al., 2024). Based on visual inspection of cluster separation and biological interpretability, a resolution of 0.6 was selected for downstream analysis.

For visualization and assessment of cluster structure, UMAP was performed using the same 30 PCs. UMAP plots were generated and colored by original sample identity (`orig.ident`) to assess potential batch effects, and by cluster identity to evaluate the biological organization of the data.

All analyses were executed on a high-performance computing (HPC) environment using [R Script](code/cluster/UMAP.R) submitted via [SLURM](code/cluster/final_UMAP.sh) job scheduling.

### 5. Annotation
To assign biological identities to clusters generated at a resolution of `0.6`, differential expression analysis was first performed using the `FindAllMarkers` function in the Seurat package. A Wilcoxon rank-sum test was used to identify positive marker genes, with a minimum `log₂ fold-change threshold of 0.25` and a requirement that genes be expressed in at least 25% of cells within each cluster. Marker genes were ranked by `average log₂ fold-change`, and the top markers per cluster were selected for downstream annotation.

Clusters were manually annotated based on the expression of canonical marker genes identified through differential expression analysis. Multiple markers per cluster were considered to ensure robust and accurate cell-type assignment. Both immune (e.g., Cd3d for T cells, Cd79a for B cells) and epithelial (e.g., Epcam) markers were used to distinguish major cell populations. Identified genes were validated using publicly available databases such as **PanglaoDB** as well as relevant literature (e.g., Durante et al., 2020).

Cell-type identities were assigned by cross-referencing cluster-specific markers with known tissue-specific expression profiles. For example, Adcy3 was used to identify mature olfactory neurons, Krt15 for horizontal basal cells, and C1qa for macrophages. Clusters exhibiting similar transcriptional profiles were grouped under shared cell-type identities. For instance, clusters 0, 2, 4, and 7 were collectively annotated as olfactory sensory neurons (OSNs) based on the expression of markers such as Kirrel3 and S100a5. Additional populations, including immune cells (e.g., neutrophils, B cells, NK cells) and structural cell types (e.g., endothelial cells and fibroblasts), were similarly identified.

Feature plots were generated using Seurat’s `FeaturePlot` function, to visualize the expression of key marker genes across clusters and to validate annotation decisions. The final annotated cell populations were visualized using UMAP to assess their distribution and separation in low-dimensional space.

Annotation was executed on a high-performance computing (HPC) environment using [R Script](code/annotation/annotate.R) submitted via [SLURM](code/annotation/annotate.sh) job scheduling.


### 6. Differential Expression
Differential expression analysis was performed on macrophage populations. The processed Seurat object was loaded, and the RNA assay was used for downstream analysis. To enable condition-specific comparisons, cells were grouped into two conditions: Naive and Infected, based on the `orig.ident` metadata field. The condition label was assigned such that all non-naive samples were classified as infected. Cell identities were then set to this condition variable for differential testing.

Differential gene expression between infected and naive macrophages was then computed using the Wilcoxon rank-sum test implemented in Seurat’s `FindMarkers` function. To ensure a robust comparison, genes were required to be expressed in at least 10% of cells in either group and exhibit a minimum `log₂ fold-change threshold of 0.25` for initial testing. P-values were adjusted for multiple testing using Bonferroni correction. Genes were defined as significantly differentially expressed if they reached an `adjusted p-value < 0.05` and an absolute `log₂ fold-change > 0.5`. Genes with `positive log₂ fold change` were considered upregulated in infected cells, while those with `negative log₂ fold change` were considered downregulated.

Volcano plots were generated using `ggplot2`, with genes colored based on significance category. The top differentially expressed genes were annotated using `ggrepel` to improve label readability.

To visualize gene expression patterns across cells, feature plots were generated, with cells split by condition (Naive vs Infected). Additionally, violin plots (`VlnPlot`) were used to examine the distribution of selected differentially expressed genes across individual time points. 

All differential expression analyses were executed on a high-performance computing (HPC) environment using [R Scripts](code/differential_expression/de.R) submitted via a [shell job script](code/differential_expression/de.sh).

### 7. GSEA
Functional interpretation of differentially expressed genes (DEGs) within the macrophage cluster was performed using Gene Set Enrichment Analysis (GSEA), enabling the identification of coordinated biological processes that may not be apparent from single-gene analysis.

Gene symbols obtained from differential expression analysis were converted to `Entrez IDs` using the `bitr` function and the `org.Mm.eg.db` (v3.16.0) mouse genome annotation database. To prepare input for GSEA, all genes were ranked in descending order based on their `average log₂ fold-change`, such that genes most upregulated in the Infected condition were positioned at the top of the ranked list, while those downregulated (or relatively higher in the Naive condition) were positioned at the bottom.

Enrichment analysis was conducted using the `clusterProfiler` package (v4.6.0) in R (v4.5.2). Two primary databases were queried: Gene Ontology (GO) Biological Process, to identify broad cellular functions and physiological responses, and KEGG pathways, to identify specific metabolic and signaling pathways.

For both analyses, gene sets with sizes between 15 and 500 were retained (`minGSSize = 15`, `maxGSSize = 500`) to exclude overly narrow or excessively broad categories. Statistical significance was determined using a p-value cutoff of `0.05` following Benjamini–Hochberg (BH) correction for multiple testing.

Results were visualized using dot plots generated with the `enrichplot` package. To distinguish between activated and suppressed pathways, results were stratified by the sign of the normalized enrichment score (NES) (`split = ".sign"`), where a positive NES indicates pathways enriched in infected macrophages and a negative NES indicates pathways enriched in the naive condition.

All analyses were performed in R, and the corresponding [code](code/gsea/gsea.R) is provided in the project repository.


## Results
### Quality Control
![Figure 1](results/quality_control/QC_violin_prefilter.png)

**Figure 1:** Distribution of Quality Control Metrics Prior to Filtering. This visualization displays the number of unique genes (Gene Count), total RNA counts (RNA Count), and mitochondrial gene percentage (Mitochondrial %) for each sample before data cleaning. The wide distribution of mitochondrial reads and the presence of cells with very low gene counts highlight the necessity of filtering to remove non-viable cells and technical noise.

![Figure 2](results/quality_control/QC_violin_postfilter.png)

**Figure 2:** Distribution of Quality Control Metrics Following Filtering. The metrics are shown after applying thresholds of 500–4,000 unique genes, a threshold of 0-6500 of total RNA counts, and a mitochondrial maximum of 10%. By narrowing these distributions, the dataset is restricted to high-quality, viable cells, providing a consistent baseline for the subsequent normalization.

Initial assessment of the scRNA-seq dataset revealed substantial variability in cellular quality across samples, as reflected in pre-filtering metrics ([Figure 1](results/quality_control/QC_violin_prefilter.png)). Violin plots of quality control features showed that most cells contained between 1,000 and 3,000 detected genes (Gene Count) and exhibited consistent total RNA counts (RNA Count) across time points, indicating uniform sequencing depth. However, a subset of cells displayed elevated mitochondrial gene expression (Mitochondrial %), with values extending up to ~15%, suggesting the presence of stressed or dying cells, likely associated with Influenza A infection or tissue processing.

To prevent these low-quality cells from biasing downstream clustering and functional analyses, a stringent filtering strategy was applied. Following filtering ([Figure 2](results/quality_control/QC_violin_postfilter.png)), the dataset was restricted to high-quality cells with mitochondrial content below 10% and a gene detection range between 500 and 4,000 unique features. This ensured that subsequent analyses were performed on biologically meaningful and high-confidence cell populations.

### Annotation
![Figure 3](results/annotation/UMAP_clusters_vs_annotations.png)

**Figure 3:** Comparative UMAP Analysis: Clustering vs. Biological Annotation. The side-by-side UMAP embeddings illustrate the unsupervised clustering at a resolution of 0.6 (left) and the subsequent manual cell-type annotation (right). It identifies 39 distinct cell populations including various neuronal, epithelial, and immune lineages.

Unsupervised clustering and dimensionality reduction revealed a highly complex cellular architecture within the respiratory tissues, represented by 39 distinct clusters ([Figure 3](results/annotation/UMAP_clusters_vs_annotations.png)). The global data structure was primarily driven by biological lineage, as evidenced by the clear spatial separation between neuronal populations (e.g., mature olfactory neurons, neuronal progenitors), epithelial cells (e.g., horizontal basal cells, secretory cells), and infiltrating immune populations.

Comparison of the cluster-based UMAP with the annotated version confirmed that a resolution of 0.6 effectively captured the biological diversity of the tissue without over-clustering technical noise. This map provided the foundation for lineage-specific analyses, enabling the precise isolation of the macrophage compartment.

Macrophages, comprising distinct clusters of homeostatic, M2-like, and tissue-resident populations, were identified across clusters 5, 13, and 16 and selected for further investigation of the transcriptomic response to Influenza A infection.

![Figure 4](results/annotation/featureplot_markers.png)

**Figure 4:** Distribution and intensity of specific gene transcripts. The color gradient indicates the expression level, with darker or highlighted regions representing higher concentration of mRNA for a given marker. Genes such as C1qa, Adcy3, and Krt15 are shown to be restricted to specific spatial coordinates on the UMAP, corresponding to the Macrophage, Mature Olfactory Neuron, and Basal Cell populations, respectively.

To validate the biological accuracy of the 39 identified clusters, feature plots were generated to visualize the expression of canonical lineage-specific marker genes ([Figure 4](results/annotation/featureplot_markers.png)). The results demonstrated clear spatial localization of marker expression within their respective annotated clusters.

For example, markers such as C1qa (macrophages), Adcy3 (neurons), and Krt15 (basal cells) exhibited highly restricted expression patterns, appearing in specific clusters rather than being broadly distributed across the dataset. This distinct localization of gene expression supports the biological validity of the clustering and confirms that each cluster represents a unique cell population with a characteristic transcriptomic profile.

### Differential Expression
![Figure 5](results/differential_expression/volcano_Macrophages.png)

**Figure 5:** Differential Gene Expression in Macrophages. This volcano plot illustrates the transcriptional changes in macrophages when comparing Infected vs. Naive conditions. The x-axis represents the log2 fold change, while the y-axis shows the statistical significance (-log10 adjusted p-value). Red points represent significantly upregulated genes, including Gm42418 (the most significant) and Obp1a (the highest fold change). Blue points indicate downregulated genes, such as Cd14, Nfkbia, and Cxcl2. Gray points represent genes that did not meet the significance threshold.

The differential expression analysis of the macrophage compartment revealed a distinct transcriptomic response to infection ([Figure 5](results/differential_expression/volcano_Macrophages.png)). The most statistically significant upregulated gene was Gm42418, while Obp1a showed the most substantial positive fold change, suggesting these may be key markers of the infection state in this tissue. Additionally, several genes associated with inflammatory signaling and structural regulation, such as mt-Atp8 and Tubb4b, were significantly induced. In contrast, the downregulation of several myeloid and signaling genes, including Cd14, Nfkbia, and Dusp1, was observed. 

### GSEA

## Discussion




## References

