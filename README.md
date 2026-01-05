# Bulk RNA_Seq_Analysis
## Beginner guide to learn how to use Bioinformatics tools and interpret results and understand Differential Gene Expression Analysis using DESeq2

Overview:

This project aims to show a complete Bulk RNA-Seq workflow from Raw sequencing data up to differential genes expression (DEGS) analysis using DESeq2. The analysis uses the `airway` dataset, a publicly available RNA-seq dataset bundled with Bioconductor.

For each human airway smooth muscle cell lines, Dexamethasone vs Untreated was compared to identify genes that responded to dexamethasone treatment.
---

## 📁 Files

- `DESeq2_SA_MEGHA.Rmd` – R Markdown file containing the full analysis pipeline.
- `README.md` – Project overview and usage instructions.

---

## 🧪 Analysis Overview

The workflow includes the following steps:

1. **Dataset Loading**
   - Uses the `airway` Bioconductor package to load RNA-seq data.
   - Samples are human airway smooth muscle cells treated with or without dexamethasone.

2. **Preprocessing**
   - Filtering low-expression genes (`sumofmin10 ≥ 2`)
   - Converting treatment conditions to factor variables

3. **Differential Expression Analysis**
   - Normalization using `DESeq2`
   - Statistical testing to identify differentially expressed genes
   - Filtering by adjusted p-value < 0.05

4. **Visualization**
   - PCA plots, MA plots, volcano plots
   - Boxplots and scatterplots of individual genes (e.g., `ENSG00000106211`)

5. **Functional Enrichment Analysis**
   - Gene Ontology (GO) enrichment (BP, MF, CC)
   - Visualizations using dot plots and bar plots

---
## ▶️ How to Run

1. Open `DESeq2_SA_MEGHA.Rmd` in **RStudio**
2. Make sure the following packages are installed:

```r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "airway", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("tidyverse", "pheatmap", "ggplot2"))
