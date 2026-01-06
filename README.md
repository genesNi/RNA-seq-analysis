# Bulk RNA_Seq_Analysis
## Beginner guide to learn how to use Bioinformatics tools and interpret results and understand Differential Gene Expression Analysis using DESeq2

Overview:

This project aims to show a complete Bulk RNA-Seq workflow from Raw sequencing data up to differential genes expression (DEGS) analysis using DESeq2. The analysis uses the `airway` dataset, a publicly available RNA-seq dataset bundled with Bioconductor.
for each human airway smooth muscle cell lines, Dexamethasone vs Untreated was compared to identify genes that responded to dexamethasone treatment.


# Data Sourcing
1. Download raw sequencing data (FASTA FILES) using SRA toolkit. Raw RNA-seq datasets submitted to the NCBI Sequence Read Archive are stored in a compressed binary format called .sra. To begin any bulk RNA-seq workflow, the first step is obtaining these files from SRA using tools such as prefetch.

2. Convert .sra → .fastq.gz using fastq-dump (or fasterq-dump) The .sra format is not directly usable by downstream RNA-seq tools such as FastQC or HISAT2.These tools require reads in FASTQ format, which contains: nucleotide sequences, per-base quality scores, identifiers. Therefore, each .sra file must be converted into .fastq (or better, .fastq.gz for compression). (WORKS GREAT ON ANY SYSTEM OTHER THAN MAC M1 M2 , BUT IF THEY WORK ON THESE THINGS LATER IT WILL BE HELPFUL , IN 2025 IT DIDN'T GOT FIXED)

```
sudo apt install sra-toolkit
prefetch SRR1039508

# Converting to fastq
   fastq-dump --outdir fastq --gzip --skip-technical --readids \
   --read-filter pass --dumpbase --split-3 --clip SRR1039508.sra
```
* loop for downloading many samples

3. From your sample copy the Bioproject ID and Open European Nucleotide Archive (ENA) website and paste it, you will get all the samples already present in fastq.gz format and can be simply downloaded using donwload sample tab or get download script tab which will give whole address which can be pasted to terminal and samples can be downloaded. (BEST OPTION FOR MAC M1 USERS)

* on a normal txt file just save all your urls and then use xargs 

```
cat airway_ena_url.txt | xargs -n 1 -P 3 wget -c
```
## Output :

```
(airway) MacBook-Pro:data megha$ ls	
SRR1039508_1.fastq.gz		SRR1039516_1.fastq.gz
SRR1039508_2.fastq.gz		SRR1039516_2.fastq.gz
SRR1039509_1.fastq.gz		SRR1039517_1.fastq.gz
SRR1039509_2.fastq.gz		SRR1039517_2.fastq.gz
SRR1039512_1.fastq.gz		SRR1039520_1.fastq.gz
SRR1039512_2.fastq.gz		SRR1039520_2.fastq.gz
SRR1039513_1.fastq.gz		SRR1039521_1.fastq.gz
SRR1039513_2.fastq.gz      SRR1039521_2.fastq.gz

```

4.Another tool which works the best is fastqdl which works way faster than fasterqdump but does not works on MAC M1's yet (https://youtu.be/3A-VrGAu7d4?si=9TcyBNVbHJqONynC )

# Quality Control Using tool Fastqc
we need to know the quality of reads whether or not it needs trimming. Reading multiple fastqc files is not possible so we first get all the fastqc files and then use another tool Multiqc to get a combined report about all the samples.

![fastqc](https://github.com/genesNi/RNA-seq-analysis/blob/main/fastqc.png)

![multiqc](<img width="1726" height="935" alt="multiqc" src="https://github.com/user-attachments/assets/cc32aa87-391c-4f1c-9580-ed8170dbd021" />)


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
