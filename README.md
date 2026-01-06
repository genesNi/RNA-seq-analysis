# Bulk RNA_Seq_Analysis
## Beginner guide to learn how to use Bioinformatics tools and interpret results and understand Differential Gene Expression Analysis using DESeq2

Overview:

This project aims to show a complete Bulk RNA-Seq workflow from Raw sequencing data up to differential genes expression (DEGS) analysis using DESeq2. The analysis uses the `airway` dataset, a publicly available RNA-seq dataset bundled with Bioconductor.
for each human airway smooth muscle cell lines, Dexamethasone vs Untreated was compared to identify genes that responded to dexamethasone treatment.
GSE_ID: GSE52778.

# Operating System :
to be used for preprocessing : Linux
But starting with Windows would require WSL (Windows Subsytem for linux )
or Linux as Operating System in it (can be done by using Virtual Box which can run both windows and linux operating system within the same machine)
Youtube link: (https://youtu.be/QBfh0RkCYOI?si=PXRNkvAVrQXs97qr) (for wsl)
             ( https://youtu.be/YjG1yG2l9v0?si=82E519JZXckRBAQl) (use of Virtual box to download linux on windows)
 * on your laptop you can choose which operating system you want as Bioinformatician Linux is a great choice but as a beginner it would be difficult to deal with

# Why RAM is extremely important :
Sequenced files are very large so your datasets will take up a lot space of your RAM. 
if you are starting out with 4 or 8 GB RAM you should work with very small datasets. 
for this dataset 8 GB RAM would work fine but do remember to keep your samples compressed , but sometimes some tools might require uncompressed form then work in sections do not work with all samples together.
if you have option start out with 16 GB RAM 

# Data Sourcing
Opt 1. Download raw sequencing data (FASTA FILES) using SRA toolkit. Raw RNA-seq datasets submitted to the NCBI Sequence Read Archive are stored in a compressed binary format called .sra. To begin any bulk RNA-seq workflow, the first step is obtaining these files from SRA using tools such as prefetch.

Convert .sra → .fastq.gz using fastq-dump (or fasterq-dump) The .sra format is not directly usable by downstream RNA-seq tools such as FastQC or HISAT2.These tools require reads in FASTQ format, which contains: nucleotide sequences, per-base quality scores, identifiers. Therefore, each .sra file must be converted into .fastq (or better, .fastq.gz for compression). (WORKS GREAT ON ANY SYSTEM OTHER THAN MAC M1 M2 , BUT IF THEY WORK ON THESE THINGS LATER IT WILL BE HELPFUL , IN 2025 IT DIDN'T GOT FIXED)

```
sudo apt install sra-toolkit
prefetch SRR1039508

# Converting to fastq
   fastq-dump --outdir fastq --gzip --skip-technical --readids \
   --read-filter pass --dumpbase --split-3 --clip SRR1039508.sra
```
* loop for downloading many samples

Opt 2. From your sample copy the Bioproject ID and Open European Nucleotide Archive (ENA) website and paste it, you will get all the samples already present in fastq.gz format and can be simply downloaded using donwload sample tab or get download script tab which will give whole address which can be pasted to terminal and samples can be downloaded. (BEST OPTION FOR MAC M1 USERS)

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

Opt 3.Another tool which works the best is fastqdl which works way faster than fasterqdump but does not works on MAC M1's yet (https://youtu.be/3A-VrGAu7d4?si=9TcyBNVbHJqONynC )

# Quality Control Using tool Fastqc
we need to know the quality of reads whether or not it needs trimming.Quality control is essential before alignment because it allows detection of: 1.Poor-quality reads 2.Adapter contamination 3.Over-represented sequences 4.GC content deviations 5.Per-base quality drops.

Reading multiple fastqc files is not possible so we first get all the fastqc files and then use another tool Multiqc to get a combined report about all the samples. MultiQC aggregates all individual FastQC reports into a single unified report, making it easier to compare quality metrics across multiple samples in one place. This helps quickly identify sample-specific issues, batch effects, or global sequencing problems without opening dozens of separate HTML files.

Undersatnding Box Pot ()
Undersatnding Normal Distribution ()

![fastqc](https://github.com/genesNi/RNA-seq-analysis/blob/main/fastqc.png)

![multiqc](https://github.com/genesNi/RNA-seq-analysis/blob/main/multiqc.png)

for understanding your fastqc reports better (https://sequencing.qcfail.com )

# Read Trimming: 
Multiple tools we have  : Trimmomatic , Cutadapt, Fastp(launched in 2019) choosing a tool requires remembering your RAM, how fast it can process an output for you and how easily it can be downloaded and used.
These tools have function of trimming reads either from head or trail , remove adapters, remove any poor quality base etc 

On the basis of reading your Multiqc Report you can decide whether or not trimming is required.
Few things which can be noticed:


  * Trimmomatic Manual: (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
  * Cutadapt Manual: (https://cutadapt.readthedocs.io/en/stable/)
  * Fastp Manual: (https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp)


 
 # What is Technical Replicates:

 
 if your datasets has multiple technical replicates for each sample then you might would require to concatinate those separate files, as done in (https://github.com/tiashaghosh024-dotcom/Bulk-RNA-seq-project)

# Alignment 
After checking the quality of reads we need to align it to reference Genome or Transcriptome.
Multiple tools we have but are categories on the basis of what it aligns to, Genome (STAR, HISAT2) or Transcriptome : Aligners (Bowtie2, BWA) , Quasi-mappers (Salmon, Kallisto)

+ Remember to use tools according to the technology used whether Reads are short read or Long
+ how fast tool works and what minimum RAM it takes
+ what is the need of your experiment

##### Reference Genome have types: read (https://www.nature.com/articles/s41592-025-02850-9)
Databases for downloading reference Genomes are many, watch (https://youtu.be/eIVlSG11umQ?si=RVUrztO8cmSVrskn)
  
All of these tools require Genome Indexes either you can build index on your terminal or download Genome indexes from their websites according to the organism used in your dataset
  
* STAR : (https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
* HISAT2: (https://daehwankimlab.github.io/hisat2/manual/)
* Bowtie2: (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* Salmon: (https://app.readthedocs.org/projects/salmon/downloads/pdf/stable/)

# Quantification
After Alignment we need to count which genes expressed what number of reads.
we have multiple tools for that : htseq-counts, feature-counts, RSEM


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
