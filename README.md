<img width="2048" height="768" alt="image" src="https://github.com/user-attachments/assets/0ce817aa-d7c5-44a3-bd54-a7a36d514ffb" />


**BIOME** is a Windows desktop application for automated end-to-end 16S and 18S rRNA amplicon sequencing analysis. It combines a Python-based graphical user interface with an integrated DADA2-based R pipeline, enabling users to process raw sequencing data through quality control, ASV inference, taxonomic classification, diversity analyses, differential abundance testing, and publication-ready visualizations.

The application is built on the **Automated Microbiome R Pipeline**, enabling researchers to perform microbiome analysis without writing R scripts or using the command line.
>
> **Automated-Microbiome-R-Pipeline**
>
> https://github.com/hoang-minh-tho/Automated-Microbiome-R-Pipeline.git

# Features

- Integrated primer removal using Cutadapt
- Automated DADA2-based ASV inference pipeline
- Automatic quality filtering and read trimming
- Advanced options for read processing and parameter optimization
- Automatic taxonomic assignment using SILVA (16S) and PR2 or SILVA (18S)
- Automatic rarefaction depth selection
- Alpha and beta diversity analyses
- Differential abundance analysis
- Publication-ready figures, and tables
- Supports both 16S and 18S rRNA amplicon sequencing

# Pipeline Overview

```
Raw FASTQ + Metadata
          │
          ▼
        BIOME
          │
          ▼
Automated Microbiome R Pipeline
          │
          ▼
       Outputs
   - Tables
   - Figures
   - .rds files
   - Log files

```

# Installation

Download the latest installer or portable version from the **Releases** section.

```
BIOME-Setup-v1.0.0.exe
```

or

```
BIOME-Portable-v1.0.0.zip
```

# Usage

1. Select the FASTQ folder.
2. Select the metadata file.
3. Choose an output directory.
4. Enter primer sequences.
5. Configure analysis parameters.
6. Click **Run selected steps**.

# Outputs

The pipeline generates:
- Alpha diversity tables and boxplots
- Beta diversity distance matrices
- PCoA, NMDS, or PCA figures
- PERMANOVA and betadisper results
- Differential abundance tables
- Relative abundance plots
- Rarefaction curves
- Publication-ready figures
- Analysis logs

- Tables
- Figures
- .rds files
- Log files

