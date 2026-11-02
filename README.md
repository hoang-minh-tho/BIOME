# BIOME: Bioinformatic Interface for Optimized Multi-marker Ecology 

**BIOME** is a Windows desktop application for automated end-to-end 16S and 18S rRNA amplicon sequencing analysis. It combines a Python-based graphical user interface with an integrated DADA2-based R pipeline, enabling users to process raw sequencing data through quality control, ASV inference, taxonomic classification, diversity analyses, differential abundance testing, and publication-ready visualizations.

The application is built on the **Automated Microbiome R Pipeline**, enabling researchers to perform microbiome analysis without writing R scripts or using the command line.
>
> **Automated-Microbiome-R-Pipeline**
>
> https://github.com/grace-hoang/Automated-Microbiome-R-Pipeline

---

# Features

- Graphical desktop interface
- Automatic FASTQ file pairing
- Metadata validation
- Integrated primer removal using Cutadapt
- Automated DADA2 workflow
- Taxonomic assignment
- Alpha and beta diversity analysis
- Publication-ready figures and summary tables
- Windows installer and portable version

---

# Pipeline Overview

```
Raw FASTQ + Metadata
          │
          ▼
       Databio
          │
          ▼
Automated Microbiome R Pipeline
          │
          ▼
Results
 • ASV table
 • Taxonomy
 • Diversity analysis
 • Figures
 • Reports
```

---

# Installation

## Windows (Recommended)

Download the latest installer or portable version from the **Releases** page.

```
Databio-Setup-v1.0.0.exe
```

or

```
Databio-Portable-v1.0.0.zip
```
---

# Usage

1. Select the FASTQ folder.
2. Select the metadata file.
3. Choose an output directory.
4. Enter primer sequences if required.
5. Configure analysis parameters.
6. Click **Run Pipeline**.

---

# Outputs

The pipeline generates:

- ASV abundance table
- Taxonomic assignments
- Phyloseq object
- Alpha diversity analysis
- Beta diversity analysis
- Publication-ready figures
- Summary report
- Log files

---

# License

This project is distributed under the MIT License.

---

# Related Projects

**Automated Microbiome R Pipeline**

Standalone R implementation of the analysis workflow.

https://github.com/grace-hoang/Automated-Microbiome-R-Pipeline
