# Pipeline-for-scRNA-seq-and-Spatial-Transcriptomics-in-R

Single-Cell and Spatial Transcriptomics Analysis Pipeline
This pipeline was developed during the “Gaúcha School of Bioinformatics” in July 2025. 

All credits for Carolina Saibro-Girardi 

1. Overview

This repository provides a comprehensive R-based pipeline for the analysis of single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics data. The workflow starts from raw data processing and guides the user through quality control, normalization, clustering, differential expression analysis, and cell type annotation.

The pipeline is structured into modular R scripts, each handling a specific stage of the analysis:

    DOC1_environment_setup.R: Sets up the R environment and installs all required packages.

    DOC2_scRNAseq_processing.R: Processes a 10X PBMC dataset, from data loading to clustering.

    DOC3_spatial_processing.R: Processes a 10X Visium mouse brain dataset.

    DOC4_downstream_analysis.R: Performs differential expression and annotation on the processed PBMC data.


2. Getting Started
Prerequisites

    R (version 4.3.0)

    RStudio (recommended)

    Access to a command line/terminal.

Installation

    Clone the repository:

    git clone [URL_to_your_repository]
    cd [repository_name]

    Set up the R Environment:
    Open R or RStudio and run the first script to install all necessary packages from CRAN, Bioconductor, and GitHub.

    source("DOC1_environment_setup.R")

    Note: This step may take a significant amount of time, as it involves compiling and installing numerous packages.

3. How to Run the Pipeline

The pipeline is designed to be run sequentially. After setting up the environment, execute the scripts in order.
Step 1: Process Single-Cell RNA-seq Data

This script downloads the PBMC dataset, performs QC, normalization, dimensionality reduction, and clustering.

# Run from within R/RStudio
source("DOC2_scRNAseq_processing.R")

    Input: Raw 10X PBMC data (downloaded automatically).

    Output:

        PBMC.RDS: A Seurat object containing the processed data.

        results/PBMC_CellQC.pdf: Quality control plots.

        results/PBMC_PCA.pdf: PCA diagnostic plots.

        results/PBMC_clustering.pdf: UMAP/t-SNE plots of the final clusters.

Step 2: Process Spatial Transcriptomics Data (Optional)

This script analyzes a 10X Visium mouse brain dataset.

# Run from within R/RStudio
source("DOC3_spatial_processing.R")

    Input: Mouse brain data from the SeuratData package.

    Output:

        Brain.RDS: A Seurat object with processed spatial data.

        results/Brain_SpotQC.pdf, results/Brain_clustering.pdf, etc.

Step 3: Perform Downstream Analysis

This script loads the processed PBMC.RDS object and performs differential expression analysis and cell type annotation.

# Run from within R/RStudio
source("DOC4_downstream_analysis.R")

    Input:

        PBMC.RDS object from Step 1.

        data/PBMC_ref_CellTypeMarkers.xlsx for manual annotation.

    Output:

        top.markers_pbmc.xlsx: A table of the top marker genes for each cluster.

        irGSEA.AUCell_pbmc.xlsx: Cell type scores from AUCell analysis.

        Various plots in the results/ directory visualizing gene expression and cell annotations (Annotation_pbmc.pdf, Heatmap_top5markers_pbmc.pdf, etc.).

