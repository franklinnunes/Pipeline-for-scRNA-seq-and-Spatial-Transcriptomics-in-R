### MC13: Análise de Dados de Trancriptômica Single-Cell
### Carolina Saibro-Girardi and Natanael Flach, 2025

### DOC1: Setting the environment 

# Copy course materials to your home dir
file.copy(from = list.files('/opt/EGBMC13', full.names = T), 
          to = ".")

###############################################################################

### Part 1: interact with the RStudio server

# Check the path for your working directory
getwd()
# Optional: create and set new working directories 
# new_dir <- "~/home/cgirardi"
# setwd(new_dir)

# Check files in your working directory
list.files()

# Check the path where your packages are saved
.libPaths()

###############################################################################

### Part 2: Install packages for single-cell and spatial RNA analysis 
# Execute RStudio in admin mode to avoid permission issues 

## From CRAN
install.packages('Seurat', version = '5.1.0')  # single-cell and spatial analysis, R v4.0 or greater is required
install.packages("sctransform")                # normalization and feature selection
install.packages("tidyverse")                  # dataframe editing
install.packages('dplyr')                      # dataframe editing
install.packages('ggplot2')                    # visualization
install.packages('patchwork')                    # visualization
install.packages('gridExtra')                    # visualization
install.packages('openxlsx')                   # interaction with xlsx file format
install.packages('devtools')                   # install packages from github
install.packages("BiocManager")           # install packages from bioconductor
install.packages("doMC", repos="http://R-Forge.R-project.org")    # dependency
# irGSEA dependencies from CRAN (levou 1 minuto)
cran.packages <- c("aplot", "circlize", "cowplot","data.table",
                   "doParallel", "doRNG", "ggfun", "gghalves",
                   "ggplotify", "ggridges", "ggsci", "irlba",
                   "magrittr", "Matrix", "msigdbr", "plyr", "pointr",
                   "purrr", "RcppML", "reshape2", "reticulate",
                   "rlang", "RMTstat", "RobustRankAggreg", "roxygen2",
                   "stringr")
for (i in cran.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i, ask = F, update = F)
  }
}

## From Bioconductor
BiocManager::install("glmGamPoi")   # optimization of single-cell and spatial analysis
BiocManager::install("MAST")    # DE testing for single-cell and spatial, using MAST (levou 4 minutos com atualização de todos pacotes)
# irGSEA dependencies from CRAN 
bioconductor.packages <- c('AUCell', 'BiocParallel', 'ComplexHeatmap', 'decoupleR',
                           'fgsea', 'ggtree', 'GSEABase', 'Nebulosa', 'scde',
                           'SummarizedExperiment', 'sparseMatrixStats', 'GSVA', 'singscore', 'UCell')
for (i in bioconductor.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    BiocManager::install(i, ask = F, update = F)
  }
}

## From GitHub 
devtools::install_github('satijalab/seurat-data')    # import 10X data from repository
devtools::install_github("https://github.com/JiaLiVUMC/scMRMA")  # automated annotation of clusters 
devtools::install_github("https://github.com/chuiqin/irGSEA")   # single cell rank-based gene set enrichment analysis based on AUCell 
