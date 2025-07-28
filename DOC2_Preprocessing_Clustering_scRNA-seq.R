### MC13: Análise de Dados de Trancriptômica Single-Cell
### Carolina Saibro-Girardi and Natanael Flach, 2025

### DOC2: loading seurat object, preprocessing, and clustering for scRNA-seq

# set the working directory

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(readxl)
library(openxlsx)
library(ggplot2) 
library(future) 
library(gridExtra)

#### CREATE SEURAT OBJECT FROM THE RAW DATA ####################################

# for this tutorial, we will be analyzing the a dataset of Peripheral Blood
# Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700
# single cells that were sequenced on the Illumina NextSeq 500.
# ref: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
file <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dest_file <- basename(file)
download.file(url = file, destfile = dest_file, mode = "wb")
extract_dir <- "pbmc3k"
untar(tarfile = dest_file, exdir = extract_dir)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = file.path("pbmc3k", "filtered_gene_bc_matrices", "hg19"))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
# Check object info
pbmc


#### CELL QUALITY CONTROL AND FILTERING ########################################

# Check initial object dimensions (features,samples)
dim(pbmc)

# Create annotation of percentage of mitochondrial counts per cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Cell QC parameters:
# percent.mt: The percentage of reads that map to the mitochondrial genome indicate low-quality / dying cells
# nFeature_RNA: Low-quality cells or empty droplets will often have very few genes
# nCount_RNA: Cell doublets or multiplets may exhibit an aberrantly high gene count

# View histogram of percent.mt
hist(pbmc$percent.mt, xlab = "percent.mt", main = "Histogram of MT genes")

# plot cell QC parameters
plot.cellQC = VlnPlot(pbmc,
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      pt.size = 0, ncol = 3)
# view QC parameters pre-filtering
plot.cellQC
# plot scatter plots linking distinct QC parameters
QCdotplot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = "none")
QCdotplot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = "none")
QCdotplot3 = FeatureScatter(pbmc, feature1 = "percent.mt", feature2 = "nFeature_RNA") + theme(legend.position = "none")

# apply filters
new_seurat <- subset(pbmc, subset = nFeature_RNA > 200 &
                       nFeature_RNA < 2000 &
                       nCount_RNA < 10000 &
                       percent.mt < 5) # 5-10% is a common threshold for % mt

# check new object dimensions 
dim(new_seurat)

# adjust filtering criteria according to plot.cellQC and object new dimensions

# repeat plot cell QC parameters
new_plot.cellQC = VlnPlot(new_seurat,
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          pt.size = 0, ncol = 3)

# save results
pdf("PBMC_CellQC.pdf")
plot.cellQC + labs(caption = "Raw object")
grid.arrange(QCdotplot1,
             QCdotplot2,
             QCdotplot3,
             nrow = 2)
new_plot.cellQC + labs(caption = "Filtered object")
dev.off()

# adopt the filtered new_seurat object for further analysis
pbmc <- new_seurat
rm(new_seurat)


#### PREPROCESSING #############################################################

## Normalization, feature selection and scaling

# Use sctransform for normalization, feature selection and scaling
pbmc <- SCTransform(pbmc, do.scale = T,
                    return.only.var.genes = F,
                    variable.features.n = 3000,
                    vst.flavor = 'v2')
# note that this command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
# transformed data will be available in the SCT assay, which is set as the default after running sctransform
# ref: https://satijalab.org/seurat/reference/sctransform

# Optional: regress possible confounding sources of variation
# regression can be perfomed within SCTransform argument "vars.to.regress"

## Dimensionality reduction

# Run linear dimensionality reduction with Principal Component Analysis (PCA)
pbmc <- RunPCA(pbmc, verbose = FALSE)
# plot PCA
ElbowPlot(pbmc)
DimHeatmap(pbmc, dims = 1:15, cells = 300, balanced = TRUE)
plot.PCAfeat = VizDimLoadings(pbmc, dims = 1:15, reduction = "pca")
plot.PCA1xPCA2 = DimPlot(pbmc, reduction = "pca")
plot.PCA2xPCA3 = DimPlot(pbmc, reduction = "pca", dims = 2:3)
plot.PCA1xPCA3 = DimPlot(pbmc, reduction = "pca", dims = c(1,3))
# select the top PCs, that will be continuously used in the next steps

# save PCA results
pdf("PBMC_PCA.pdf", width=9)
ElbowPlot(pbmc)
DimHeatmap(pbmc, dims = 1:15, cells = 300, balanced = TRUE)

for (p in plot.PCAfeat) {
  print(p)  
}
gridExtra::grid.arrange(plot.PCA1xPCA2, plot.PCA2xPCA3, plot.PCA1xPCA3, nrow = 2)
dev.off()

# Run UMAP and t-SNE for non-linear dimensionality reduction
# Uniform manifold approximation and projection (UMAP) highlight the differences between cell communities
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
# t-distributed stochastic neighbor embedding (t-SNE) highlight the similarities between cell communities
pbmc <- RunTSNE(pbmc, dims = 1:10, verbose = FALSE)

# Plot
DimPlot(pbmc, reduction = "umap") 
DimPlot(pbmc, reduction = "tsne")


#### CLUSTERING ###############################################################

# Run K-nearest neighbor (KNN) graph analysis
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# Apply Louvain algorithm to find clusters
# resolution parameter: sets the ‘granularity’ of the downstream clustering,
# with increased values leading to a greater number of clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# save object
saveRDS(pbmc, file = 'PBMC.RDS')


#### VISUALIZATION #############################################################

# Plot umap and tsne
umap <- DimPlot(pbmc, reduction = "umap")
tsne <- DimPlot(pbmc, reduction = "tsne")

# Plot distributions
data <- table(pbmc$seurat_clusters) # make df with distributions
df <- as.data.frame(data)
colnames(df) <- c('Cluster','Cells')
# total nr of cells
p1 <- ggplot(df, aes(x = "", y = Cells, fill = Cluster)) + geom_col()
p1
# % of cells
p2 <- ggplot(df, aes(x = "", y = Cells, fill = Cluster)) +
  geom_bar(position="fill", stat="identity")
p2
# save results
pdf("PBMC_clustering.pdf")
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")
grid.arrange(p1, p2,
             nrow = 1)
dev.off()


