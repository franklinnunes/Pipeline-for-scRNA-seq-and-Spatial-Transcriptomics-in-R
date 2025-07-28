### MC13: Análise de Dados de Trancriptômica Single-Cell
### Carolina Saibro-Girardi and Natanael Flach, 2025

### DOC3: loading seurat object, preprocessing, and clusterization for spatial RNA-seq

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(gridExtra)
library(ggplot2)
library(SeuratData)
library(future)

#### SPATIAL TRANSCRIPTOMICS: MOUSE BRAIN

# ref: https://satijalab.org/seurat/articles/spatial_vignette

#### CREATE SEURAT OBJECT FROM THE RAW DATA ####################################

# Here, we will be using a recently released dataset of sagital mouse brain slices
# generated using the Visium v1 chemistry. There are two serial anterior sections,
# and two (matched) serial posterior sections.
# ref: https://support.10xgenomics.com/spatial-gene-expression/datasets

# Import data using SeuratData
InstallData("stxBrain")
InstalledData()

# Load dataset for anterior and posterior slices into a seurat object
brain1 <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "posterior1")
# Check object info
brain1
brain2

# Merge objects
brain <- merge(brain1, brain2)
brain # Check object info


#### QUALITY CONTROL AND FILTERING ############################################

# Check initial object dimensions ([1] features samples)
dim(brain)

# Create annotation of percentage of mitochondrial counts per cell
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-") # For animal genes, use Mt- instead of MT-

# Plot nCount
vln1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0) + NoLegend()
anterior1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial", images = "anterior1")
posterior1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial", images = "posterior1")
wrap_plots(vln1, anterior1, posterior1)
# Plot nFeature
vln2 <- VlnPlot(brain, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
anterior2 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial", images = "anterior1")
posterior2 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial", images = "posterior1")
wrap_plots(vln2, anterior2, posterior2)
# Plot percent.mt
vln3 <- VlnPlot(brain, features = "percent.mt", pt.size = 0) + NoLegend()
anterior3 <- SpatialFeaturePlot(brain, features = "percent.mt", images = "anterior1")
posterior3 <- SpatialFeaturePlot(brain, features = "percent.mt", images = "posterior1")
wrap_plots(vln3, anterior3, posterior3)

# Filter out spots with aberrant QC metrics, which is usually unnecessary

# Save results
pdf("Brain_SpotQC.pdf", width = 8.5, height = 11)
grid.arrange(vln1,
             anterior1,
             posterior1,
             nrow = 2,
             top = "Counts per spot")
grid.arrange(vln2,
             anterior2,
             posterior2,
             nrow = 2,
             top = "Features per spot")
grid.arrange(vln3,
             anterior3,
             posterior3,
             nrow = 2,
             top = "Mt- percent per spot")
dev.off()


#### PREPROCESSING, DIM. REDUCTION, AND CLUSTERING #############################

brain <- SCTransform(brain, assay = "Spatial",
                     do.scale = T, return.only.var.genes = F, verbose = FALSE)
brain <- RunPCA(brain, verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:10)
brain <- RunTSNE(brain, dims = 1:10)
brain <- FindNeighbors(brain, dims = 1:10)
brain <- FindClusters(brain, verbose = FALSE)

# save object
saveRDS(brain, file = 'Brain.RDS')


#### VISUALIZATION #############################################################

# Plot umap and tsne
umap_brain <- DimPlot(brain, reduction = "umap")
tsne_brain <- DimPlot(brain, reduction = "tsne")

# Plot spatial
slice <- SpatialDimPlot(brain)

# Plot distributions
data <- table(brain$seurat_clusters, brain$orig.ident) # make df with distributions
df <- as.data.frame(data)
colnames(df) <- c('Cluster','Cells','Slice')
# total nr of cells
p1 <- ggplot(df, aes(fill=Cluster, y=Cells, x=Slice)) +
  geom_bar(position="stack", stat="identity")
# % of cells
p2 <- ggplot(df, aes(fill=Cluster, y=Cells, x=Slice)) +
  geom_bar(position="fill", stat="identity")

# save results
pdf("Brain_clustering.pdf")
umap_brain
tsne_brain
slice
grid.arrange(p1, p2,
             nrow = 1)
dev.off()

# Plot gene expression
genes <- c('Hpca','Nrgn',  # neuronal
           'Plp1',         # oligodendrocyte
           'Aif1',         # microglia
           'Gfap',         # astrocyte
           'Ttr',          # choroid plexus cell
           'Mgp')          # fibroblast
brain_genes <- list()
for (i in 1:length(unique(genes))) {
  brain_genes[[i]] <- SpatialFeaturePlot(brain, features = genes[i], ncol = 2)
}    

pdf("Brain_Genes.pdf")
for (i in brain_genes) {
  print(i)
  }
dev.off()

