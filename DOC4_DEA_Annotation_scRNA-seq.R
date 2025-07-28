### MC13: Análise de Dados de Trancriptômica Single-Cell
### Carolina Saibro-Girardi and Natanael Flach, 2025

### DOC4: downstream analysis of single-cell datasets (post cluster analysis)

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(openxlsx)
library(RobustRankAggreg)
library(ggplot2)

# set the working directory

# load the seurat object
pbmc <- readRDS("~/PBMC.RDS")
pbmc


#### DIFFERENTIAL EXPRESSION TESTING ###########################################

# Use MAST as test for DE testing
# ref: https://rglab.github.io/MAST/
library(MAST)

### Generate table with top markers for each cluster

# Define thresholds additional to adj.p-val<0.05
FC <- 0.5      # limit testing to genes which show at least X-fold difference (log-scale) between the two groups
PCT <- 0.25    # only test genes detected in a minimum fraction of cells in either of the two populations
# Run FindAllMarkers
markers <- FindAllMarkers(pbmc,
                          min.cells.group = 1,
                          min.cells.feature = 1,
                          logfc.threshold = FC,
                          min.pct = PCT,
                          only.pos = T,
                          test.use = 'MAST')
# save table
openxlsx::write.xlsx(markers, file = 'top.markers_pbmc.xlsx', rowNames=T)

### Visualization of gene markers

# Select top markers
markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.3) %>%
  top_n(n = 5, wt = avg_log2FC) -> top.markers # select markers of interest

# Heatmap with top markers
# dsample <- subset(pbmc, downsample = 500) # downsample the data the make the plot easier, especially for larger datasets
DoHeatmap(pbmc, features = top.markers$gene) + NoLegend()
ggsave('Heatmap_top5markers_pbmc.pdf')

# Select one marker for each cluster
markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.5) %>%
  top_n(n = 1, wt = avg_log2FC) -> top.markers # select markers of interest
# Violin plot (stacked)
VlnPlot(pbmc, features = top.markers$gene, stack=TRUE, fill.by="ident", flip=TRUE)
ggsave('Vlnplot_topmarker_pbmc.pdf')
# Dot plot (stacked)
DotPlot(pbmc, features = top.markers$gene) +
  coord_flip() + theme(axis.text.x = element_text(angle = 90))
ggsave('Dotplot_topmarker_pbmc.pdf')


### Optional: Generate table with DE for all clusters, including all genes
# this can be a time-consuming task on larger datasets

# Define thresholds additional to adj.p-val<0.05
FC <- 0     # limit testing to genes which show at least X-fold difference (log-scale) between the two groups
PCT <- 0    # only test genes detected in a minimum fraction of cells in either of the two populations
# Run FindMarkers to compare ident.1 vs ident.2 (all)
de_all <- list()
for (i in 1:length(unique(pbmc$seurat_clusters))) {
  de_all[[i]] <- FindMarkers(pbmc,
                             ident.1 = (i-1),
                             ident.2 = NULL,
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             logfc.threshold = FC,
                             min.pct = PCT,
                             only.pos = F,
                             test.use = 'MAST')
}
# add suffix to colnames for each cluster
a <- 0:(length(unique(pbmc$seurat_clusters))-1)
for (i in 1:length(a)) {
  colnames(de_all[[i]]) <- paste(colnames(de_all[[i]]), "cluster", a[i], sep="_")
}
# merge all DE into one table
DEtable <- de_all[[1]]
for (i in 2:length(unique(pbmc$seurat_clusters))) {
  DEtable <- merge(DEtable, de_all[[i]], by = 0, all = T)
  rownames(DEtable) <- DEtable[, 1]
  DEtable <- DEtable[,-1]
}

#save table
openxlsx::write.xlsx(DEtable, file = 'DEtable.complete_pbmc.xlsx', rowNames=T)


#### ANNOTATION ################################################################

#### Automated annotation

# Automated annotation of clusters using single cell Multi-Resolution Marker-based Annotation Algorithm (scMRMA) and Panglao database
# ref: https://github.com/JiaLiVUMC/scMRMA 
library(scMRMA)

autom.ann <- scMRMA(input = pbmc,
                    species = "Hs",
                    db = "panglaodb",
                    p = 0.05,
                    normalizedData = F,
                    selfDB = NULL,
                    selfClusters = NULL,
                    k=20)

# salve as metadata in seurat object
pbmc[["scMRMA"]] <- autom.ann$multiR$annotationResult[colnames(pbmc), ncol(autom.ann$multiR$annotationResult)]

# save umap
pdf("scMRMA_pbmc.pdf")
DimPlot(pbmc,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE)
dev.off()

#### Manual annotation

# Apply user-defined lists of reference marker genes to identify cell types

# load the list of reference marker genes
df <- read_excel('/opt/EGBMC13/PBMC_ref_CellTypeMarkers.xlsx')  # import from EGBMC13 dir in the server
ref_markers <- list(Monocyte = df[which(df$celltype == 'Monocyte'),]$gene,
                    Dendritic = df[which(df$celltype == 'Dendritic'),]$gene,
                    NK = df[which(df$celltype == 'NK'),]$gene,
                    T.CD8 = df[which(df$celltype == 'T-CD8'),]$gene,
                    B.cell = df[which(df$celltype == 'B-cell'),]$gene,
                    Naive.T.CD4 = df[which(df$celltype == 'Naive-T-CD4'),]$gene,
                    Mem.T.CD4 = df[which(df$celltype == 'Mem-T-CD4'),]$gene,
                    Platelet = df[which(df$celltype == 'Platelet'),]$gene)

# Use a scoring tool to assign a signature based on the list of marker genes

# irGSEA to use AUCell or other scoring tools
# ref: https://github.com/chuiqin/irGSEA
# install
devtools::install_github("https://github.com/chuiqin/irGSEA", upgrade = "never")   
library(irGSEA) # load

# AUCell
# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#follow-up-examples
# ref: https://doi.org/10.1038/nmeth.4463
library(AUCell)

# Compute scores
seurat_scored <- irGSEA.score(object = pbmc,
                              custom = T, geneset = ref_markers, msigdb = F,
                              method = 'AUCell',
                              aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                              kcdf = 'Gaussian')
# Compute differential gene sets using robust rank aggregation algorithm (RRA) and Wilcoxon test
result.dge <- irGSEA.integrate(seurat_scored,
                               group.by = "seurat_clusters",
                               method = c('AUCell'))
# Save table
openxlsx::write.xlsx(result.dge, file="irGSEA.AUCell_pbmc.xlsx")

# Plot scores
# set list of genes from signature
ct <- c(names(ref_markers)) # underscore is not accepted
# plot scatter
scatter <- irGSEA.density.scatterplot(object = seurat_scored,
                                      method = c('AUCell'),
                                      show.geneset = ct,
                                      reduction = "umap")
# plot violins
library(scales)
hue <- hue_pal()(length(unique(pbmc$seurat_clusters))) # set colors according to nr. of clusters (same colors as seurat base)
halfvln <- list()
for (i in 1:length(unique(ct))) {
  halfvln[[i]] <- irGSEA.halfvlnplot(object = seurat_scored, 
                                     method = "AUCell", 
                                     show.geneset = ct[i], 
                                     group.by = "seurat_clusters", 
                                     color.cluster = hue)
}
# Save plots
pdf("AUCell_pbmc.pdf", widt)
print(scatter)
for (i in halfvln) {
  print(i)
}
dev.off()


#### Save annotation as metadata

# Create new celltype metadata 
pbmc[["celltype"]] <- NA
# Assign celltype annotation
pbmc@meta.data[which(pbmc$seurat_clusters %in% c(1,6)),]$celltype <- "Monocyte"
pbmc@meta.data[which(pbmc$seurat_clusters == 0),]$celltype <- "Naive T CD4"
pbmc@meta.data[which(pbmc$seurat_clusters == 2),]$celltype <- "Memory T CD4"
pbmc@meta.data[which(pbmc$seurat_clusters == 3),]$celltype <- "B cell"
pbmc@meta.data[which(pbmc$seurat_clusters == 4),]$celltype <- "T CD8"
pbmc@meta.data[which(pbmc$seurat_clusters == 5),]$celltype <- "NK"
pbmc@meta.data[which(pbmc$seurat_clusters == 7),]$celltype <- "Dendritic"
pbmc@meta.data[which(pbmc$seurat_clusters == 8),]$celltype <- "Platelet"
# Check the number of cells per celltype
table(pbmc$celltype)

# Plot annotation
ct_pie <- pie(table(pbmc$celltype), col = hue)
ct_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "celltype") + NoLegend()
ct_umap

# Save plots
pdf("Annotation_pbmc.pdf")
ct_pie <- pie(table(pbmc$celltype), col = hue)
ct_umap
dev.off()


#### VISUALIZE GENE EXPRESSION ################################################

## Plot the expression level of selective genes

# Set gene list
gene_list <- c('CD14', 'CD8A', 'MS4A1','ACTB')

# dotplot
genes_dotplot <- DotPlot(pbmc, features = gene_list) + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90))
# vlnplot
genes_vlnplot <- VlnPlot(pbmc, features = gene_list,
                         stack=TRUE, fill.by="ident", flip=TRUE) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))
# umap
genes_umap <- FeaturePlot(pbmc, reduction = "umap", features = gene_list) 

# Save plots
pdf("GeneListExpression_pbmc.pdf")
gridExtra::grid.arrange(genes_dotplot, genes_vlnplot, nrow=2)
genes_umap
dev.off()
