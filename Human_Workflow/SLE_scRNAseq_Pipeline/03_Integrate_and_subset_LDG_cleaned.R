# 1. Load libraries and set the integration workspace
library(Seurat)
library(dplyr)
library(magrittr)
library(gridExtra)
library(ggplot2)
library(SingleR)
library(scRNAseq)
library(scater)
library(knitr)
library(stringr)
library(reshape2)
library(celldex)
library(ggrepel)
library(RColorBrewer)
library(scCustomize)

rm(list = ls())

project <- "NIAMS-34"
workdir <- "/data/Amblerwg_Schaughency/projects/0624_SLE/Integrate/"
setwd(workdir)
metadata_file <- "../RDSlist_metadata_batch.csv"

save_plot <- function(filename, width, height, plot_expr) {
  png(filename, width = width, height = height, res = 150)
  print(plot_expr)
  dev.off()
}

add_batch_metadata <- function(seur, metadata_df, cluster_name) {
  metadata_sub <- metadata_df[metadata_df$Cluster == cluster_name, c("Sample.Name", "BatchID")]
  seur$BatchID <- NA
  metadata_all <- data.frame(seur@meta.data[, c("orig.ident", "BatchID")])
  metadata_all <- merge(metadata_all, metadata_sub, by.x = "orig.ident", by.y = "Sample.Name", all.x = TRUE, all.y = FALSE, sort = FALSE)
  metadata_all <- metadata_all[, -2]
  rownames(metadata_all) <- rownames(seur@meta.data)
  colnames(metadata_all)[2] <- "BatchID"
  AddMetaData(seur, metadata = metadata_all)
}

integrate_batch_split <- function(seur, output_pre, output_post, plot_prefix) {
  seur.list <- SplitObject(seur, split.by = "BatchID")
  seur.list <- lapply(seur.list, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x
  })
  anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)
  integrate.obj <- IntegrateData(anchorset = anchors, dims = 1:30)
  DefaultAssay(integrate.obj) <- "integrated"
  saveRDS(integrate.obj, output_pre)

  integrate.obj <- ScaleData(integrate.obj, verbose = FALSE)
  integrate.obj <- RunPCA(integrate.obj, npcs = 30, verbose = FALSE)
  integrate.obj <- RunUMAP(integrate.obj, reduction = "pca", dims = 1:30)
  integrate.obj <- FindNeighbors(integrate.obj, reduction = "pca", dims = 1:30)
  integrate.obj <- FindClusters(integrate.obj, resolution = 0.8)
  saveRDS(integrate.obj, output_post)
  save_plot(paste0(plot_prefix, "_Post.png"), 1600, 1000, DimPlot(integrate.obj, split.by = "BatchID", ncol = 3))
  integrate.obj
}

# 2. Integrate the neutrophil and PBMC objects by batch
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

seur.neutrophil <- readRDS('../Merged_neutrophils_analysis/neutrophils_seur_clusters.rds')
seur.neutrophil <- add_batch_metadata(seur.neutrophil, metadata, "neutrophil")
saveRDS(seur.neutrophil, '../Merged_neutrophils_analysis/neutrophils_seur_clusters.rds')
save_plot('UMAP_neutrophils_SeuratInt_Pre.png', 1600, 1000, DimPlot(seur.neutrophil, split.by = "BatchID", ncol = 3))
integrate.neutrophil <- integrate_batch_split(seur.neutrophil, 'integrated_SeuratInt_neutrophil.rds', 'integrated_SeuratInt_neutrophil_clustered.rds', 'UMAP_neutrophils_SeuratInt')
rm(seur.neutrophil, integrate.neutrophil)

seur.PBMC <- readRDS('../Merged_PBMC_analysis/PBMC_seur_clusters.rds')
seur.PBMC <- add_batch_metadata(seur.PBMC, metadata, 'PBMC')
saveRDS(seur.PBMC, '../Merged_PBMC_analysis/PBMC_seur_clusters.rds')
save_plot('UMAP_PBMC_SeuratInt_Pre.png', 1600, 1000, DimPlot(seur.PBMC, split.by = 'BatchID', ncol = 3))
integrate.PBMC <- integrate_batch_split(seur.PBMC, 'integrated_SeuratInt_PBMC.rds', 'integrated_SeuratInt_PBMC_clustered.rds', 'UMAP_PBMC_SeuratInt')
rm(seur.PBMC, integrate.PBMC)

# 3. Subset LDG cells from the integrated PBMC object and recluster that branch
pbmc.integrated <- readRDS('integrated_SeuratInt_PBMC_clustered.rds')
DimPlot(pbmc.integrated, label = TRUE)
FeaturePlot(pbmc.integrated, features = 'FCGR3B')
FeaturePlot(pbmc.integrated, features = 'S100A9')
FeaturePlot(pbmc.integrated, features = 'S100A8')
FeaturePlot(pbmc.integrated, features = 'LYZ')
FeaturePlot(pbmc.integrated, features = 'CD14')
FeaturePlot(pbmc.integrated, features = 'CSF3R')
FeaturePlot(pbmc.integrated, features = 'LCN2')

neutro1 <- subset(pbmc.integrated, idents = c(20, 34))
DimPlot(neutro1, label = TRUE)
FindMarkers(neutro1, ident.1 = 34)
FeaturePlot(neutro1, features = 'MPO')

neutro1 <- FindVariableFeatures(neutro1, selection.method = 'vst', nfeatures = 2000)
neutro1 <- ScaleData(neutro1, features = rownames(neutro1))
neutro1 <- RunPCA(neutro1, features = VariableFeatures(object = neutro1))
neutro1 <- FindNeighbors(neutro1, dims = 1:50)
neutro1 <- FindClusters(neutro1, resolution = 0.5, algorithm = 3, verbose = FALSE)
neutro1 <- RunUMAP(neutro1, reduction = 'pca', dims = 1:50, assay = 'RNA')
DimPlot(neutro1, label = TRUE)

cluster8.markers <- FindMarkers(neutro1, ident.1 = 8)
cluster4.markers <- FindMarkers(neutro1, ident.1 = 4)
cluster7.markers <- FindMarkers(neutro1, ident.1 = 7)
neutro2 <- subset(neutro1, idents = c(0, 1, 2, 3, 4, 5, 7))
saveRDS(neutro2, file = 'LDG_0624_will.rds')

# 4. Integrate the neutrophil and LDG objects into the final neutrophil total object
will <- readRDS('LDG_0624_will.rds')
neutro <- readRDS('integrated_SeuratInt_neutrophil_clustered.rds')
DefaultAssay(will) <- 'RNA'
DefaultAssay(neutro) <- 'RNA'

seur.list <- c(will, SplitObject(neutro, split.by = 'orig.ident'))
seur.list <- lapply(seur.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
  x
})

anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)
integrate.neutrophil <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(integrate.neutrophil) <- 'integrated'
saveRDS(integrate.neutrophil, 'integrated_SeuratInt_neutrophil_total.rds')

integrate.neutrophil <- ScaleData(integrate.neutrophil, verbose = FALSE)
integrate.neutrophil <- RunPCA(integrate.neutrophil, npcs = 30, verbose = FALSE)
integrate.neutrophil <- RunUMAP(integrate.neutrophil, reduction = 'pca', dims = 1:30)
integrate.neutrophil <- FindNeighbors(integrate.neutrophil, reduction = 'pca', dims = 1:30)
integrate.neutrophil <- FindClusters(integrate.neutrophil, resolution = 0.8)
saveRDS(integrate.neutrophil, 'integrated_SeuratInt_neutrophil_clustered_total.rds')
save_plot('UMAP_neutrophils_SeuratInt_Post_total.png', 1600, 1000, DimPlot(integrate.neutrophil, split.by = 'BatchID', ncol = 3))
