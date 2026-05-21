suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(magrittr)
  library(knitr)
  library(scRNAseq)
  library(scater)
  library(stringr)
  library(reshape2)
  library(celldex)
  library(ggrepel)
  library(RColorBrewer)
  library(scCustomize)
  library(patchwork)
  library(Nebulosa)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(EnhancedVolcano)
})

# 1. Integrate the neutrophil subsets by batch, then cluster the integrated object.

source('AIS_shared_helpers.R')

project <- 'NIAMS-34'
workdir <- '/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Integrate/'
setwd(workdir)

seur.PBMC <- readRDS('/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge_Analysis/PBMC_seur_clusters.rds')
seur.neu <- subset(seur.PBMC, idents = c(0, 1, 10, 12))
saveRDS(seur.neu, 'seur_neu_clusters.rds')

save_png('UMAP_neu_SeuratInt_Pre.png', 10.7, 6.7, print(DimPlot(seur.neu, split.by = 'BatchID', ncol = 3)))

seur.list <- SplitObject(seur.neu, split.by = 'BatchID')
seur.list <- lapply(seur.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
  x
})

anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)
integrate.neu <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(integrate.neu) <- 'integrated'
saveRDS(integrate.neu, 'integrated_SeuratInt_neu.rds')

integrate.neu <- ScaleData(integrate.neu, verbose = FALSE)
integrate.neu <- RunPCA(integrate.neu, npcs = 30, verbose = FALSE)
integrate.neu <- RunUMAP(integrate.neu, reduction = 'pca', dims = 1:30)
integrate.neu <- FindNeighbors(integrate.neu, reduction = 'pca', dims = 1:30)
integrate.neu <- FindClusters(integrate.neu, resolution = 0.5)
saveRDS(integrate.neu, 'integrated_SeuratInt_neu_clustered.rds')

save_png('UMAP_neu_SeuratInt_Post.png', 10.7, 6.7, print(DimPlot(integrate.neu, split.by = 'BatchID', ncol = 3)))
write_session_info('Integrate_SeuratInt')
