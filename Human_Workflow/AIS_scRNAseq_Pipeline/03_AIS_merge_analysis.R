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

# 1. Recluster the merged PBMC object, add metadata, and generate the main QC summary plots.

source('AIS_shared_helpers.R')

project <- 'NIAMS-34'
workdir <- '/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge_Analysis/'
setwd(workdir)

seur <- readRDS('/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge/merged_PBMC.rds')
DefaultAssay(seur) <- 'RNA'

# 1a. Standard normalization, variable feature selection, PCA, clustering, and UMAP.
seur <- NormalizeData(seur)
seur <- FindVariableFeatures(seur, selection.method = 'vst', nfeatures = 2000)

save_png(paste0(project, '_features.png'), 15, 6.7, {
  plot1 <- VariableFeaturePlot(seur)
  plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(seur), 10), repel = TRUE)
  print(plot1 + plot2)
})

seur <- ScaleData(seur, features = rownames(seur))
seur <- RunPCA(seur, features = VariableFeatures(object = seur))

save_png(paste0(project, '_elbow.png'), 6, 6, print(ElbowPlot(seur, ndims = 50)))
seur <- FindNeighbors(seur, dims = 1:50)
seur <- FindClusters(seur, resolution = 0.8, algorithm = 3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:50, assay = 'RNA')
saveRDS(seur, 'PBMC_seur_clusters.rds')

# 1b. Attach sample-level metadata by matching orig.ident to Sample.Name.
metadata <- read.csv('RDSlist_metadata.csv', stringsAsFactors = FALSE)
metadata_pbmc <- metadata[metadata$Cluster == 'PBMC', c('Sample.Name', 'Patient.ID', 'Sex', 'Age', 'Identifier', 'Identifier2', 'BatchID')]
metadata_pbmc <- metadata_pbmc[match(seur$orig.ident, metadata_pbmc$Sample.Name), ]
rownames(metadata_pbmc) <- colnames(seur)
seur <- AddMetaData(seur, metadata = metadata_pbmc[, c('Sex', 'Patient.ID', 'Age', 'Identifier', 'Identifier2', 'BatchID')])

# 1c. UMAPs and cell-count summaries for manual review.
save_png(paste0(project, '_umap.png'), 6, 5.3, print(DimPlot(seur, reduction = 'umap', order = TRUE) + ggtitle('')))
save_png(paste0(project, '_umap_labeled.png'), 6, 5.3, print(DimPlot(seur, reduction = 'umap', label = TRUE, order = TRUE) + ggtitle('')))
save_png(paste0(project, '_umap_sample.png'), 6, 5.3, print(DimPlot(seur, reduction = 'umap', group.by = 'orig.ident', order = TRUE)))
save_png(paste0(project, '_umap_sample_split.png'), 6, 5.3, print(DimPlot(seur, reduction = 'umap', split.by = 'orig.ident', ncol = 3, order = TRUE)))
save_png(paste0(project, '_umap_Celltype.png'), 6, 5.3, print(DimPlot(seur, reduction = 'umap', group.by = 'predicted.celltype.l1', order = TRUE)))

save_png(paste0(project, '_seurat_sample_counts.png'), 8, 20, {
  seurat_clusters_table <- data.frame(table(seur$seurat_clusters), stringsAsFactors = FALSE)
  names(seurat_clusters_table) <- c('cluster', 'count')
  grid.table(seurat_clusters_table, rows = NULL, theme = ttheme_default(base_size = 10))
})

# 1d. Module score and marker gene summaries used downstream.
gene_list <- list(c('IFI27', 'IFI6', 'IFI44', 'IFI44L', 'USP18', 'LY6E', 'OAS1',
                    'ISG15', 'IFIT1', 'OAS3', 'HERC5', 'MX1', 'LAMP3', 'EPSTI1',
                    'IFIT3', 'OAS2', 'RTP4', 'PLSCR1', 'SPATS2L', 'RSAD2', 'SIGLEC1'))

seur <- AddModuleScore(seur, features = gene_list, name = 'IFN_scores1')
save_png(paste0(project, '_moduleScore.png'), 6, 5.3, print(FeaturePlot(seur, features = 'IFN_scores11')))
save_png(paste0(project, '_moduleScore_split.png'), 6, 5.3, print(FeaturePlot(seur, features = 'IFN_scores11', split.by = 'Identifier2')))
save_png(paste0(project, '_IFN_scores_violin.png'), 13.3, 13.3, print(VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = 0)))
save_png(paste0(project, '_IFN_scores_split_violin.png'), 16.7, 13.3, print(VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = 0, split.by = 'Identifier2')))

seur$combined <- 'Sample'
Idents(seur) <- 'combined'
save_png(paste0(project, '_IFN_scores_split_violin_noCluster.png'), 6.7, 13.3, print(VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = NULL, split.by = 'Identifier2')))
Idents(seur) <- 'seurat_clusters'

save_png(paste0(project, '_FCGR3BScore.png'), 6, 5.3, print(FeaturePlot(seur, features = 'FCGR3B')))
save_png(paste0(project, '_FCGR3BScore_split.png'), 6, 5.3, print(FeaturePlot(seur, features = 'FCGR3B', split.by = 'Identifier2')))
save_png(paste0(project, '_FCGR3BScore_split_Vln.png'), 10, 3.3, print(VlnPlot(seur, features = 'FCGR3B', split.by = 'Identifier2', pt.size = 0)))

saveRDS(seur, 'PBMC_seur_clusters.rds')
write_session_info(project)
