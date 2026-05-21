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

# 1. Read each raw sample, apply the light QC filter, and run azimuth annotation.

source('AIS_shared_helpers.R')

project <- 'NIAMS-34'
sample_list <- paste0('GEXh_WA0', c(74, 81, 83, 71), '_outs')
qcDir <- '/data/Amblerwg_Schaughency/projects/0524_AIS/SampleQC_noLowerBounds_Raw/'
dataDir <- '/data/Amblerwg_Schaughency/rawdata/will/'

for (sample in sample_list) {
  message('Processing: ', sample)
  data_path <- file.path(dataDir, sample, 'raw_feature_bc_matrix')
  data <- Read10X(data_path)

  seur.raw <- if (inherits(data, 'dgCMatrix')) {
    CreateSeuratObject(counts = data, project = project)
  } else {
    CreateSeuratObject(counts = data[['Gene Expression']], project = project)
  }

  # Keep the same minimum-gene filter used in the original script.
  seur <- subset(seur.raw, subset = nFeature_RNA > 200)
  seur[['percent.mito']] <- PercentageFeatureSet(seur, pattern = '^MT-')
  rm(seur.raw)

  workdir <- file.path(qcDir, sample)
  make_dir(workdir)
  setwd(workdir)

  # 1a. QC plots after the initial filter.
  save_png('PostFilter_Gene_Plot.png', 10, 5, {
    print(FeatureScatter(seur, feature1 = 'nCount_RNA', feature2 = 'percent.mito') +
            FeatureScatter(seur, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA'))
  })

  save_png('PostFilter_VlnPlot_RNA.png', 7, 7, {
    print(VlnPlot(seur, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), ncol = 3))
  })

  # 1b. Standard Seurat workflow.
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
  save_png('TopVariableGenes.png', 12, 7, {
    plot1 <- VariableFeaturePlot(seur)
    top10 <- head(VariableFeatures(seur), 10)
    print(plot1 + LabelPoints(plot = plot1, points = top10, repel = TRUE))
  })

  seur <- ScaleData(seur, features = rownames(seur))
  seur <- RunPCA(seur)
  save_png('Heatmap_15PC.png', 8, 12, print(DimHeatmap(seur, dims = 1:15, cells = 500, balanced = TRUE)))
  save_png('ElbowPlot.png', 8, 6, print(ElbowPlot(seur, ndims = 50)))

  seur <- FindNeighbors(seur, dims = 1:30)
  seur <- FindClusters(seur, resolution = 0.8, algorithm = 3, verbose = FALSE)
  seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA')

  # 1c. Azimuth reference mapping and final QC outputs.
  seur <- RunAzimuth(seur, reference = '/data/OpenOmics/references/cyte-seek/hg38/Azimuth/pbmcref/')

  save_png('UMAP_RNA.png', 6, 5.3, print(DimPlot(seur, reduction = 'umap', label = TRUE) + ggtitle('RNA')))
  save_png('UMAP_azimuth.png', 6, 5.3, print(DimPlot(seur, reduction = 'umap', group.by = 'predicted.celltype.l1') + ggtitle('RNA')))

  saveRDS(seur, 'seur_cluster.rds')
  write_session_info(paste0(project, '_', sample, '_QC'))
}
