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

# 1. Merge the sample-level QC objects by cluster and save one merged RDS per cluster.

source('AIS_shared_helpers.R')

project <- 'NIAMS-34'
qcDir <- '/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/SampleQC_noLowerBounds_Raw/'
workdir <- '/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge/'
setwd(workdir)

rds.path <- read.csv('RDSlist.csv', stringsAsFactors = FALSE)
rds.path$FullPath <- file.path(qcDir, rds.path$Path)

for (cluster.name in unique(rds.path$Cluster)) {
  cluster.paths <- rds.path[rds.path$Cluster == cluster.name, ]
  seur_list <- lapply(cluster.paths$FullPath, readRDS)

  for (j in seq_along(seur_list)) {
    seur_list[[j]]$orig.ident <- cluster.paths$Sample.Name[j]
  }

  if (length(seur_list) == 1L) {
    combined <- seur_list[[1]]
  } else {
    combined <- merge(
      x = seur_list[[1]],
      y = seur_list[-1],
      add.cell.ids = sapply(seur_list, function(x) unique(x$orig.ident))
    )
  }

  cellcounts <- as.data.frame(table(combined$orig.ident), stringsAsFactors = FALSE)
  names(cellcounts) <- c('Sample', 'Counts')

  save_png(paste0('CellCounts_', cluster.name, '.png'), 8, max(2, nrow(cellcounts) * 1), {
    grid.table(cellcounts, rows = NULL)
  })

  saveRDS(combined, paste0('merged_', cluster.name, '.rds'))
}

write_session_info('Merge')
