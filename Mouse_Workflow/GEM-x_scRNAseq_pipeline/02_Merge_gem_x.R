# ---- Libraries ----
library(Seurat)
library(magrittr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(SingleR)
library(scRNAseq)
library(scater)
library(knitr)
library(stringr)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(scCustomize)

rm(list = ls())

# ---- Project setup ----
project <- "NIAMS-34"
qcDir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate"
workdir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_gem_x/"
setwd(workdir)

metadata <- read.csv("Metadata.csv")
metadata$FullPath <- file.path(qcDir, metadata$Path)
cluster.list <- unique(metadata$Batch_ID)

save_cell_counts <- function(obj, cluster_name) {
  cellcounts <- data.frame(table(obj@meta.data$orig.ident))
  colnames(cellcounts) <- c("Sample", "Counts")

  png(paste0("CellCounts_", cluster_name, ".png"), width = 800, height = (nrow(cellcounts) + 1) * 100, res = 300)
  grid.table(cellcounts, rows = NULL)
  dev.off()
}

for (cluster.name in cluster.list) {
  cluster.paths <- metadata[metadata$Batch_ID == cluster.name, , drop = FALSE]
  seur_list <- lapply(cluster.paths$FullPath, readRDS)

  for (j in seq_along(seur_list)) {
    seur_list[[j]]$orig.ident <- cluster.paths$Sample.Name[j]
    names(seur_list)[j] <- cluster.paths$Sample.Name[j]
  }

  combined <- merge(
    seur_list[[1]],
    y = seur_list[-1],
    add.cell.ids = vapply(seur_list, function(x) unique(x$orig.ident), character(1))
  )

  save_cell_counts(combined, cluster.name)
  saveRDS(combined, paste0("merged_", cluster.name, ".rds"))

  rm(combined, seur_list)
}

writeLines(capture.output(sessionInfo()), "Merge_sessionInfo.txt")
