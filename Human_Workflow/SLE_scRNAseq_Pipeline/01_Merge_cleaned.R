# 1. Load libraries and prepare the merge workspace
library(Seurat)
library(dplyr)
library(magrittr)
library(gridExtra)
library(ggplot2)
library(reshape2)

rm(list = ls())

project <- "NIAMS-34"
qcDir <- "/data/Amblerwg_Schaughency/projects/0624_SLE/SampleQC_noLowerBounds_Raw/"
workdir <- "/data/Amblerwg_Schaughency/projects/0624_SLE/Merge/"
setwd(workdir)

# 2. Merge all QC-filtered Seurat objects listed in RDSlist.csv
rds_path <- read.csv("RDSlist.csv", stringsAsFactors = FALSE)
rds_path$FullPath <- file.path(qcDir, rds_path$Path)
cluster_list <- unique(rds_path$Cluster)

merge_cluster_group <- function(cluster_name) {
  cluster_paths <- rds_path[rds_path$Cluster == cluster_name, ]
  seur_list <- lapply(cluster_paths$FullPath, readRDS)

  for (j in seq_along(seur_list)) {
    seur_list[[j]]$orig.ident <- cluster_paths$Sample.Name[j]
    names(seur_list)[j] <- cluster_paths$Sample.Name[j]
  }

  combined <- merge(
    x = seur_list[[1]],
    y = seur_list[2:length(seur_list)],
    add.cell.ids = vapply(seur_list, function(x) unique(x$orig.ident), character(1))
  )

  cellcounts <- data.frame(table(combined@meta.data$orig.ident))
  colnames(cellcounts) <- c("Sample", "Counts")

  png(paste0("CellCounts_", cluster_name, ".png"), width = 800, height = (length(seur_list) + 1) * 100, res = 300)
  grid.table(cellcounts, rows = NULL)
  dev.off()

  saveRDS(combined, paste0("merged_", cluster_name, ".rds"))
}

for (cluster_name in cluster_list) {
  merge_cluster_group(cluster_name)
}

# 3. Save provenance information for reproducibility
writeLines(capture.output(sessionInfo()), "Merge_sessionInfo.txt")
