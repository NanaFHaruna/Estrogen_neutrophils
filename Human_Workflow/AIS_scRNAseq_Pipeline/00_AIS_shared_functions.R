# Shared utilities for the AIS workflow.
# Source this file at the top of the downstream analysis scripts.

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

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_png <- function(filename, width, height, expr, res = 300) {
  png(filename, width = width, height = height, units = 'in', res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(filename)
}

write_session_info <- function(prefix) {
  writeLines(capture.output(sessionInfo()), paste0(prefix, '_sessionInfo.txt'))
  invisible(TRUE)
}

rename_clusters <- function(object, new_labels) {
  names(new_labels) <- levels(object)
  RenameIdents(object, new_labels)
}

vln_box <- function(object, feature, group.by, pt.size = 0, ylim = NULL) {
  p <- VlnPlot(object, features = feature, group.by = group.by, pt.size = pt.size) +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

run_deg <- function(object, id_col, downsample_n, ident.1 = 'HC', ident.2 = 'AIS', test.use = 'MAST') {
  Idents(object) <- id_col
  if (!is.null(downsample_n)) {
    object <- subset(object, downsample = downsample_n)
  }
  deg <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, test.use = test.use)
  deg
}

deg_to_df <- function(deg) {
  deg <- as.data.frame(deg)
  deg$gene <- rownames(deg)
  rownames(deg) <- NULL
  deg[, c('gene', setdiff(names(deg), 'gene'))]
}

read_deg_for_volcano <- function(path) {
  deg <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ('X' %in% names(deg)) {
    names(deg)[names(deg) == 'X'] <- 'gene'
  }
  deg
}

make_volcano_prep <- function(deg, cutoff = 5) {
  deg <- deg_to_df(deg)
  deg$gene <- if ('gene' %in% names(deg)) deg$gene else rownames(deg)
  rownames(deg) <- deg$gene
  up <- rownames(subset(deg, avg_log2FC >= cutoff))
  down <- rownames(subset(deg, avg_log2FC <= -cutoff))
  capped <- deg
  capped$avg_log2FC <- pmax(pmin(capped$avg_log2FC, cutoff), -cutoff)

  shape <- rep(19, nrow(capped))
  names(shape) <- rep('group1', nrow(capped))
  shape[rownames(capped) %in% up] <- -9658
  names(shape)[rownames(capped) %in% up] <- 'group2'
  shape[rownames(capped) %in% down] <- -9668
  names(shape)[rownames(capped) %in% down] <- 'group3'

  size <- rep(2.0, nrow(capped))
  size[rownames(capped) %in% up] <- 8
  size[rownames(capped) %in% down] <- 8

  list(raw = deg, capped = capped, up = up, down = down, shape = shape, size = size)
}
