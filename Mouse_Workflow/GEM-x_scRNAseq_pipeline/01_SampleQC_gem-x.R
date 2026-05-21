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
library(celldex)
library(ggrepel)
library(RColorBrewer)
library(scCustomize)

# ---- Project setup ----
project <- "NIAMS-34"
sample_list <- c(
  "CD451-WT-E2-Blood", "CD451-WT-Veh-Blood", "CD452-ERKO-E2-Blood", "CD452-ERKO-Veh-Blood",
  "ERKO_Blood_KO", "ERKO_Blood_WT", "Mrp8Esrfl_Blood_KO", "Mrp8Esrfl_Blood_WT",
  "Pellet_Blood_Estrogen", "Pellet_Blood_Oop", "Pellet_Blood_Sham",
  "SQ_Blood_Estrogen", "SQ_Blood_Oop", "SQ_Blood_Sham"
)
qcDir <- "/data/Amblerwg_Schaughency/will/Gonadectomy_final/Sample_QC_Blood/"
dataDir <- "/data/Amblerwg_Schaughency/rawdata/will/202503_mouse_blood"

# ---- Small helpers ----
create_sample_object <- function(data, project) {
  if (inherits(data, "dgCMatrix")) {
    CreateSeuratObject(counts = data, project = project)
  } else {
    CreateSeuratObject(counts = data$`Gene Expression`, project = project)
  }
}

plot_pre_filter <- function(seur, file_prefix) {
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^mt-")

  plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
  plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

  png(paste0(file_prefix, "PreFilter_Gene_Plot.png"), height = 5, width = 10, units = "in", res = 300)
  print(plot1 + plot2)
  dev.off()

  thresh <- list(
    nFeature_RNA_low  = median(log1p(seur$nFeature_RNA)) - 3 * mad(log1p(seur$nFeature_RNA)),
    nFeature_RNA_high = median(log1p(seur$nFeature_RNA)) + 3 * mad(log1p(seur$nFeature_RNA)),
    nCount_RNA_low    = median(log1p(seur$nCount_RNA)) - 3 * mad(log1p(seur$nCount_RNA)),
    nCount_RNA_high   = median(log1p(seur$nCount_RNA)) + 3 * mad(log1p(seur$nCount_RNA)),
    mt_high           = median(log1p(seur$percent.mito)) + 3 * mad(log1p(seur$percent.mito))
  )

  cellsToRemove.Feature_RNA <- colnames(seur)[which(log1p(seur$nFeature_RNA) < thresh$nFeature_RNA_low | log1p(seur$nFeature_RNA) > thresh$nFeature_RNA_high)]
  cellsToRemove.Count_RNA   <- colnames(seur)[which(log1p(seur$nCount_RNA) < thresh$nCount_RNA_low | log1p(seur$nCount_RNA) > thresh$nCount_RNA_high)]
  cellsToRemove.mito        <- colnames(seur)[which(log1p(seur$percent.mito) > thresh$mt_high)]
  thresh$numCellsRemove     <- length(unique(c(cellsToRemove.Feature_RNA, cellsToRemove.Count_RNA, cellsToRemove.mito)))

  plot1 <- VlnPlot(seur, features = "nFeature_RNA") + NoLegend() +
    geom_hline(yintercept = exp(thresh$nFeature_RNA_low) - 1, linetype = "dashed") +
    geom_hline(yintercept = exp(thresh$nFeature_RNA_high) - 1, linetype = "dashed")
  plot2 <- VlnPlot(seur, features = "nCount_RNA") + NoLegend() +
    geom_hline(yintercept = exp(thresh$nCount_RNA_low) - 1, linetype = "dashed") +
    geom_hline(yintercept = exp(thresh$nCount_RNA_high) - 1, linetype = "dashed")
  plot3 <- VlnPlot(seur, features = "percent.mito", ncol = 3) + NoLegend() +
    geom_hline(yintercept = exp(thresh$mt_high) - 1, linetype = "dashed")

  png(paste0(file_prefix, "PreFilter_VlnPlot_RNA.png"), height = 7, width = 7, units = "in", res = 300)
  grid.arrange(plot1, plot2, plot3, nrow = 1)
  dev.off()

  seur <- subset(
    seur,
    cells = unique(c(cellsToRemove.Feature_RNA, cellsToRemove.Count_RNA, cellsToRemove.mito)),
    invert = TRUE
  )

  thresh$nFeature_RNA_low  <- expm1(thresh$nFeature_RNA_low)
  thresh$nFeature_RNA_high <- expm1(thresh$nFeature_RNA_high)
  thresh$nCount_RNA_low    <- expm1(thresh$nCount_RNA_low)
  thresh$nCount_RNA_high   <- expm1(thresh$nCount_RNA_high)
  thresh$mt_high           <- expm1(thresh$mt_high)

  write.table(
    data.frame(Threshold = names(thresh), Value = unlist(thresh), row.names = NULL),
    "cell_filter_info.csv", quote = FALSE, row.names = FALSE, sep = ","
  )

  plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  png(paste0(file_prefix, "PostFilter_Gene_Plot.png"), height = 5, width = 10, units = "in", res = 300)
  print(plot1 + plot2)
  dev.off()

  png(paste0(file_prefix, "PostFilter_VlnPlot_RNA.png"), height = 7, width = 7, units = "in", res = 300)
  print(VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()

  seur
}

process_sample <- function(sample, project, dataDir, qcDir) {
  message("Processing sample: ", sample)

  data_path <- file.path(dataDir, sample, "count", "sample_filtered_feature_bc_matrix")
  data <- Read10X(data_path)
  seur <- create_sample_object(data, project)

  workdir <- file.path(qcDir, sample)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  setwd(workdir)

  seur <- plot_pre_filter(seur, file_prefix = "")

  # ---- RNA normalization and clustering ----
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)

  plot1 <- VariableFeaturePlot(seur)
  top10 <- head(VariableFeatures(seur), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  png("TopVariableGenes.png", height = 7, width = 12, units = "in", res = 300)
  print(plot1 + plot2)
  dev.off()

  seur <- ScaleData(seur, features = rownames(seur))
  seur <- RunPCA(seur)

  png("Heatmap_15PC.png", height = 12, width = 8, units = "in", res = 300)
  print(DimHeatmap(seur, dims = 1:15, cells = 500, balanced = TRUE))
  dev.off()

  png("ElbowPlot.png", height = 6, width = 8, units = "in", res = 300)
  print(ElbowPlot(seur, ndims = 50))
  dev.off()

  seur <- FindNeighbors(seur, dims = 1:30)
  seur <- FindClusters(seur, resolution = 0.8, algorithm = 3, verbose = FALSE)
  seur <- RunUMAP(seur, reduction = "pca", dims = 1:30, assay = "RNA")

  png("UMAP_RNA.png", width = 1800, height = 1600, res = 300)
  print(DimPlot(seur, reduction = "umap", label = TRUE) + ggtitle("RNA"))
  dev.off()

  saveRDS(seur, "seur_cluster.rds")
  writeLines(capture.output(sessionInfo()), paste0(project, "_", sample, "_QC_sessionInfo.txt"))

  invisible(seur)
}

# ---- Run all samples ----
for (sample in sample_list) {
  process_sample(sample, project, dataDir, qcDir)
}
