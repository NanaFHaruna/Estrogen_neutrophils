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
library(celldex, lib.loc = "/data/Amblerwg_Schaughency/R_libs")
library(ggrepel)
library(RColorBrewer)
library(scCustomize, lib.loc = "/data/Amblerwg_Schaughency/R4.3.0_lib")

rm(list = ls())

# ---- Project setup ----
project <- "NIAMS-34"
workdir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/"
setwd(workdir)
read.csv("Metadata.csv")

# ---- Read merged cluster objects ----
input_files <- c(
  SQ = "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x/SQ_1/SQ_seur_clusters.rds",
  Pellet = "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x/Pellet_2/Pellet_seur_clusters.rds",
  ERKO = "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x/ERKO_4/ERKO_seur_clusters.rds",
  Chimera2_E = "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x/Chimera2_E_5/Chimera2_E_seur_clusters.rds",
  Chimera2_Veh = "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x/Chimera2_Veh_6/Chimera2_Veh_seur_clusters.rds"
)
seur.list <- lapply(input_files, readRDS)

# ---- Integration ----
seur.list <- lapply(seur.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})

anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)
integrate.PBMC <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(integrate.PBMC) <- "integrated"
saveRDS(integrate.PBMC, "integrated_SeuratInt_gem_x.rds")

# ---- Clustering/UMAP on the integrated assay ----
integrate.PBMC <- ScaleData(integrate.PBMC, verbose = FALSE)
integrate.PBMC <- RunPCA(integrate.PBMC, npcs = 30, verbose = FALSE)
integrate.PBMC <- RunUMAP(integrate.PBMC, reduction = "pca", dims = 1:30)
integrate.PBMC <- FindNeighbors(integrate.PBMC, reduction = "pca", dims = 1:30)
integrate.PBMC <- FindClusters(integrate.PBMC, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

saveRDS(integrate.PBMC, "integrated_SeuratInt_gem_x_clustered.rds")

png("UMAP_gem_x_SeuratInt_Post.png", width = 1600, height = 1000, res = 150)
print(DimPlot(integrate.PBMC, split.by = "Batch_ID", ncol = 4))
dev.off()

writeLines(capture.output(sessionInfo()), "Integrate_SeuratInt_sessionInfo.txt")
