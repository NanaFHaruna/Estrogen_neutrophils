# ---- Libraries ----
library(Seurat, lib.loc = "/data/Amblerwg_Schaughency/R4.3.0_lib")
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
library(patchwork)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)

# ---- Input ----
setwd("/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/")
seur <- readRDS("integrated_SeuratInt_gem_x_clustered.rds")
DefaultAssay(seur) <- "RNA"

# ---- Remove pellet sham samples ----
DimPlot(seur, label = TRUE)
Idents(seur) <- "orig.ident"
seur <- subset(seur, idents = c("Pellet_Blood_Sham", "Pellet_ckit_Sham", "Pellet_BM_Sham", "Pellet_Spleen_Sham"), invert = TRUE)

# ---- Restore compartment labels stored in the Sex metadata column ----
Idents(seur) <- "seurat_clusters"
DimPlot(seur, group.by = "Sex")

compartment_map <- data.frame(
  orig.ident = c(
    "CD451-WT-Veh-Blood", "CD451-WT-Veh-BM", "CD451-WT-Veh-ckit",
    "CD452-ERKO-Veh-Blood", "CD452-ERKO-Veh-BM", "CD452-ERKO-Veh-ckit"
  ),
  Sex = c("Blood", "BM", "ckit", "Blood", "BM", "ckit"),
  stringsAsFactors = FALSE
)

seur$Sex <- ifelse(
  is.na(seur$Sex),
  compartment_map$Sex[match(seur$orig.ident, compartment_map$orig.ident)],
  seur$Sex
)

DimPlot(seur, split.by = "Sex")
DimPlot(seur, split.by = "orig.ident", ncol = 4)
seur@meta.data$Compartment <- seur@meta.data$Sex

# ---- Split by compartment ----
Idents(seur) <- "Compartment"
blood_spleen <- subset(seur, idents = c("Blood", "Spleen"))
bm_ckit <- subset(seur, idents = c("ckit", "BM"))
Idents(seur) <- "seurat_clusters"
Idents(blood_spleen) <- "seurat_clusters"
Idents(bm_ckit) <- "seurat_clusters"

# ---- Blood/spleen cleanup ----
FeaturePlot(blood_spleen, features = "Pf4")
FeaturePlot(blood_spleen, features = "Hba-a1")
FeaturePlot(blood_spleen, features = "Irf8")

table(blood_spleen@meta.data$seurat_clusters)
blood_spleen_subset <- subset(blood_spleen, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,13,14,16,22,29,35,36))

blood_spleen_subset <- FindVariableFeatures(blood_spleen_subset)
blood_spleen_subset <- ScaleData(blood_spleen_subset)
blood_spleen_subset <- RunPCA(blood_spleen_subset)
blood_spleen_subset <- FindNeighbors(blood_spleen_subset, dims = 1:30)
blood_spleen_subset <- FindClusters(blood_spleen_subset, resolution = 0.4)
blood_spleen_subset <- RunUMAP(blood_spleen_subset, dims = 1:30)
DimPlot(blood_spleen_subset, reduction = "umap", label = TRUE)

# Check and remove the residual platelet cluster.
down_blood_spleen_subset <- subset(blood_spleen_subset, downsample = 5000)
cluster6.markers <- FindMarkers(down_blood_spleen_subset, ident.1 = 8)
blood_spleen_subset2 <- subset(blood_spleen_subset, idents = c(0,1,2,3,4,5,6,7,9))

# Recluster the cleaned subset.
sub <- blood_spleen_subset2
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- FindNeighbors(sub, reduction = "pca", dims = 1:30, graph.name = "sub_snn")
sub <- FindClusters(sub, graph.name = "sub_snn", resolution = 0.4)
DimPlot(sub, label = TRUE)

saveRDS(sub, "Blood_spleen_clean_gem_x2.rds")

# ---- BM/ckit ----
DimPlot(bm_ckit, label = TRUE)
saveRDS(bm_ckit, "BM_ckit_gem_x2.rds")
