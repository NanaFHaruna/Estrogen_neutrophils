# ---- Load Libraries ----
# Loading the specific libraries required for integration and downstream analysis
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

## ---- Project setup ----
# Setting the directory for integration outputs
project <- "NIAMS-34"
workdir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_next_gem/"
setwd(workdir)

# Load metadata to reference sample/batch info
metadata <- read.csv("Metadata.csv")

## ---- Read Merged Objects ----
# Importing the batch-specific processed files generated in the previous step.
# These files have already undergone normalization and initial clustering.
One   <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/One/One_seur_clusters.rds')
Two   <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/Two/Two_seur_clusters.rds')
Three <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/Three/Three_seur_clusters.rds')
Four  <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/Four/Four_seur_clusters.rds')
Five  <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/Five/Five_seur_clusters.rds')

## ---- Integration Setup ----
# Combining the individual batch objects into a single list for processing
seur.list <- c(One, Two, Three, Four, Five) 

# Removing individual objects from memory to free up RAM for the integration steps
rm(One, Two, Three, Four, Five)

# Standardizing the list: Re-running Normalization and Variable Feature detection
# as per the original workflow to ensure consistency across the objects.
seur.list <- lapply(X = seur.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Finding 'anchors' which are the mathematical points of overlap between batches
anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:30)

# Creating the master integrated object 'integrate.PBMC'
integrate.PBMC <- IntegrateData(anchorset = anchors, dims = 1:30)

# Cleanup the list and anchors now that the integrated object is created
rm(seur.list, anchors)
gc() # Forced garbage collection to ensure stability

# Set the default assay to 'integrated' for the downstream dimensionality reduction
DefaultAssay(integrate.PBMC) <- "integrated"

# Initial save of the raw integrated object (checkpoint)
saveRDS(integrate.PBMC, 'integrated_SeuratInt_next_gem.rds')

# ---- Standard Workflow for Visualization and Clustering ----
# Scaling and PCA are performed on the 'integrated' assay values
integrate.PBMC <- ScaleData(integrate.PBMC, verbose = FALSE)
integrate.PBMC <- RunPCA(integrate.PBMC, npcs = 30, verbose = FALSE)

# UMAP Generation: Creating a 2D coordinate system for the integrated data
integrate.PBMC <- RunUMAP(integrate.PBMC, reduction = "pca", dims = 1:30)

# Graph-based Clustering: Identifying groups of similar cells across batches
integrate.PBMC <- FindNeighbors(integrate.PBMC, reduction = "pca", dims = 1:30)
integrate.PBMC <- FindClusters(integrate.PBMC, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

# Save the final clustered master object
saveRDS(integrate.PBMC, 'integrated_SeuratInt_next_gem_clustered.rds')

# ---- Visualization Check ----
# Generating UMAP plots to verify that batches are correctly integrated.
# This splits the UMAP by Batch_ID to check for batch-specific clusters.
png("UMAP_next_gem_SeuratInt_Post.png", width = 1600, height = 1000, res = 150)
print(DimPlot(integrate.PBMC, split.by = "Batch_ID", ncol = 4))
dev.off()

# Clean up memory
rm(integrate.PBMC)
gc()

# Export session info for reproducibility
out_file <- paste0("Integrate_SeuratInt_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), out_file)
