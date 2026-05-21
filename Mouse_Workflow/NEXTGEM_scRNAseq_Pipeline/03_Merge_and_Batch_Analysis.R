# ---- Load Libraries ----
# Loading essential Seurat, plotting, and data manipulation tools
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(scater)
library(SeuratData)

rm(list = ls())

## ---- Project Setup ----
# Defining the directory paths where corrected HTO files (Inputs) 
# and processed Batch files (Outputs) will live.
project <- "NIAMS-34"
qcDir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Correct_HTO_next_gem/"
workdir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_next_gem/"
setwd(workdir)

# Load master metadata file which serves as the "map" for the loop.
# It tells the script which samples belong to which batch (1, 2, 3...).
metadata_master <- read.csv("Metadata.csv")
batch_list <- unique(metadata_master$Batch_ID)

# ---- Start Combined Loop (Merging & Processing) ----
# This loop handles all batches 
# automatically while keeping the original internal logic intact.
for (b in batch_list) {
  message("Processing Batch: ", b)
  
  # 1. [Group Samples] Filter metadata for the specific Batch ID
  batch_info <- metadata_master[metadata_master$Batch_ID == b, ]
  batch_info$FullPath <- paste0(qcDir, batch_info$Path)
  
  # 2. [Load Data] Import individual corrected HTO objects into a list
  seur_list <- lapply(batch_info$FullPath, readRDS)
  names(seur_list) <- batch_info$Sample.Name
  
  # Ensure the internal identity matches the sample name for plotting later
  for (j in 1:length(seur_list)) {
    seur_list[[j]]$orig.ident <- names(seur_list)[j]
  }

  # 3. [Merge - Step 1] Combine samples into a single "Batch" object
  # add.cell.ids prepends the sample name to cell barcodes to prevent overlap
  combined <- merge(
    x = seur_list[], 
    y = seur_list[2:length(seur_list)], 
    add.cell.ids = names(seur_list)
  )
  
  # 4. [Processing - Step 2] Normalization and Variable Feature Detection
  # Identifies the 2000 most biologically informative genes
  DefaultAssay(combined) <- "RNA"
  combined <- NormalizeData(combined)
  combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  
  # Plotting top 10 most variable genes for QC
  top10 <- head(VariableFeatures(combined), 10)
  plot1 <- VariableFeaturePlot(combined)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  filename <- paste0("Batch_", b, "_features.png")
  png(filename, width=4500, height=2000, res = 300)
  print(plot1 + plot2)
  dev.off()
  
  # 5. [Analysis] Scaling, PCA, and Elbow Plot
  # Prepares data for dimensionality reduction by centering gene expression
  all.genes <- rownames(combined)
  combined <- ScaleData(combined, features = all.genes)
  combined <- RunPCA(combined, features = VariableFeatures(object = combined))
  
  # ElbowPlot helps decide how many Principal Components (PCs) are needed
  filename <- paste0("Batch_", b, "_elbow.png")
  png(filename, width=1800, height=1600, res = 300)
  print(ElbowPlot(combined, ndims=50))
  dev.off()
  
  # 6. [Visualization] Clustering & UMAP
  # Finds clusters based on PCA and creates a 2D map for viewing
  combined <- FindNeighbors(combined, dims = 1:50)
  combined <- FindClusters(combined, resolution = 0.8, algorithm=3, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = 'pca', dims = 1:50, assay = 'RNA')

  # 7. [Metadata Bridge] Attach Biological Traits
  # This block aligns the CSV metadata (Sex, Condition) with the Seurat object.
  # We subset the master metadata to match ONLY the samples currently being merged.
  metadata.PBMC <- batch_info[, c('Experiment', 'Sex', 'Condition', 'Chemistry', 'Sample.Name')]
  
  # Initialize the alignment data frame
  metadata.all <- data.frame(combined@meta.data[,c('orig.ident', 'Sex')])
  
  # Merge original metadata with CSV info based on Sample Name
  metadata.all <- merge(metadata.all, metadata.PBMC, by.x = 'orig.ident', by.y = 'Sample.Name', all.x = T, all.y = F, sort = F)
  metadata.all <- metadata.all[,-2] # Drop the dummy column created during initialization
  
  # Verify that the cell order hasn't changed during the merge
  if(sum((combined@meta.data$orig.ident) == metadata.all$orig.ident ) != nrow(combined@meta.data)) {
    warning("CRITICAL: Metadata alignment check failed for Batch ", b)
  }
  
  rownames(metadata.all) <- rownames(combined@meta.data)
  
  # Inject the metadata into the Seurat object slot
  combined <- AddMetaData(combined, metadata = metadata.all)
  
  # 8. [UMAP Outputs] Generate various standard plots for inspection
  ## 1: Standard UMAP
  filename <- paste0("Batch_", b, "_umap.png")
  png(filename, width=1800, height=1600, res = 300)
  print(DimPlot(combined, reduction = "umap", order = T) + ggtitle(''))
  dev.off()

  ## 2: Labeled Clusters
  filename <- paste0("Batch_", b, "_umap_labeled.png")
  png(filename, width=1800, height=1600, res = 300)
  print(DimPlot(combined, reduction = "umap", label = TRUE, order = T) + ggtitle(''))
  dev.off()

  ## 3: Colored by Sample
  filename <- paste0("Batch_", b, "_umap_sample.png")
  png(filename, width=1800, height=1600, res = 300)
  print(DimPlot(combined, reduction = "umap", group.by = "orig.ident", order = T))
  dev.off()

  ## 4: Split by Sample (to check for sample-specific clusters)
  filename <- paste0("Batch_", b, "_umap_sample_split.png")
  png(filename, width=1800, height=1600, res = 300)
  print(DimPlot(combined, reduction = "umap", split.by = "orig.ident", ncol = 3, order = T))
  dev.off()

  ## 9. [Final Output] Save the processed Batch object
  # These files are used as the input for the Integration script.
  saveRDS(combined, paste0("Batch_", b, "_seur_clusters.rds"))
  
  ## Cleanup RAM: Important for loops to avoid "Out of Memory" errors
  rm(combined, seur_list, metadata.all, metadata.PBMC)
  gc()
}

# Session Info for reproducibility
out_file <- paste0(project, "_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), out_file)
