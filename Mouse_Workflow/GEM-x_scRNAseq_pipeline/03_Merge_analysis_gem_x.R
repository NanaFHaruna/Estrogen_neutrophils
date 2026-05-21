# ---- Libraries ----
library(Seurat)
library(gridExtra)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(scater)
library(SeuratData)

rm(list = ls())

# ---- Project setup ----
project <- "NIAMS-34"
workdir <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge_analysis_gem_x"
setwd(workdir)

# ---- Load the merged blood object ----
seur <- readRDS("/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Merge/merged_Blood.rds")
DefaultAssay(seur) <- "RNA"

# ---- Standard RNA workflow ----
seur <- NormalizeData(seur)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(seur)
top10 <- head(VariableFeatures(seur), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(paste0(project, "_features.png"), width = 4500, height = 2000, res = 300)
print(plot1 + plot2)
dev.off()

seur <- ScaleData(seur, features = rownames(seur))
seur <- RunPCA(seur, features = VariableFeatures(object = seur))

png(paste0(project, "_elbow.png"), width = 1800, height = 1600, res = 300)
print(ElbowPlot(seur, ndims = 50))
dev.off()

seur <- FindNeighbors(seur, dims = 1:50)
seur <- FindClusters(seur, resolution = 0.8, algorithm = 3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = "pca", dims = 1:50, assay = "RNA")

saveRDS(seur, "Blood_seur_clusters.rds")

# ---- Add sample metadata ----
metadata <- read.csv("Metadata.csv")
metadata.PBMC <- metadata[metadata$Cluster2 == "Blood", c("Sample.Name", "Sex", "CellType", "Condition", "Chemistry", "Experiment", "Batch_ID")]

# Start with existing orig.ident/Sex columns and then merge in the sample annotations.
seur$Sex <- NA
metadata.all <- data.frame(seur@meta.data[, c("orig.ident", "Sex")])
metadata.all <- merge(metadata.all, metadata.PBMC, by.x = "orig.ident", by.y = "Sample.Name", all.x = TRUE, sort = FALSE)
metadata.all <- metadata.all[, -2]
rownames(metadata.all) <- rownames(seur@meta.data)
colnames(metadata.all)[3] <- "Sex"
seur <- AddMetaData(seur, metadata = metadata.all)

# ---- UMAP views ----
png(paste0(project, "_umap.png"), width = 1800, height = 1600, res = 300)
print(DimPlot(seur, reduction = "umap", order = TRUE) + ggtitle(""))
dev.off()

png(paste0(project, "_umap_labeled.png"), width = 1800, height = 1600, res = 300)
print(DimPlot(seur, reduction = "umap", label = TRUE, order = TRUE) + ggtitle(""))
dev.off()

png(paste0(project, "_umap_sample.png"), width = 1800, height = 1600, res = 300)
print(DimPlot(seur, reduction = "umap", group.by = "orig.ident", order = TRUE))
dev.off()

png(paste0(project, "_umap_sample_split.png"), width = 1800, height = 1600, res = 300)
print(DimPlot(seur, reduction = "umap", split.by = "orig.ident", ncol = 3, order = TRUE))
dev.off()

saveRDS(seur, "Blood_seur_clusters.rds")
writeLines(capture.output(sessionInfo()), paste0(project, "_sessionInfo.txt"))
