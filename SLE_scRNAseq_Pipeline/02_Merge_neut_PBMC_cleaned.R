# 1. Load libraries and define shared helpers
library(Seurat)
library(dplyr)
library(magrittr)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(reshape2)
library(stringr)

rm(list = ls())

project <- "NIAMS-34"
root_dir <- "/data/Amblerwg_Schaughency/projects/0624_SLE"
metadata_csv <- file.path(root_dir, "RDSlist_metadata.csv")

save_plot <- function(filename, width, height, plot_expr) {
  png(filename, width = width, height = height, res = 300)
  print(plot_expr)
  dev.off()
}

save_table_png <- function(tbl, filename, width, height, base_size = 10) {
  png(filename, width = width, height = height, res = 300)
  grid.table(tbl, theme = ttheme_default(base_size = base_size), rows = NULL)
  dev.off()
}

add_sample_metadata <- function(seur, metadata_file, cluster_name) {
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  metadata_sub <- metadata[metadata$Cluster == cluster_name, c("Sample.Name", "Patient.ID", "Sex", "Age", "Identifier")]

  seur$Sex <- NA
  metadata_all <- data.frame(seur@meta.data[, c("orig.ident", "Sex")])
  metadata_all <- merge(metadata_all, metadata_sub, by.x = "orig.ident", by.y = "Sample.Name", all.x = TRUE, all.y = FALSE, sort = FALSE)
  metadata_all <- metadata_all[, -2]

  rownames(metadata_all) <- rownames(seur@meta.data)
  colnames(metadata_all)[3] <- "Sex"
  AddMetaData(seur, metadata = metadata_all)
}

run_merged_dataset <- function(input_rds, workdir, cluster_name, output_rds) {
  setwd(workdir)

  # 2. Load the merged object and run the shared Seurat QC workflow
  seur <- readRDS(input_rds)
  DefaultAssay(seur) <- "RNA"
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)

  top10 <- head(VariableFeatures(seur), 10)
  plot1 <- VariableFeaturePlot(seur)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  save_plot(paste0(project, "_features.png"), 4500, 2000, plot1 + plot2)

  all.genes <- rownames(seur)
  seur <- ScaleData(seur, features = all.genes)
  seur <- RunPCA(seur, features = VariableFeatures(object = seur))
  save_plot(paste0(project, "_elbow.png"), 1800, 1600, ElbowPlot(seur, ndims = 50))

  seur <- FindNeighbors(seur, dims = 1:50)
  seur <- FindClusters(seur, resolution = 0.8, algorithm = 3, verbose = FALSE)
  seur <- RunUMAP(seur, reduction = "pca", dims = 1:50, assay = "RNA")
  saveRDS(seur, output_rds)

  # 3. Reload the clustered object and attach sample metadata
  seur <- readRDS(output_rds)
  seur <- add_sample_metadata(seur, metadata_csv, cluster_name)
  saveRDS(seur, output_rds)

  # 4. Visualize the clustered object and summarize composition
  save_plot(paste0(project, "_umap.png"), 1800, 1600, DimPlot(seur, reduction = "umap", order = TRUE) + ggtitle(""))
  save_plot(paste0(project, "_umap_labeled.png"), 1800, 1600, DimPlot(seur, reduction = "umap", label = TRUE, order = TRUE) + ggtitle(""))
  save_plot(paste0(project, "_umap_sample.png"), 1800, 1600, DimPlot(seur, reduction = "umap", group.by = "orig.ident", order = TRUE))
  save_plot(paste0(project, "_umap_sample_split.png"), 1800, 1600, DimPlot(seur, reduction = "umap", split.by = "orig.ident", ncol = 3, order = TRUE))
  save_plot(paste0(project, "_umap_Sex.png"), 1800, 1600, DimPlot(seur, reduction = "umap", group.by = "Sex", order = TRUE))
  save_plot(paste0(project, "_umap_Sex_split.png"), 1800, 1600, DimPlot(seur, reduction = "umap", split.by = "Sex", ncol = 3, order = TRUE))
  save_plot(paste0(project, "_umap_Sex_splitbyCelltype.png"), 1800, 1600, DimPlot(seur, reduction = "umap", split.by = "Sex", group.by = "predicted.celltype.l1", ncol = 3, order = TRUE))
  save_plot(paste0(project, "_umap_Celltype.png"), 1800, 1600, DimPlot(seur, reduction = "umap", group.by = "predicted.celltype.l1", order = TRUE))

  seurat_clusters_table <- data.frame(table(seur$seurat_clusters))
  save_table_png(seurat_clusters_table, paste0(project, "_seurat_sample_counts.png"), 800, 2000)

  cluster_sample_counts <- data.frame(table(seur$seurat_clusters, seur$Sex))
  cluster_sample_counts <- dcast(cluster_sample_counts, Var1 ~ Var2)
  rownames(cluster_sample_counts) <- cluster_sample_counts[, 1]
  cluster_sample_counts <- cluster_sample_counts[, -1]
  for (i in c(1, 2)) {
    cell_total <- sum(cluster_sample_counts[[i]])
    new_col <- paste0(colnames(cluster_sample_counts)[i], "_freq")
    cluster_sample_counts[[new_col]] <- round((cluster_sample_counts[[i]] / cell_total) * 100, 2)
  }
  cluster_sample_counts <- cbind(cluster = rownames(cluster_sample_counts), cluster_sample_counts)
  save_table_png(cluster_sample_counts, paste0(project, "_seurat_sex_counts.png"), 1200, 2000)

  cluster_sample_counts <- data.frame(table(seur$predicted.celltype.l1, seur$Sex))
  cluster_sample_counts <- dcast(cluster_sample_counts, Var1 ~ Var2)
  rownames(cluster_sample_counts) <- cluster_sample_counts[, 1]
  cluster_sample_counts <- cluster_sample_counts[, -1]
  for (i in c(1, 2)) {
    cell_total <- sum(cluster_sample_counts[[i]])
    new_col <- paste0(colnames(cluster_sample_counts)[i], "_freq")
    cluster_sample_counts[[new_col]] <- round((cluster_sample_counts[[i]] / cell_total) * 100, 2)
  }
  cluster_sample_counts <- cbind(cluster = rownames(cluster_sample_counts), cluster_sample_counts)
  save_table_png(cluster_sample_counts, paste0(project, "_celltype_sex_counts.png"), 1200, 800)

  # 5. Run module-score and marker-expression summaries
  gene_list <- list(c("IFI27", "IFI6", "IFI44", "IFI44L", "USP18", "LY6E", "OAS1",
                      "ISG15", "IFIT1", "OAS3", "HERC5", "MX1", "LAMP3", "EPSTI1",
                      "IFIT3", "OAS2", "RTP4", "PLSCR1", "SPATS2L", "RSAD2", "SIGLEC1"))

  score_name <- if (cluster_name == "PBMC") "IFN_scores1" else "IFN_scores11"
  seur <- AddModuleScore(object = seur, features = gene_list, name = if (cluster_name == "PBMC") 'IFN_scores' else 'IFN_scores1')

  save_plot(paste0(project, "_moduleScore.png"), 1800, 1600, FeaturePlot(seur, features = score_name))
  save_plot(paste0(project, "_moduleScore_split.png"), 1800, 1600, FeaturePlot(seur, features = score_name, split.by = "Sex"))
  save_plot(paste0(project, "_IFN_scores_violin.png"), 4000, 4000, VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = 0))
  save_plot(paste0(project, "_IFN_scores_split_violin.png"), 5000, 4000, VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = 0, split.by = "Sex"))

  seur$combined <- "Sample"
  seur <- SetIdent(seur, value = seur@meta.data$combined)
  save_plot(paste0(project, "_IFN_scores_split_violin_noCluster.png"), 2000, 4000, VlnPlot(seur, features = gene_list[[1]], stack = TRUE, flip = TRUE, pt.size = NULL, split.by = "Sex"))
  seur <- SetIdent(seur, value = seur@meta.data$seurat_clusters)

  save_plot(paste0(project, "_FCGR3BScore.png"), 1800, 1600, FeaturePlot(seur, features = "FCGR3B"))
  save_plot(paste0(project, "_FCGR3BScore_split.png"), 1800, 1600, FeaturePlot(seur, features = "FCGR3B", split.by = "Sex"))
  save_plot(paste0(project, "_FCGR3BScore_split_Vln.png"), 3000, 1000, VlnPlot(seur, features = "FCGR3B", split.by = "Sex", pt.size = 0))

  immature_neutrophils <- c('LCN2', 'LTF', 'MPO', 'CEACAM6', 'OLR1', 'PI3', 'AZU1', 'RETN', 'SLC2A5', 'SLPI', 'BPI', 'TCN1', 'DEFA1B', 'CAMP', 'MS4A3')
  save_plot(paste0(project, "_immature_neutrophils.png"), 4000, 4000, FeaturePlot(seur, features = immature_neutrophils))
  save_plot(paste0(project, "_immature_neutrophils_violin.png"), 5000, 4000, VlnPlot(seur, features = immature_neutrophils, stack = TRUE, flip = TRUE, pt.size = 0))
  save_plot(paste0(project, "_immature_neutrophils_split_violin.png"), 7000, 4000, VlnPlot(seur, features = immature_neutrophils, stack = TRUE, flip = TRUE, pt.size = 0, split.by = "Sex"))

  # 6. Save the final clustered object and run the differential-expression comparison
  saveRDS(seur, output_rds)
  seur.DEG <- subset(seur, downsample = 200)
  DEG.M.vs.F <- FindMarkers(seur.DEG, ident.1 = 'Male', ident.2 = 'Female', group.by = 'Sex', test.use = 'MAST')
  write.csv(DEG.M.vs.F, 'DEG.M.vs.F.csv')

  markers12 <- c('CAMP', 'FCGR3B', 'LTF', 'MMP9', 'ELANE', 'S100A8', 'S100A9', 'DEFA4', 'AZU1', 'MPO', 'BPI', 'LCN2')
  save_plot(paste0(project, "_Markers12_FeaturePlot.png"), 4000, 2400, FeaturePlot(seur, features = markers12, ncol = 3))
  save_plot(paste0(project, "_Markers12_VlnPlot.png"), 4000, 2500, VlnPlot(seur, features = markers12, stack = TRUE, flip = TRUE, pt.size = 0))
  save_plot(paste0(project, "_Markers12_SexSplit_VlnPlot.png"), 4000, 2500, VlnPlot(seur, features = markers12, stack = TRUE, flip = TRUE, pt.size = 0, split.by = 'Sex'))

  saveRDS(seur, output_rds)
  writeLines(capture.output(sessionInfo()), paste0(project, '_sessionInfo.txt'))
}

# 2. Rebuild the neutrophil and PBMC merged analyses with the same workflow
run_merged_dataset(
  input_rds = "/data/Amblerwg_Schaughency/projects/0624_SLE/Merge/merged_neutrophil.rds",
  workdir = "/data/Amblerwg_Schaughency/projects/0624_SLE/Merged_neutrophils_analysis/",
  cluster_name = "neutrophil",
  output_rds = "neutrophils_seur_clusters.rds"
)

run_merged_dataset(
  input_rds = "/data/Amblerwg_Schaughency/projects/0624_SLE/Merge/merged_PBMC.rds",
  workdir = "/data/Amblerwg_Schaughency/projects/0624_SLE/Merged_PBMC_analysis/",
  cluster_name = "PBMC",
  output_rds = "PBMC_seur_clusters.rds"
)
