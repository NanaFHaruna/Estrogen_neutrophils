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

# 1. Reproduce the downstream AIS plot and DEG summaries with shared helper functions.
#    The original script contained a number of legacy blocks that depended on undefined
#    objects (AIS2, AIS_orig, AIS3_downsampled_by_hormone). Those are removed here so
#    the file stays runnable and focused on the analyses supported by the loaded inputs.

source('AIS_shared_helpers.R')

project <- 'NIAMS-34'
workdir <- '/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge_Analysis/'
setwd(workdir)

# 1a. Load the two objects used for the downstream visualizations.
seur.PBMC <- readRDS('/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Merge_Analysis/PBMC_seur_clusters.rds')
AIS <- readRDS('/data/Amblerwg_Schaughency/will/AIS/202412_secondintegration/Integrate/integrated_SeuratInt_neu_clustered.rds')

# 1b. PBMC overview and neutrophil-focused filters.
save_png(paste0(project, '_PBMC_marker_featureplot.png'), 8, 6, {
  p <- FeaturePlot(seur.PBMC, features = c('CSF3R', 'FCGR3B', 'PPBP', 'IGKC'), combine = FALSE)
  print(wrap_plots(p, ncol = 2))
})

Idents(seur.PBMC) <- 'seurat_clusters'
save_png(paste0(project, '_PBMC_clusters.png'), 6, 5.3, print(DimPlot(seur.PBMC, label = TRUE)))

seur.PBMC1 <- subset(seur.PBMC, idents = c(4, 23, 9, 13, 15), invert = TRUE)
Idents(seur.PBMC1) <- 'predicted.celltype.l1'
save_png(paste0(project, '_PBMC_celltypes_filtered.png'), 6, 5.3, print(DimPlot(seur.PBMC1, label = FALSE)))
save_png(paste0(project, '_PBMC_IFN_scores.png'), 6, 5.3, print(FeaturePlot(seur.PBMC1, features = 'IFN_scores11')))

# 1c. Re-label the integrated neutrophil clusters and inspect the IFN signal.
AIS <- subset(AIS, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
AIS <- rename_clusters(AIS, c('Nh1/2', 'Nh2', 'Nh1', 'Nh3', 'Nh3', 'Nh0', 'Nh0', 'Nh4', 'Nh3'))
save_png(paste0(project, '_AIS_IFN_scores.png'), 6, 5.3, print(FeaturePlot(AIS, features = 'IFN_scores11')))

IFN_1_cluster <- subset(AIS, idents = 'Nh3')
save_png(paste0(project, '_AIS_Nh3_IFN_scores_split.png'), 8, 5.3, print(FeaturePlot(IFN_1_cluster, features = 'IFN_scores11', split.by = 'Identifier2', combine = FALSE)))
save_png(paste0(project, '_AIS_Nh3_IFN_violin.png'), 8, 5.3, print(vln_box(IFN_1_cluster, 'IFN_scores11', 'Identifier2', pt.size = 0, ylim = c(0, 0.75))))
save_png(paste0(project, '_AIS_Nh3_OAZ1_feature.png'), 8, 5.3, print(FeaturePlot(IFN_1_cluster, features = 'OAZ1', split.by = 'Identifier2')))
save_png(paste0(project, '_AIS_Nh3_OAZ1_violin.png'), 8, 5.3, print(vln_box(IFN_1_cluster, 'OAZ1', 'Identifier2', pt.size = 0, ylim = c(0, 1700))))
save_png(paste0(project, '_AIS_Nh3_density.png'), 8, 5.3, {
  p <- UMAPPlot(IFN_1_cluster, split.by = 'Identifier2', combine = FALSE, ncol = 2)
  print(p[[1]] + stat_density_2d(aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level)),
                                 geom = 'density_2d_filled', colour = 'ivory', alpha = 0.2,
                                 contour_var = 'ndensity'))
})

# 1d. Compare IFN-I/II and IFN-I-all states within the AIS neutrophil clusters.
IFN_1_cluster$Identifier <- factor(IFN_1_cluster$Identifier, levels = c('HC', 'AIS'))
IFN_1_cluster$Identifier2 <- factor(IFN_1_cluster$Identifier2, levels = c('HC', 'AIS', 'AIS_plus_hormone', 'AIS_no_hormone'))
save_png(paste0(project, '_AIS_Nh3_IFN_scores_by_ID.png'), 8, 5.3, print(vln_box(IFN_1_cluster, 'IFN_scores11', 'Identifier', pt.size = 0)))
save_png(paste0(project, '_AIS_Nh3_IFN_scores_by_ID2.png'), 8, 5.3, print(vln_box(IFN_1_cluster, 'IFN_scores11', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_Nh3_IFN_scores_by_patient.png'), 8, 5.3, print(vln_box(IFN_1_cluster, 'IFN_scores11', 'Patient.ID', pt.size = 0)))
save_png(paste0(project, '_AIS_Nh3_IFN_scores_downsampled.png'), 8, 5.3, print(vln_box(subset(IFN_1_cluster, downsample = 300), 'IFN_scores11', 'Identifier2', pt.size = 0)))

IFN_I_all_cluster <- subset(AIS, idents = c('Nh3-IFN-I', 'Nh3-IFN-II', 'imm-Nh3'))
save_png(paste0(project, '_AIS_IFN_Iall_scores.png'), 8, 5.3, print(vln_box(IFN_I_all_cluster, 'IFN_scores11', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_Iall_scores_downsampled.png'), 8, 5.3, print(vln_box(subset(IFN_I_all_cluster, downsample = 800), 'IFN_scores11', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_Iall_scores_by_patient.png'), 8, 5.3, print(vln_box(IFN_I_all_cluster, 'IFN_scores11', 'Patient.ID', pt.size = 0)))

IFN_II_cluster <- subset(AIS, idents = 'Nh3-IFN-II')
save_png(paste0(project, '_AIS_IFN_II_scores.png'), 8, 5.3, print(vln_box(IFN_II_cluster, 'IFN_scores11', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_II_scores_downsampled.png'), 8, 5.3, print(vln_box(subset(IFN_II_cluster, downsample = 300), 'IFN_scores11', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_II_scores_by_patient.png'), 8, 5.3, print(vln_box(IFN_II_cluster, 'IFN_scores11', 'Patient.ID', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_II_GBP5.png'), 8, 5.3, print(vln_box(IFN_II_cluster, 'GBP5', 'Identifier2', pt.size = 0)))
save_png(paste0(project, '_AIS_IFN_II_STAT1_downsampled.png'), 8, 5.3, print(vln_box(subset(IFN_II_cluster, downsample = 500), 'STAT1', 'Identifier2', pt.size = 0)))

# 1e. Differential expression using the same identifiers and MAST test as the original script.
DEG_hormone <- run_deg(AIS, id_col = 'Identifier2', downsample_n = 1000)
write.csv(DEG_hormone, '202412DEG.AISvsAIS_plus_hormone.csv')

DEG_nohormone <- run_deg(AIS, id_col = 'Identifier', downsample_n = 2000)
write.csv(DEG_nohormone, '2026DEG.HCvsAIS_allbutnotincluding_gonadectomy_withouthormone.csv')

data1 <- deg_to_df(DEG_hormone)
data3 <- read_deg_for_volcano('2026DEG.HCvsAIS_allbutnotincluding_gonadectomy_withouthormone.csv')

# 1f. Volcano plots, kept in the same style but wrapped in a reusable prep step.
vol_prep <- make_volcano_prep(DEG_hormone, cutoff = 5)
p1 <- EnhancedVolcano(
  data1,
  lab = data1$gene,
  x = 'avg_log2FC', y = 'p_val_adj',
  selectLab = c('IFI440', 'HERC5', 'EPSTI1', 'IFI6', 'ISG15', 'LYGE', 'PLSCR1', 'OAS2', 'OAS3', 'IFIT1', 'LAMP3', 'USP18', 'SIGLEC1', 'IFI27', 'TRP4', 'SPATS2L', 'MX1', 'RSAD2', 'IL18RAP', 'EIF2AK2', 'S100A9', 'DDX60L', 'TRIM22', 'FKBP5', 'CMPK2', 'PRH1'),
  pCutoff = 0.1, FCcutoff = 0.25,
  drawConnectors = TRUE, legendPosition = 'none'
)

p2 <- EnhancedVolcano(
  vol_prep$capped,
  lab = vol_prep$capped$gene,
  x = 'avg_log2FC', y = 'p_val_adj',
  selectLab = c('IFI440', 'HERC5', 'EPSTI1', 'IFI6', 'ISG15', 'LYGE', 'PLSCR1', 'OAS2', 'OAS3', 'IFIT1', 'LAMP3', 'USP18', 'SIGLEC1', 'IFI27', 'TRP4', 'SPATS2L', 'MX1', 'RSAD2', 'OAS3', 'IL18RAP', 'EIF2AK2', 'DDX60L', 'DDX60', 'TRIM22', 'FKBP5', 'CMPK2', 'PRH1', 'IFIT3', 'NLRC5', 'FGD4', 'LOXL1', 'PADI4', 'TEX2', 'CERS3', 'IGSF11', 'FOSB', 'S100A8', 'S100A9', 'FCGR3B'),
  pCutoff = 0.05, FCcutoff = 0.25,
  drawConnectors = TRUE,
  shapeCustom = vol_prep$shape,
  pointSize = vol_prep$size,
  legendPosition = 'none',
  xlim = c(-5, 5), labSize = 5, labFace = 'bold'
)

save_png(paste0(project, '_volcano_AIS_plus_hormone.png'), 8, 6, print(p1))
save_png(paste0(project, '_volcano_AIS_plus_hormone_custom.png'), 8, 6, print(p2))

save_png(paste0(project, '_volcano_AIS_no_hormone.png'), 8, 6, {
  print(EnhancedVolcano(
    data3,
    lab = data3$gene,
    x = 'avg_log2FC', y = 'p_val_adj',
    selectLab = c('CLEC12A', 'HLA-A', 'MNDA', 'TXNIP', 'NCF1', 'S100A9', 'NCOA4', 'B2M', 'CXCR2', 'IRF2', 'KDM6B', 'FOS', 'CEBPB', 'KDM6A', 'DDX60L', 'IRF1', 'HIF1A', 'NUP98', 'FKBP5', 'SAT1', 'MME', 'CXCL8', 'PDE4B', 'GBP5', 'SAMSN1', 'RABGEF1', 'BCL2A1', 'CCND3', 'STAT5B', 'TNFRSF1A', 'IKZF1', 'CASP4'),
    pCutoff = 0.1, FCcutoff = 0.25,
    xlim = c(-2, 2), drawConnectors = TRUE, max.overlaps = Inf
  ))
})

save_png(paste0(project, '_volcano_AIS_vs_plus_hormone.png'), 8, 6, {
  print(EnhancedVolcano(
    data1,
    lab = data1$gene,
    x = 'avg_log2FC', y = 'p_val_adj',
    selectLab = c('CLEC12A', 'TMTC1', 'FGF13', 'IFI44', 'TNFSF13B', 'DOCK4', 'HERC5', 'TRIM22', 'MX1', 'DDX60L', 'EIF2AK2', 'MX2', 'TET2', 'NAMPT', 'ALPK1', 'TLR2', 'B2M', 'CXCL8', 'FOS', 'NCF1', 'PLCG2'),
    pCutoff = 0.1, FCcutoff = 0.25,
    xlim = c(-2, 2), drawConnectors = TRUE, max.overlaps = Inf
  ))
})

# 1g. Marker gene and summary plots.
neutro.markers <- FindAllMarkers(AIS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- neutro.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
save_png(paste0(project, '_AIS_top10_heatmap.png'), 12, 8, print(DoHeatmap(subset(AIS, downsample = 1000), features = top10$gene) + NoLegend()))

# 1h. Lightweight summaries retained from the original script.
mylist <- sapply(levels(IFN_I_all_cluster$Patient.ID), function(x) mean(IFN_I_all_cluster$IFN_scores11[IFN_I_all_cluster$Patient.ID == x]))
mean_IFN_Iall_AIS <- mean(IFN_I_all_cluster$IFN_scores11[IFN_I_all_cluster$Identifier == 'AIS'])
mean_IFN_Iall_HC <- mean(IFN_I_all_cluster$IFN_scores11[IFN_I_all_cluster$Identifier == 'HC'])
wilcox_IFN_AIS_vs_HC <- wilcox.test(AIS@meta.data[AIS$Identifier == 'AIS', 'IFN_scores11'], AIS@meta.data[AIS$Identifier == 'HC', 'IFN_scores11'])
ifn_ratio <- prop.table(table(Idents(AIS), AIS$Identifier), margin = 2)

write_session_info(project)
