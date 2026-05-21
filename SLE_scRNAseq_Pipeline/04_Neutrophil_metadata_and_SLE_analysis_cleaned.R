# 1. Load libraries and define a compact sample-annotation lookup
library(Seurat)
library(dplyr)
library(magrittr)
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
library(EnhancedVolcano)

rm(list = ls())

project <- "NIAMS-34"
workdir <- "/data/Amblerwg_Schaughency/projects/0624_SLE/Integrate/"
setwd(workdir)

save_plot <- function(filename, width, height, plot_expr) {
  png(filename, width = width, height = height, res = 300)
  print(plot_expr)
  dev.off()
}

sample_annotations <- data.frame(
  orig.ident = c(
    'WA_001-CMO301','WA_001-CMO302','WA_002-CMO303','WA_002-CMO304','WA_003-CMO305','WA_003-CMO306',
    'WA_027-CMO301','WA_027-CMO302','WA_028-CMO303','WA_028-CMO304','WA029','WA030','WA031','WA032',
    'WA048','WA049','WA050','WA051','WA052','WA053','WA054','WA055','WA064','WA065','WA066','WA067',
    'WA070','WA072','WA078','WA093','WA079','WA094'
  ),
  Categories = c(
    'SLE M LDG','SLE M NDN','SLE F LDG','SLE F NDN','SLE M LDG','SLE M NDN',
    'SLE M LDG','SLE M NDN','SLE F LDG','SLE F NDN','SLE F NDN','SLE M NDN','SLE F LDG','SLE M LDG',
    'SLE F LDG','SLE F NDN','SLE M LDG','SLE M NDN','SLE M NDN','SLE M LDG','SLE F NDN','SLE F LDG','HC M NDN','HC F NDN','HC F LDG','HC F LDG',
    'HC F NDN','HC F NDN','SLE F NDN','SLE M NDN','SLE F LDG','SLE M LDG'
  ),
  Density = c(
    'LDG','NDN','LDG','NDN','LDG','NDN',
    'LDG','NDN','LDG','NDN','NDN','NDN','LDG','LDG',
    'LDG','NDN','LDG','NDN','NDN','LDG','NDN','LDG','NDN','NDN','LDG','LDG',
    'NDN','NDN','NDN','NDN','LDG','LDG'
  ),
  Condition = c('SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'SLE', 'HC', 'HC', 'HC', 'HC', 'HC', 'HC', 'SLE', 'SLE', 'SLE', 'SLE'),
  Age_category = c(
    '<65','<65','<65','<65','<65','<65',
    '<65','<65','<65','<65','>65','>65','>65','>65',
    '<65','<65','<65','<65','<65','<65','<65','<65','<65','<65','<65','<65',
    '<65','<65','>65','>65','>65','>65'
  ),
  stringsAsFactors = FALSE
)

# 2. Load the integrated neutrophil object and attach the metadata fields used throughout the analysis
Neu <- readRDS('integrated_SeuratInt_neutrophil_clustered_total.rds')
DimPlot(Neu, label = TRUE)
DefaultAssay(Neu) <- 'integrated'
Neu <- FindClusters(Neu, resolution = c(0.2, 0.3, 0.4, 0.5))

annotation_idx <- match(Neu$orig.ident, sample_annotations$orig.ident)
Neu$Categories <- sample_annotations$Categories[annotation_idx]
Neu$Density <- sample_annotations$Density[annotation_idx]
Neu$Condition <- sample_annotations$Condition[annotation_idx]
Neu$Age_category <- sample_annotations$Age_category[annotation_idx]

Idents(Neu) <- 'integrated_snn_res.0.4'

table(Neu$Categories)
table(Idents(Neu), Neu$Categories)
prop.table(table(Idents(Neu), Neu$Categories), margin = 2)

# 3. Drop the low-information cluster and add the IFN module score
Neu2 <- subset(Neu, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
DimPlot(Neu2, label = TRUE)
FeaturePlot(Neu, features = 'S100A9', split.by = 'Density')
DefaultAssay(Neu2) <- 'RNA'

gene_list <- list(c('IFI27', 'IFI6', 'IFI44', 'IFI44L', 'USP18', 'LY6E', 'OAS1', 'ISG15', 'IFIT1', 'OAS3', 'HERC5', 'MX1', 'LAMP3', 'EPSTI1', 'IFIT3', 'OAS2', 'RTP4', 'PLSCR1', 'SPATS2L', 'RSAD2', 'SIGLEC1'))
Neu2 <- AddModuleScore(object = Neu2, features = gene_list, name = 'IFN_scores_all')
FeaturePlot(Neu2, features = 'IFN_scores_all1', split.by = 'Categories')
FeaturePlot(Neu2, features = 'CAMP', split.by = 'Condition')
DimPlot(Neu2, label = TRUE)
saveRDS(Neu2, '0724_total_neutrophils_subset.rds')

# 4. Reload the saved subset and recreate the marker/heatmap summaries from the original script
allNeu <- readRDS('0724_total_neutrophils_subset.rds')
DefaultAssay(allNeu) <- 'RNA'
allNeu <- ScaleData(allNeu, features = rownames(allNeu))

neutro.markers <- FindAllMarkers(allNeu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- neutro.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(subset(allNeu, downsample = 1000), features = top10$gene) + NoLegend()

# 5. Build the SLE subsets used in the downstream comparisons
SLEneu <- subset(x = allNeu, cells = WhichCells(allNeu, expression = Condition == 'SLE'))
NDN <- subset(allNeu, cells = WhichCells(allNeu, expression = Density == 'NDN'))
youngSLE <- subset(x = SLEneu, cells = WhichCells(SLEneu, expression = Age_category == '<65'))
olderSLE <- subset(x = SLEneu, cells = WhichCells(SLEneu, expression = Age_category == '>65'))
LDG <- subset(x = youngSLE, cells = WhichCells(youngSLE, expression = Density == 'LDG'))

DimPlot(allNeu, split.by = 'Condition')
prop.table(table(Idents(allNeu), allNeu$Categories), margin = 2)
prop.table(table(Idents(youngSLE), youngSLE$Categories), margin = 2)
prop.table(table(Idents(olderSLE), olderSLE$Categories), margin = 2)

# 6. Run the IFN score and sex comparisons from the original analysis
youngSLE <- AddModuleScore(object = youngSLE, features = gene_list, name = 'IFN_scores_all')
VlnPlot(youngSLE, features = 'IFN_scores_all1', group.by = 'Density', pt.size = 0) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)

Idents(youngSLE) <- 'Sex'
DEG <- subset(youngSLE, downsample = 1000)
SLEFvsMmarkers <- FindMarkers(DEG, ident.1 = 'Female', ident.2 = 'Male', test.use = 'MAST')
write.csv(SLEFvsMmarkers, '0724_youngFvsMmarkers.csv')

VlnPlot(NDN, features = 'IFN_scores11', group.by = 'Categories', pt.size = 0, raster = FALSE)
mean(NDN$IFN_scores11[NDN$Categories == 'SLE F NDN'])
mean(NDN$IFN_scores11[NDN$Categories == 'SLE M NDN'])
wilcox.test(NDN@meta.data[NDN$Categories == 'SLE F NDN', 'IFN_scores11'], NDN@meta.data[NDN$Categories == 'SLE M NDN', 'IFN_scores11'])
wilcox.test(NDN@meta.data[NDN$Categories == 'HC F NDN', 'IFN_scores11'], NDN@meta.data[NDN$Categories == 'HC M NDN', 'IFN_scores11'])
wilcox.test(LDG@meta.data[LDG$Sex == 'Female', 'IFN_scores11'], LDG@meta.data[LDG$Sex == 'Male', 'IFN_scores11'])
wilcox.test(NDN@meta.data[NDN$Condition == 'SLE', 'IFN_scores11'], NDN@meta.data[NDN$Condition == 'HC', 'IFN_scores11'])
wilcox.test(LDG@meta.data[LDG$Categories == 'SLE F LDG', 'IFN_scores_all1'], LDG@meta.data[LDG$Categories == 'SLE M LDG', 'IFN_scores_all1'])

Idents(youngSLE) <- 'seurat_clusters'
DimPlot(LDG, label = TRUE)

# 7. Save provenance information
writeLines(capture.output(sessionInfo()), 'SLE_analysis_sessionInfo.txt')
