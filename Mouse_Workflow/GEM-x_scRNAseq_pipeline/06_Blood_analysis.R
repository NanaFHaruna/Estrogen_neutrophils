# Streamlined blood/spleen analysis.
# Keeps the same outputs, but removes repeated subset / DEG / volcano / score blocks.

source('/mnt/data/gem_x_analysis_helpers.R')

setwd('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/')
seur <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/Blood_spleen_clean_gem_x2.rds')
DefaultAssay(seur) <- 'RNA'

# -----------------------------------------------------------------------------
# Initial cluster identification and annotation
# -----------------------------------------------------------------------------
DimPlot(seur, label = TRUE)
table(seur@meta.data$seurat_clusters)

seur_down <- subset(seur, downsample = 10000)
DimPlot(seur_down, label = TRUE)
seur_down <- ScaleData(seur_down, features = rownames(seur_down))

neutro.markers <- FindAllMarkers(seur_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(neutro.markers, 'original_markers_blood&spleen_clustersv2.csv')

top10 <- neutro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(subset(seur, downsample = 1000), features = top10$gene) + NoLegend()

seur <- rename_clusters(seur, cluster_map_8)
DimPlot(seur, label = TRUE)
DimPlot(seur, cols = sample_palette)
seur@active.ident <- factor(seur@active.ident, levels = c('G4', 'G5a', 'G5a/c', 'G5b', 'G5c'))
DotPlot(seur, features = c('Mmp8', 'Retnlg', 'S100a8', 'S100a9', 'Wfdc21', 'Egr1', 'Wfdc17', 'Ifit3', 'Rsad2', 'Isg15', 'Gngt2', 'Gm2a', 'Itga4')) + RotatedAxis()

# -----------------------------------------------------------------------------
# Split blood from spleen and define the main downstream objects
# -----------------------------------------------------------------------------
Idents(seur) <- 'Compartment'
blood <- subset(seur, idents = 'Blood')
Idents(blood) <- 'seurat_clusters'

SQ <- subset(blood, idents = 'SQ')
Pellet <- subset(blood, idents = 'Pellet')
ERKO <- subset(blood, idents = 'ERKO')
Idents(blood) <- 'orig.ident'
Chimera_long <- subset(blood, idents = c('CD451-WT-E2-Blood', 'CD452-ERKO-E2-Blood', 'CD451-WT-Veh-Blood', 'CD452-ERKO-Veh-Blood'))

Chimera_long$Condition <- dplyr::case_when(
  Chimera_long$orig.ident == 'CD451-WT-Veh-Blood' ~ 'WT-Veh',
  Chimera_long$orig.ident == 'CD452-ERKO-Veh-Blood' ~ 'ERKO-Veh',
  Chimera_long$orig.ident == 'CD451-WT-E2-Blood' ~ 'WT-E2',
  Chimera_long$orig.ident == 'CD452-ERKO-E2-Blood' ~ 'ERKO-E2',
  TRUE ~ NA_character_
)

# -----------------------------------------------------------------------------
# SQ
# -----------------------------------------------------------------------------
SQ_equal <- prepare_subset(SQ, 'Condition', 4450, cluster_map_8, c('Sham', 'Oop', 'E2'))
plot_branch_panels(SQ_equal, 'Condition', 'Rsad2')

IFNcluster <- subset(SQ_equal, idents = 'G5b')
make_ifn_heatmap(IFNcluster, ifn_i_genes[[1]], column_split = c(1, 2, 3), column_gap_mm = 1, show_colnames = TRUE, fontsize = 12)

sq_deg <- run_deg_block(
  SQ_equal, 'Condition', 'E2', 'Oop', 2000,
  'SQ_EvsVeh_DEG.csv', xlim = c(-1.5, 1.5), ylim = c(0, 180)
)
plot_deg_volcano(sq_deg$deg_tbl, c(-1.5, 1.5), c(0, 180), blood_select_main, lab_size = 5, add_connectors = TRUE)

SQ_equal <- add_ifn_scores(SQ_equal)
stats_SQ <- compare_scores(
  SQ_equal, 'Condition', c('IFNIscores501', 'HallmarkIFN1'),
  list(c('E2', 'Oop'), c('E2', 'Sham'), c('Sham', 'Oop'))
)
print(stats_SQ)

# -----------------------------------------------------------------------------
# Pellet
# -----------------------------------------------------------------------------
Pellet_equal <- prepare_subset(Pellet, 'Condition', 5475, cluster_map_8, c('Oop', 'E2'))
plot_branch_panels(Pellet_equal, 'Condition', 'Rsad2')

pellet_deg <- run_deg_block(
  Pellet_equal, 'Condition', 'E2', 'Oop', 2000,
  'Pellet_EvsVeh_DEG.csv', xlim = c(-1.5, 1.5), ylim = c(0, 155)
)
plot_deg_volcano(pellet_deg$deg_tbl, c(-1.5, 1.5), c(0, 155), blood_select_main, lab_size = 4, add_connectors = TRUE)

Pellet_equal <- add_ifn_scores(Pellet_equal)
stats_Pellet <- compare_scores(
  Pellet_equal, 'Condition', c('IFNIscores501', 'HallmarkIFN1'),
  list(c('E2', 'Oop'))
)
print(stats_Pellet)

# -----------------------------------------------------------------------------
# ERKO
# -----------------------------------------------------------------------------
ERKO_equal <- prepare_subset(ERKO, 'Condition', 6206, cluster_map_7, c('WT', 'KO'))
plot_branch_panels(ERKO_equal, 'Condition', 'Rsad2')

ERKO_deg <- run_deg_block(
  ERKO_equal, 'Condition', 'WT', 'KO', 2000,
  'ERKO_WTvsKO_DEG.csv', xlim = c(-1, 1.75), ylim = c(0, 210)
)
plot_deg_volcano(ERKO_deg$deg_tbl, c(-1, 1.75), c(0, 210), blood_select_main, lab_size = 4, add_connectors = TRUE)
plot_deg_volcano(ERKO_deg$deg_tbl, c(-0.5, 1.75), c(0, 210), blood_select_alt, lab_size = 4, add_connectors = TRUE)

ERKO_equal <- add_ifn_scores(ERKO_equal)
stats_ERKO <- compare_scores(
  ERKO_equal, 'Condition', c('IFNIscores501', 'HallmarkIFN1'),
  list(c('WT', 'KO'))
)
print(stats_ERKO)

# -----------------------------------------------------------------------------
# Chimera long
# -----------------------------------------------------------------------------
Chimera_long_equal <- prepare_subset(
  Chimera_long, 'orig.ident', 7427, cluster_map_7,
  c('CD451-WT-Veh-Blood', 'CD452-ERKO-Veh-Blood', 'CD451-WT-E2-Blood', 'CD452-ERKO-E2-Blood')
)
plot_branch_panels(Chimera_long_equal, 'orig.ident', 'Rsad2')

make_ifn_heatmap(Chimera_long_equal, ifn_i_genes[[1]], column_split = c(4, 2, 3, 1), column_gap_mm = 1, show_colnames = TRUE, fontsize = 12)

chimera_long_deg <- run_deg_block(
  Chimera_long_equal, 'orig.ident', 'CD451-WT-E2-Blood', 'CD452-ERKO-E2-Blood', 2000,
  'Chimera_long_E2_WTvsKO_DEG.csv', xlim = c(-0.75, 1.5), ylim = c(0, 100)
)
plot_deg_volcano(chimera_long_deg$deg_tbl, c(-0.75, 1.5), c(0, 100), blood_select_main, lab_size = 4, add_connectors = TRUE)
plot_deg_volcano(chimera_long_deg$deg_tbl, c(-0.5, 1.75), c(0, 100), blood_select_alt, lab_size = 4, add_connectors = TRUE)

Chimera_long_equal <- add_ifn_scores(Chimera_long_equal)
stats_Chimera_long <- compare_scores(
  Chimera_long_equal, 'orig.ident', c('IFNIscores501', 'HallmarkIFN1'),
  list(
    c('CD451-WT-Veh-Blood', 'CD452-ERKO-Veh-Blood'),
    c('CD451-WT-E2-Blood', 'CD452-ERKO-E2-Blood')
  )
)
print(stats_Chimera_long)


