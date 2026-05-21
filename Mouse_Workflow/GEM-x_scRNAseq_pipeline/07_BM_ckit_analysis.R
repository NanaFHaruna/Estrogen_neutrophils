# Streamlined BM/cKit analysis.
# Same branch logic and CSV outputs, but repeated DEG / volcano / module-score code is collapsed into helpers.

source('/mnt/data/gem_x_analysis_helpers.R')

setwd('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/Chimera/')
seur <- readRDS('/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/BM_ckit_gem_x2.rds')
DefaultAssay(seur) <- 'RNA'

# -----------------------------------------------------------------------------
# Initial annotation and cluster renaming
# -----------------------------------------------------------------------------
Idents(seur) <- 'Chemistry'
seur_down <- subset(seur, downsample = 30000)
Idents(seur) <- 'seurat_clusters'
Idents(seur_down) <- 'seurat_clusters'

DimPlot(seur, label = TRUE)
DimPlot(seur_down, label = TRUE)
write.csv(FindMarkers(seur_down, ident.1 = 39), 'cluster39_markers.csv')
write.csv(FindMarkers(seur_down, ident.1 = 42), 'cluster42_markers.csv')

# Marker panels used for cluster annotation in the original script.
marker_panels <- list(
  c('Procr', 'Hlf'),
  c('Dntt', 'Flt3', 'Bst2', 'Ccr9'),
  c('Cd79a', 'Cd74', 'Igkc', 'Cd19', 'Il7r', 'H2-Aa'),
  c('Hexb', 'Prtn3', 'Elane', 'Mpo'),
  c('Cebpe', 'Gfi1', 'Per3'),
  c('Irf8', 'Irf5', 'Ly86', 'Ccr2'),
  c('Mki67', 'Pcna', 'Top2a', 'Mcm6'),
  c('S100a8', 'Ly6g', 'Mmp9', 'Ltf'),
  c('Hba-a1', 'Pf4', 'Gata2', 'Prss34')
)
walk(marker_panels, ~ FeaturePlot(seur, features = .x))

bm_cluster_map <- c(
  'G5a', 'G5a', 'G5a/b', 'G4', 'G5c', 'G3', 'G4', 'G5c', 'G5c', 'G5b', 'G3', 'G2', 'Ly6C-GMP', 'G4', 'G3', 'proNeus', 'G4', 'Ly6C-GMP', 'pro-B cell',
  'preNeus', 'HSPC', 'proNeus', 'G4/G5b', 'CLP', 'pre/immature B cell', 'proNeus', 'mature B cell', 'pro-B cell', 'CDP', 'preNeus',
  'Erythroid lin', 'Erythroid lin', 'Monocyte/DC', 'Mk', 'Baso/Mast', 'G5b', 'G4/G5b', 'T cell lin', 'preNeus', 'proNeus', 'pro-B cell', 'pre/immature B cell',
  'HSPC', 'G4/G5b', 'G4'
)
seur <- rename_clusters(seur, bm_cluster_map)
bm_select_sq <- c('Thbs1', 'Myl10', 'Ngp', 'Chil3', 'S100a8', 'Hspa8', 'Vcam1', 'Tex2', 'Hk1', 'H2-Q7', 'Tap2', 'Psmb8', 'Psmb9', 'Ifitm2', 'Padi4')
bm_select_erko <- c('Thbs1', 'Myl10', 'Fgd4', 'Tex2', 'S100a8', 'Hspa8', 'Vcam1', 'Tex2', 'Hk1', 'H2-Q7', 'Tap2', 'Psmb8', 'Psmb9', 'Ifitm2', 'Padi4', 'Irf1', 'Ern1', 'Ifi27', 'H2-K1', 'B2m', 'Angpt1', 'Cpa3', 'Cxcl10', 'Stat1', 'Myc')
LabelClusters(DimPlot(seur) + NoLegend(), id = 'ident', color = 'black', size = 5, repel = TRUE, fontface = 'bold', box.padding = 1)

# -----------------------------------------------------------------------------
# Split by experiment, then analyze the relevant downstream branches
# -----------------------------------------------------------------------------
Idents(seur) <- 'Experiment'
SQ <- subset(seur, idents = 'SQ')
ERKO <- subset(seur, idents = 'ERKO')
Idents(seur) <- 'orig.ident'
Chimera_long <- subset(
  seur,
  idents = c(
    'CD451-WT-E2-BM', 'CD452-ERKO-E2-BM', 'CD451-WT-Veh-BM', 'CD452-ERKO-Veh-BM',
    'CD451-WT-E2-ckit', 'CD452-ERKO-E2-ckit', 'CD451-WT-Veh-ckit', 'CD452-ERKO-Veh-ckit'
  )
)

Chimera_long$Condition[Chimera_long$orig.ident %in% c('CD451-WT-E2-BM', 'CD451-WT-E2-ckit')] <- 'WT-E2'
Chimera_long$Condition[Chimera_long$orig.ident %in% c('CD451-WT-Veh-BM', 'CD451-WT-Veh-ckit')] <- 'WT-Veh'
Chimera_long$Condition[Chimera_long$orig.ident %in% c('CD452-ERKO-E2-BM', 'CD452-ERKO-E2-ckit')] <- 'ERKO-E2'
Chimera_long$Condition[Chimera_long$orig.ident %in% c('CD452-ERKO-Veh-BM', 'CD452-ERKO-Veh-ckit')] <- 'ERKO-Veh'

# -----------------------------------------------------------------------------
# Small reusable runners for the repeated branch analyses
# -----------------------------------------------------------------------------
run_single_deg_branch <- function(obj, group_col, ident1, ident2, csv_file, xlim, ylim,
                                 down_vln = NULL, vln_features = NULL, select_labs = NULL, lab_size = 6,
                                 vln_ylim_1 = NULL, vln_ylim_2 = NULL) {
  Idents(obj) <- group_col
  if (!is.null(down_vln) && !is.null(vln_features)) {
    obj_vln <- subset(obj, downsample = down_vln)
    if (length(vln_features) >= 1) {
      p <- VlnPlot(obj_vln, features = vln_features[1], pt.size = 0) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
      if (!is.null(vln_ylim_1)) p <- p + ylim(vln_ylim_1[1], vln_ylim_1[2])
      print(p)
    }
    if (length(vln_features) >= 2) {
      p <- VlnPlot(obj_vln, features = vln_features[2], pt.size = 0) +
        geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
      if (!is.null(vln_ylim_2)) p <- p + ylim(vln_ylim_2[1], vln_ylim_2[2])
      print(p)
    }
  }

  deg <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, test.use = 'MAST')
  write.csv(deg, csv_file)
  deg_tbl <- format_deg(deg)
  plot_deg_volcano(deg_tbl, xlim, ylim)
  if (!is.null(select_labs)) {
    plot_deg_volcano(deg_tbl, xlim, ylim, select_labs, lab_size = lab_size, add_connectors = TRUE)
  }
  list(deg = deg, deg_tbl = deg_tbl)
}

run_pair_deg_branch <- function(obj, group_col, ident1_e2, ident2_e2, csv_e2, ident1_veh, ident2_veh, csv_veh,
                                xlim_e2, ylim_e2, xlim_veh, ylim_veh, select_labs = NULL, lab_size = 6,
                                down_vln = NULL, vln_features = NULL, extra_gene = NULL, extra_gene_group = NULL) {
  Idents(obj) <- group_col
  obj_down <- subset(obj, downsample = 1000)

  if (!is.null(down_vln) && !is.null(vln_features)) {
    obj_vln <- subset(obj, downsample = down_vln)
    if (length(vln_features) >= 1) print(VlnPlot(obj_vln, features = vln_features[1], pt.size = 0) + geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1))
    if (length(vln_features) >= 2) print(VlnPlot(obj_vln, features = vln_features[2], pt.size = 0) + geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1))
  }
  if (!is.null(extra_gene)) {
    gene_de <- FindMarkers(obj_down, ident.1 = ident1_e2, ident.2 = ident2_e2, group.by = extra_gene_group, features = extra_gene, test.use = 'wilcox')
    invisible(gene_de)
  }

  deg_e2 <- FindMarkers(obj, ident.1 = ident1_e2, ident.2 = ident2_e2, test.use = 'MAST')
  write.csv(deg_e2, csv_e2)
  deg_e2_tbl <- format_deg(deg_e2)
  plot_deg_volcano(deg_e2_tbl, xlim_e2, ylim_e2)
  if (!is.null(select_labs)) plot_deg_volcano(deg_e2_tbl, xlim_e2, ylim_e2, select_labs, lab_size = lab_size, add_connectors = TRUE)

  deg_veh <- FindMarkers(obj_down, ident.1 = ident1_veh, ident.2 = ident2_veh, test.use = 'MAST')
  write.csv(deg_veh, csv_veh)
  deg_veh_tbl <- format_deg(deg_veh)
  plot_deg_volcano(deg_veh_tbl, xlim_veh, ylim_veh)
  if (!is.null(select_labs)) plot_deg_volcano(deg_veh_tbl, xlim_veh, ylim_veh, select_labs, lab_size = lab_size, add_connectors = TRUE)

  list(deg_e2 = deg_e2, deg_e2_tbl = deg_e2_tbl, deg_veh = deg_veh, deg_veh_tbl = deg_veh_tbl)
}

# -----------------------------------------------------------------------------
# SQ: proNeus
# -----------------------------------------------------------------------------
SQ_equal <- prepare_subset(SQ, 'Condition', 11626, bm_cluster_map, c('Sham', 'Oop', 'E2'))
plot_branch_panels(SQ_equal, 'Condition', 'Rsad2')

SQ_proNeu <- subset(SQ_equal, idents = 'proNeus')
Idents(SQ_proNeu) <- 'Condition'
SQ_proNeuvln <- subset(SQ_proNeu, downsample = 882)
VlnPlot(SQ_proNeuvln, features = 'Stat1', pt.size = 0) + geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
VlnPlot(SQ_proNeuvln, features = 'Irf9', pt.size = 0)

sq_proneu_deg <- run_single_deg_branch(
  SQ_proNeu, 'Condition', 'E2', 'Oop', 'SQ_proNeu_EvsVeh_DEG.csv',
  xlim = c(-0.75, 2), ylim = c(0, 240),
  select_labs = bm_select_sq,
  lab_size = 6
)
SQ_equal <- add_ifn_scores(SQ_equal)
print(compare_scores(SQ_equal, 'Condition', c('IFNIscores501', 'HallmarkIFN1'), list(c('KO', 'WT'))))

# -----------------------------------------------------------------------------
# ERKO: proNeus
# -----------------------------------------------------------------------------
ERKO_equal <- prepare_subset(ERKO, 'Condition', 18789, bm_cluster_map, c('WT', 'KO'))
plot_branch_panels(ERKO_equal, 'Condition', 'Rsad2')

ERKO_proNeu <- subset(ERKO_equal, idents = 'proNeus')
Idents(ERKO_proNeu) <- 'Condition'
ERKO_proNeuvln <- subset(ERKO_proNeu, downsample = 1401)
VlnPlot(ERKO_proNeuvln, features = 'Stat1', pt.size = 0) + geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
VlnPlot(ERKO_proNeuvln, features = 'Irf9', pt.size = 0)

ERKO_proNeu_deg <- run_single_deg_branch(
  ERKO_proNeu, 'Condition', 'WT', 'KO', 'ERKO_proNeu_EvsVeh_DEG.csv',
  xlim = c(-1.2, 1.2), ylim = c(0, 290),
  select_labs = bm_select_erko,
  lab_size = 6
)
ERKO_equal <- add_ifn_scores(ERKO_equal)
print(compare_scores(ERKO_equal, 'Condition', c('IFNIscores501', 'HallmarkIFN1'), list(c('KO', 'WT'))))

# -----------------------------------------------------------------------------
# Chimera long: shared preprocessing and then branch-by-branch DEG comparisons
# -----------------------------------------------------------------------------
Chimera_long_equal <- prepare_subset(
  Chimera_long, 'Condition', 17554, bm_cluster_map,
  c('WT-Veh', 'ERKO-Veh', 'WT-E2', 'ERKO-E2')
)
plot_branch_panels(Chimera_long_equal, 'Condition', 'Rsad2')

Chimera_long_proNeu <- subset(Chimera_long_equal, idents = 'proNeus')
Idents(Chimera_long_proNeu) <- 'Condition'
Chimera_long_proNeuvln <- subset(Chimera_long_proNeu, downsample = 1384)
VlnPlot(Chimera_long_proNeuvln, features = 'Padi4', pt.size = 0) + geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 1)
VlnPlot(Chimera_long_proNeuvln, features = 'Stat1', pt.size = 0)

invisible(FindMarkers(Chimera_long_proNeuvln, ident.1 = 'WT-E2', ident.2 = 'ERKO-E2', group.by = 'Condition', features = 'Tap2', test.use = 'wilcox'))

chimera_proneu_deg <- run_pair_deg_branch(
  Chimera_long_proNeu, 'Condition',
  'WT-E2', 'ERKO-E2', 'Chimera_long_proNeu_E_WTvsKO_DEG.csv',
  'WT-Veh', 'ERKO-Veh', 'Chimera_long_proNeu_Veh_WTvsKO_DEG.csv',
  xlim_e2 = c(-1, 2), ylim_e2 = c(0, 320),
  xlim_veh = c(-1, 2), ylim_veh = c(0, 330),
  select_labs = bm_select_main,
  lab_size = 6,
  down_vln = 1384, vln_features = c('Padi4', 'Stat1')
)
Chimera_long_proNeu <- add_ifn_scores(Chimera_long_proNeu)
print(compare_scores(Chimera_long_proNeu, 'Condition', c('IFNIscores501', 'HallmarkIFN1'), list(c('WT-Veh', 'ERKO-Veh'), c('WT-E2', 'ERKO-E2'))))

chimera_branch_cfgs <- list(
  list(subset_ident = 'Ly6C-GMP', csv_prefix = 'Chimera_long_Ly6C_GMP', x_e2 = c(-1, 2), y_e2 = c(0, 320), x_veh = c(-1, 2), y_veh = c(0, 330)),
  list(subset_ident = 'preNeus',   csv_prefix = 'Chimera_long_preNeus',   x_e2 = c(-1, 2), y_e2 = c(0, 320), x_veh = c(-1, 2), y_veh = c(0, 330)),
  list(subset_ident = 'HSPC',      csv_prefix = 'Chimera_long_HSPC',      x_e2 = c(-1, 2), y_e2 = c(0, 320), x_veh = c(-1, 2), y_veh = c(0, 330)),
  list(subset_ident = 'CDP',       csv_prefix = 'Chimera_long_CDP',       x_e2 = c(-1, 2), y_e2 = c(0, 320), x_veh = c(-1, 2), y_veh = c(0, 330)),
  list(subset_ident = 'CLP',       csv_prefix = 'Chimera_long_CLP',       x_e2 = c(-1, 2), y_e2 = c(0, 320), x_veh = c(-1, 2), y_veh = c(0, 330))
)

for (cfg in chimera_branch_cfgs) {
  sub <- subset(Chimera_long_equal, idents = cfg$subset_ident)
  Idents(sub) <- 'Condition'
  run_pair_deg_branch(
    sub, 'Condition',
    'WT-E2', 'ERKO-E2', paste0(cfg$csv_prefix, '_E_WTvsKO_DEG.csv'),
    'WT-Veh', 'ERKO-Veh', paste0(cfg$csv_prefix, '_Veh_WTvsKO_DEG.csv'),
    xlim_e2 = cfg$x_e2, ylim_e2 = cfg$y_e2,
    xlim_veh = cfg$x_veh, ylim_veh = cfg$y_veh,
    select_labs = bm_select_main,
    lab_size = 6
  )
}
