# Shared utilities for the streamlined gem_x analyses.
# Keep the logic the same; just remove repeated setup and plotting code.

suppressPackageStartupMessages({
  library(Seurat, lib.loc = "/data/Amblerwg_Schaughency/R4.3.0_lib")
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(EnhancedVolcano)
  library(ggrepel)
  library(scCustomize, lib.loc = "/data/Amblerwg_Schaughency/R4.3.0_lib")
  library(ComplexHeatmap)
  library(circlize)
})

sample_palette <- c("#A3A500", "#F8766D", "#00B0F6", "#E76BF3", "#00BF7D")

cluster_map_8 <- c("G5a", "G5c", "G5a", "G5b", "G5a", "G4", "G5a/c", "G5a")
cluster_map_7 <- c("G5a", "G5c", "G5a", "G5b", "G5a", "G4", "G5a/c")

ifn_hallmark_genes <- list(c(
  "Trim25","Il15","Stat2","Procr","Wars1","Psma3","Il4ra","Ly6e","Trim21","Elf1","Irf9","Psme2","Psme1","Ifit3","Casp8","Ifi27","Trim26","Csf1","Samd9l","Usp18","Psmb9",
  "Psmb8","Uba7","Cxcl10","Eif2ak2","C1s1","Isg15","Irf7","Cxcl11","Nub1","Txnip","Adar","Ripk2","Ifitm3","Gmpr","Lap3","Herc6","Lpar6","Ube2l6","Rtp4","Epsti1","Tmem140","Ifitm1","Bst2","Ifi35",
  "Ifih1","Pnpt1","Ogfr","Parp14","Ccrl2","Mvb12a","Cmtr1","Batf2","Trim14","Sp110","Trafd1","Gbp3","Nmi","Isg20","Rsad2","Dhx58","Parp9","Ifitm2","Rnf31","Ifi30","Tdrd7","Parp12","Slc25a28",
  "Oasl1","Oas1a","Ddx60","Helz2","Lamp3","Ifi44","Ncoa7","Tent5a","Trim12c","B2m","Cnp","Plscr1","Ifi44l","Cd74","Casp1","Il7","Irf1","Irf2","Cd47","Mov10","Mx2","Sell","Tap1","Ifit2","Lgals3bp","Cmpk2"
))

ifn_i_genes <- list(c(
  "Cxcl10", "Ifit3", "Rsad2", "Ifit1", "Gbp2", "Ifit3b", "Ifit2", "Ifit1bl1", "Herc6", "Usp18", "Ifi204", "Ifi47", "Isg20", "Pnp", "Cmpk2", "Znfx1", "Irf7", "Slfn5", "Ifi211",
  "Plac8", "Parp14", "Nt5c3", "Sp100", "Ctss", "Ifit3b", "Gbp5", "Gbp3", "Stat2", "Ifi209", "Fcgr1", "Slfn4", "Tor3a", "Ddx60", "H2-T23", "Socs1", "Eif2ak2", "Trim30c",
  "Oas3", "Oasl1", "Psmb9", "Psmb8", "Acer3", "Xaf1", "Sp110", "Trim30a"
))

blood_select_main <- c(
  "Ifit3", "Rsad2", "Ifit1", "Ifit3b", "Irf7", "Isg15", "Cd101", "Thbs1", "Padi4", "Acod1", "Nfkbia", "Nfkbiz", "Il1r2", "Wfdc17", "Cd14", "H2-K1", "H2-K7", "Cd244a",
  "Nr4a1", "Trem1", "Hsp90aa1", "Stat1", "Lrg1", "Lcn2", "Slnf4", "Ifi27l2a", "Oasl2", "S100a9", "Slfn1", "Slfn4", "Mmp8", "Ly6e", "Antxr2", "Fosb", "Ptsg2", "Il1b"
)

blood_select_alt <- c(
  "Ifit3", "Rsad2", "Ifit1", "Ifit3b", "Irf7", "Isg15", "Cd101", "Thbs1", "Padi4", "Nfkbiz", "Wfdc17", "H2-K1", "H2-K7", "Stat1", "Lrg1", "Lcn2", "Slnf4",
  "Ifi27l2a", "Oasl2", "S100a9", "Slfn1", "Slfn4", "Mmp8", "Ly6e", "Fosb", "Clec12a"
)

bm_select_main <- c(
  "Thbs1", "Myl10", "Fgd4", "Tex2", "S100a8", "Hspa8", "Vcam1", "Tex2", "Hk1", "H2-Q7", "Tap2", "Psmb8", "Psmb9", "Ifitm2", "Padi4", "Irf1", "Ern1", "Ifi27", "H2-K1",
  "B2m", "Angpt1", "Cpa3", "Cxcl10", "Stat1", "Myc", "Pgr"
)

bm_select_alt <- c(
  "Thbs1", "Myl10", "Fgd4", "Tex2", "S100a8", "Hspa8", "Vcam1", "Tex2", "Hk1", "H2-Q7", "Tap2", "Psmb8", "Psmb9", "Ifitm2", "Padi4"
)

# ---- Small helpers ---------------------------------------------------------

rename_clusters <- function(obj, cluster_map) {
  old_levels <- levels(obj)
  keep <- cluster_map[seq_along(old_levels)]
  names(keep) <- old_levels
  RenameIdents(obj, keep)
}

prepare_subset <- function(obj, group_col, downsample_n, cluster_map = NULL, group_levels = NULL) {
  Idents(obj) <- group_col
  out <- subset(obj, downsample = downsample_n)
  if (!is.null(cluster_map)) {
    Idents(out) <- "seurat_clusters"
    out <- rename_clusters(out, cluster_map)
  }
  if (!is.null(group_levels)) {
    out[[group_col]] <- factor(out[[group_col]], levels = group_levels)
  }
  out
}

plot_branch_panels <- function(obj, split_col, gene = "Rsad2", ncol = 2, use_palette = TRUE) {
  print(DimPlot(obj, split.by = split_col, ncol = ncol))
  if (use_palette) {
    print(DimPlot(obj, cols = sample_palette, split.by = split_col, ncol = ncol))
  }
  print(FeaturePlot(obj, features = gene, split.by = split_col, order = TRUE))
  print(plot_density(obj, gene) + facet_grid(as.formula(paste(".", split_col, sep = " ~ "))))
  umap <- UMAPPlot(obj, cols = sample_palette, split.by = split_col, combine = FALSE, ncol = ncol)
  print(
    umap[[1]] +
      stat_density_2d(
        aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"),
        geom = "density_2d_filled", colour = "ivory", alpha = 0.2, contour_var = "ndensity"
      )
  )
}

format_deg <- function(deg) {
  deg_tbl <- as.data.frame(deg, stringsAsFactors = FALSE)
  deg_tbl <- tibble::rownames_to_column(deg_tbl, var = "names")
  deg_tbl
}

plot_deg_volcano <- function(deg_tbl, xlim, ylim, select_labs = NULL, lab_size = 4, add_connectors = FALSE) {
  if (is.null(select_labs)) {
    EnhancedVolcano(deg_tbl, lab = deg_tbl$names, x = "avg_log2FC", y = "p_val_adj", pCutoff = 0.1,
                    FCcutoff = 0.25, xlim = xlim, ylim = ylim)
  } else {
    EnhancedVolcano(deg_tbl, lab = deg_tbl$names, selectLab = select_labs, x = "avg_log2FC", y = "p_val_adj",
                    pCutoff = 0.1, FCcutoff = 0.25, xlim = xlim, ylim = ylim,
                    max.overlaps = Inf, maxoverlapsConnectors = Inf,
                    drawConnectors = add_connectors, labSize = lab_size)
  }
}

run_deg_block <- function(obj, group_col, ident1, ident2, downsample_n, csv_file, xlim, ylim,
                          select_labs = NULL, lab_size = 4, add_connectors = FALSE) {
  Idents(obj) <- group_col
  obj_down <- subset(obj, downsample = downsample_n)
  Idents(obj_down) <- group_col
  deg <- FindMarkers(obj_down, ident.1 = ident1, ident.2 = ident2, test.use = "MAST")
  write.csv(deg, csv_file)
  deg_tbl <- format_deg(deg)
  plot_deg_volcano(deg_tbl, xlim = xlim, ylim = ylim, select_labs = select_labs,
                   lab_size = lab_size, add_connectors = add_connectors)
  invisible(list(object_down = obj_down, deg = deg, deg_tbl = deg_tbl))
}

add_ifn_scores <- function(obj) {
  obj <- AddModuleScore(obj, features = ifn_i_genes, name = "IFNIscores50")
  obj <- AddModuleScore(obj, features = ifn_hallmark_genes, name = "HallmarkIFN")
  obj
}

compare_scores <- function(obj, group_col, score_cols, comparisons) {
  md <- obj@meta.data
  purrr::map_dfr(score_cols, function(score) {
    purrr::map_dfr(comparisons, function(pair) {
      x <- md[md[[group_col]] == pair[1], score]
      y <- md[md[[group_col]] == pair[2], score]
      tibble(
        score = score,
        group1 = pair[1],
        group2 = pair[2],
        mean1 = mean(x, na.rm = TRUE),
        mean2 = mean(y, na.rm = TRUE),
        p_value = wilcox.test(x, y)$p.value
      )
    })
  })
}

make_ifn_heatmap <- function(obj, features, column_split, column_gap_mm = 1, show_colnames = TRUE, fontsize = 12) {
  mat <- AverageExpression(obj, features = features, group.by = "Condition")$RNA
  mat <- t(scale(t(as.matrix(mat))))
  col_fun <- colorRamp2(c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)), c("blue", "white", "red"))
  ComplexHeatmap::pheatmap(
    mat,
    column_split = column_split,
    column_gap = grid::unit(column_gap_mm, "mm"),
    cellwidth = 30,
    heatmap_legend_param = list(title = "Average Score (row scaled)"),
    col = col_fun,
    cluster_cols = FALSE,
    show_colnames = show_colnames,
    fontsize = fontsize
  )
}
