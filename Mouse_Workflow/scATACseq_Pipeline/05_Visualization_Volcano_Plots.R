# ----Load Libraries----
library(Seurat)
library(chromVAR)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)# or your chosen BSgenome
library(JASPAR2022)
library(TFBSTools)
library(GenomicRanges)
library(GenomeInfoDb)
library(dplyr)
library(scCustomize)
library(RSQLite)
#library(magrittr)
library(gridExtra)
library(ggplot2)
library(SingleR)
library(scRNAseq)
library(scater)
library(knitr)
library(stringr)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)

# --- 2. DEFINE REUSABLE PLOTTING FUNCTION ---
# This function handles filtering, axis limits, and styling to ensure all plots look consistent.
plot_chromvar_volcano <- function(file_path, title, select_motifs, x_limits, y_limits, fc_cutoff = 0.5) {
  
  # Load specific CSV result
  data <- read.csv(file_path)
  
  # Filter: keep motifs where at least 10% of cells in either group show activity
  df_filtered <- data %>% 
    filter(pct.1 >= 0.1 | pct.2 >= 0.1)
  
  # Generate Volcano Plot
  EnhancedVolcano(
    df_filtered,
    lab = df_filtered$motif_name,
    x = 'avg_diff',
    y = 'p_val',
    pCutoff = 0.00005,
    FCcutoff = fc_cutoff,
    xlim = x_limits,
    ylim = y_limits,
    title = title,
    subtitle = "Estrogen (E) vs Vehicle (Veh)",
    selectLab = select_motifs,
    xlab = 'Mean Difference',
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    maxoverlapsConnectors = Inf
  )
}

# --- 3. GENERATE ALL VOLCANO PLOTS ---

# 1. HSC
plot_chromvar_volcano(
  file_path = "chromvar_da_HSC_jaspar_with_names.csv",
  title = 'HSC: E vs Veh',
  x_limits = c(-1.5, 3.5), y_limits = c(0, 22), fc_cutoff = 0,
  select_motifs = c("GATA1::TAL1", "STAT1", "Stat4", "IRF7", "GFI1", "Gfi1B", "GATA1", "GATA2", "KLF12", "KLF6", "KLF7")
)

# 2. MPP3 (Not Dividing)
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP3notdiv_jaspar_with_names.csv",
  title = 'MPP3notdiv: E vs Veh',
  x_limits = c(-3.5, 3), y_limits = c(0, 300),
  select_motifs = c("CEBPD", "CEBPE", "CEBPG", "CEBPA", "CEBPB","Spi1", "STAT1::STAT2", "ESR1", "Irf1", "IRF3", "IRF9", "ATF4", "CTCF", "GATA2", "GATA3", "GATA1::TAL1", "IKZF1", "BATF::JUN")
)

# 3. MPP3 (Dividing)
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP3div_jaspar_with_names.csv",
  title = 'MPP3div: E vs Veh',
  x_limits = c(-3.5, 3.5), y_limits = c(0, 350),
  select_motifs = c("CEBPD", "CEBPE", "CEBPG", "CEBPA", "CEBPB","Spi1", "STAT1::STAT2", "ESR1", "Irf1", "IRF3", "IRF9", "ATF4", "CTCF", "GATA2", "GATA3", "Gata3", "GATA5", "GATA1::TAL1", "IKZF1", "BATF::JUN", "Stat2", "IRF4", "IRF8", "IRF5", "Runx1", "GATA4")
)

# 4. MPP4 (Not Dividing)
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP4notdiv_jaspar_with_names.csv",
  title = 'MPP4notdiv: E vs Veh',
  x_limits = c(-1, 1), y_limits = c(0, 50),
  select_motifs = c("ESR1", "NFIL3", "ESR2", "STAT1::STAT2", "Irf1", "Bcl11B")
)

# 5. MPP4 (Dividing)
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP4div_jaspar_with_names.csv",
  title = 'MPP4div: E vs Veh',
  x_limits = c(-3, 4), y_limits = c(0, 110),
  select_motifs = c("ESR1", "NFIL3", "ESR2", "STAT1::STAT2", "Irf1", "Bcl11B", "STAT1:STAT2", "IRF9", "IRF4", "IRF3", "IRF8", "Spi1", "ELF3", "ELF1", "GATA2", "GATA4", "GATA1::TAL1", "KLF6", "KLF2", "Gata3", "TRPS1")
)

# 6. MPP (All)
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP_jaspar_with_names.csv",
  title = 'MPP: E vs Veh',
  x_limits = c(-1, 3), y_limits = c(0, 55), fc_cutoff = 0,
  select_motifs = c("GATA2", "Gata3", "GAT1::TAL1", "GATA4", "SP1", "GATA1", "KLF15", "TRPS1", "MEF2A", "MEF2C", "MEF2D", "MEF2B", "Gfi1B", "NFKB1", "NFKB2", "SCRT2")
)

# 7. MEP
plot_chromvar_volcano(
  file_path = "chromvar_da_MEP_jaspar_with_names.csv",
  title = 'MEP: E vs Veh',
  x_limits = c(-3, 5), y_limits = c(0, 75),
  select_motifs = c("GATA2", "GATA3", "Gata3", "GATA5", "GATA1::TAL1", "GATA4", "TRPS1", "SOX2", "PRDM9", "ETV2::FOXI1", "ZNF189", "KLF7")
)

# 8. Mast_baso
plot_chromvar_volcano(
  file_path = "chromvar_da_Mast_baso_jaspar_with_names.csv",
  title = 'Mast_baso: E vs Veh',
  x_limits = c(-1, 1), y_limits = c(0, 50),
  select_motifs = c("ESR1", "NFIL3", "ESR2", "STAT1::STAT2", "Irf1", "Bcl11B")
)

# 9. MPP4 ALL
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP4_ALL_jaspar_with_names.csv",
  title = 'MPP4: E vs Veh',
  x_limits = c(-3, 4), y_limits = c(0, 100), fc_cutoff = 0,
  select_motifs = c("ESR1", "STAT1::STAT2", "Irf1", "STAT1:STAT2", "IRF9", "IRF4", "IRF3", "IRF8", "GATA2", "GATA4", "GATA1::TAL1", "KLF6", "KLF2", "Gata3")
)

# 10. MPP3 ALL
plot_chromvar_volcano(
  file_path = "chromvar_da_MPP3_ALL_jaspar_with_names.csv",
  title = 'MPP3: E vs Veh',
  x_limits = c(-3, 3), y_limits = c(0, 350), fc_cutoff = 0,
  select_motifs = c("ESR1", "STAT1::STAT2", "CEBPD", "CEBPA", "CEBPB", "CEBPE", "Spi1", "ATF4", "STAT1", "Stat2", "IKZF1", "IRF3", "IRF7", "IRF9", "CTCF", "GATA1::TAL1", "GATA2")
)








