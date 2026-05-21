# ---- Load Libraries ----
library(Seurat, lib.loc = "/data/Amblerwg_Schaughency/R4.3.0_lib")
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
library(patchwork)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)


# ---- 1. Initial Subsetting & Sub-Clustering ----
# Note: Loading the master integrated object and subsetting on clusters 
# identified as neutrophils and progenitors.
seur <- readRDS("/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_next_gem/integrated_SeuratInt_next_gem_clustered.rds")

Neu <- subset(seur, idents = c(23,31,20,28,18,10,45,15,9,8,7,19,1))

Neu <- FindVariableFeatures(object = Neu) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)

# Visualizing initial clusters
DimPlot(object = Neu, reduction = "umap", label = T)

# ---- 2. Identify and Remove Contaminants ----
# Note: Running FindMarkers on a downsampled set to identify "junk" clusters (10, 12, 13, 14).
down_neu <- subset(Neu, downsample = 5000)
cluster12.markers <- FindMarkers(down_neu, ident.1 = 12)
cluster10.markers <- FindMarkers(down_neu, ident.1 = 10)
cluster14.markers <- FindMarkers(down_neu, ident.1 = 14)
cluster13.markers <- FindMarkers(down_neu, ident.1 = 13)

write.csv(cluster12.markers, 'clustered_Neu+prog_cluster12.markers.csv')
write.csv(cluster10.markers, 'clustered_Neu+prog_cluster10.markers.csv')
write.csv(cluster14.markers, 'clustered_Neu+prog_cluster14.markers.csv')
write.csv(cluster13.markers, 'clustered_Neu+prog_cluster13.markers.csv')

# Subset on "Clean" Clusters only
Neu_clean <- subset(Neu, idents = c(11,8,7,6,9,2,4,1,3,0,5))
saveRDS(Neu_clean, '/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_next_gem/Clean_neutrophil_all.rds')

# ---- 3. Split by Tissue Compartment ----
DefaultAssay(Neu_clean) <- "RNA"
Idents(Neu_clean) <- "CellType"

bm <- subset(Neu_clean, idents = "BM")
blood_spleen <- subset(Neu_clean, idents = c("Blood", "Spleen"))

# Restore cluster identities for subsets
Idents(Neu_clean) <- "seurat_clusters"
Idents(blood_spleen) <- "seurat_clusters"
Idents(bm) <- "seurat_clusters"

# ---- 4. Mature and Post-Mitotic Refinement ----
# Note: Separating dividing cells from post-mitotic and mature neutrophils.
Idents(blood_spleen) <- "seurat_clusters"
blood_spleen_neu1 <- subset(blood_spleen, idents = c(0,1,2,3,4,5))

blood_spleen_neu1 <- FindVariableFeatures(object = blood_spleen_neu1) %>%
  ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.4) %>% RunUMAP(dims = 1:30)

saveRDS(blood_spleen_neu1, '/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_next_gem/blood+spleen_wo_dividingcells.rds')

# ---- 5. Sex and Condition Splitting (Mature & Post-Mitotic) ----
# Note: Creating tissue-specific objects for Male/Female and re-leveling Conditions.
Idents(blood_spleen_neu1) <- "Sex"
female1 <- subset(blood_spleen_neu1, idents = "Female")
Male1 <- subset(blood_spleen_neu1, idents = "Male")

Idents(female1) <- "CellType"
bloodF1 <- subset(female1, idents = "Blood")
spleenF1 <- subset(female1, idents = "Spleen")

Idents(Male1) <- "CellType"
bloodM1 <- subset(Male1, idents = "Blood")
spleenM1 <- subset(Male1, idents = "Spleen")

bloodF1$Condition <- factor(bloodF1$Condition, levels = c("Sham", "OOP"))
spleenF1$Condition <- factor(spleenF1$Condition, levels = c("Sham", "OOP"))
bloodM1$Condition <- factor(bloodM1$Condition, levels = c("Sham", "ORX"))
spleenM1$Condition <- factor(spleenM1$Condition, levels = c("Sham", "ORX"))

# ---- 6. Downsampling Post-Mitotic Neutrophils ----
# Female Blood
Idents(bloodF1) <- "Condition"
bloodF1 <- subset(bloodF1, downsample = 4343)
Idents(bloodF1) <- "seurat_clusters"

# Female Spleen
Idents(spleenF1) <- "Condition"
spleenF1 <- subset(spleenF1, downsample = 3074)
Idents(spleenF1) <- "seurat_clusters"

# Male Blood
Idents(bloodM1) <- "Condition"
bloodM1 <- subset(bloodM1, downsample = 2314)
Idents(bloodM1) <- "seurat_clusters"

# Male Spleen
Idents(spleenM1) <- "Condition"
spleenM1 <- subset(spleenM1, downsample = 1690)
Idents(spleenM1) <- "seurat_clusters"

# ---- 7. Mature Neutrophil Refinement ----
# Subsetting further on specific mature clusters (0,1,3,5)
Idents(blood_spleen) <- "seurat_clusters"
blood_spleen_neu <- subset(blood_spleen, idents = c(0,1,3,5))

blood_spleen_neu <- FindVariableFeatures(object = blood_spleen_neu) %>%
  ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.4) %>% RunUMAP(dims = 1:30)

saveRDS(blood_spleen_neu, '/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_next_gem/blood+spleen_mature.rds')

# Splitting Mature by Sex/Tissue
DefaultAssay(blood_spleen_neu) <- "RNA"
Idents(blood_spleen_neu) <- "Sex"
female <- subset(blood_spleen_neu, idents = "Female")
Male <- subset(blood_spleen_neu, idents = "Male")

Idents(female) <- "CellType"
bloodF <- subset(female, idents = "Blood")
spleenF <- subset(female, idents = "Spleen")

Idents(Male) <- "CellType"
bloodM <- subset(Male, idents = "Blood")
spleenM <- subset(Male, idents = "Spleen")

# ---- 8. Downsampling Mature Neutrophils ----
# Female Blood Mature
Idents(bloodF) <- "Condition"
bloodF <- subset(bloodF, downsample = 4332)
# Female Spleen Mature
Idents(spleenF) <- "Condition"
spleenF <- subset(spleenF, downsample = 3074)
# Male Blood Mature
Idents(bloodM) <- "Condition"
bloodM <- subset(bloodM, downsample = 2314)
# Male Spleen Mature
Idents(spleenM) <- "Condition"
spleenM <- subset(spleenM, downsample = 1690)

# Re-ordering Condition factors for Mature sets
bloodF$Condition <- factor(bloodF$Condition, levels = c("Sham", "OOP"))
spleenF$Condition <- factor(spleenF$Condition, levels = c("Sham", "OOP"))
bloodM$Condition <- factor(bloodM$Condition, levels = c("Sham", "ORX"))
spleenM$Condition <- factor(spleenM$Condition, levels = c("Sham", "ORX"))

# ---- 9. Bone Marrow (BM) Processing ----
Idents(bm) <- "Sex"
femalebm <- subset(bm, idents = "Female")
malebm <- subset(bm, idents = "Male")

# Female BM Downsampling
femalebm$Condition <- factor(femalebm$Condition, levels = c("Sham", "OOP"))
Idents(femalebm) <- "Condition"
femalebm <- subset(femalebm, downsample = 14815)
Idents(femalebm) <- "seurat_clusters"

# Male BM Downsampling
malebm$Condition <- factor(malebm$Condition, levels = c("Sham", "ORX"))
Idents(malebm) <- "Condition"
malebm <- subset(malebm, downsample = 11714)
Idents(malebm) <- "seurat_clusters"

# ---- 10. Marker Identification and Heatmaps ----
workingseur <- bloodF 
IFN_cluster <- subset(workingseur, idents = 3)
workingseur <- IFN_cluster
DefaultAssay(workingseur) <- "RNA"

# Type 1 IFN Gene Set
ifn_type1 <- c("Cxcl10", "Ifit3", "Rsad2", "Ifit1", "Gbp2", "Ifit3b", "Ifit2", "Ifit1bl1", 
               "Herc6", "Usp18", "Ifi204", "Ifi47", "Isg20", "Pnp", "Cmpk2", "Znfx1", 
               "Irf7", "Slfn5", "Ifi211", "Plac8", "Parp14", "Nt5c3", "Sp100", "Ctss", 
               "Ifit3b", "Gbp5", "Gbp3", "Stat2", "Ifi209", "Fcgr1", "Slfn4", "Tor3a", 
               "Ddx60", "H2-T23", "Socs1" ,"Eif2ak2", "Trim30c", "Oas3", "Oasl1", 
               "Psmb9", "Psmb8", "Acer3", "Xaf1", "Sp110", "Trim30a")

# Type II IFN Gene Set
ifn_type2 <- c("Gbp2", "Ifi47", "Gbp5", "Cd274", "Fcgr4", "Irf1", "Irgm1", "Socs1", 
               "Gbp7", "Wfdc17", "Igtp", "Ifi204", "Rsad2", "Parp14", "Psmb9", "Herc6", 
               "Nampt", "H2-T23", "Isg15", "Slc31a1", "Samhd1", "Psmb8", "Tnfaip2", 
               "Casp4", "Stat1", "Ifitm3", "Oasl2", "H2-D1", "Grina")

# TNF Gene Set
tnf_pathway <- c("Sod2", "Lcn2", "Chil1", "Cd14", "Cybb", "Bcl2a1a", "Prdx5", "Wfdc21", 
                 "Marcksl1", "AA467197", "Icam1", "Ehd1", "N4bp1", "Dusp16", "Il1rn", 
                 "Mapkapk2", "Cyfip1", "Upp1", "Rab20", "Ikbke", "Acod1", "Fth1", 
                 "Plac8", "Ctsb", "Cxcl2", "C3", "Vasp", "Nfkbia", "Samsn1", "Fos", 
                 "Cotl1", "Mmp9", "Taldo1", "Fxyd5", "Stk17b", "Gsn", "Fgl2", "Rassf3", 
                 "Dhrs7", "Cd300ld", "Hist1h1c", "Rgs2", "Ccl6")

# IL1b Gene Set
il1b_pathway <- c("Cxcl3", "Wfdc17", "Ifitm6", "Lcn2", "Thbs1", "Cd177", "Rnasel", "Il1r2", 
                  "Id1", "Xbp1", "Tarm1", "Cd14", "Ier3", "Plaur", "Fcgr2b", "Pnp", 
                  "Wfdc21", "Basp1", "Lipg", "Cish", "Glipr2", "Ccl6", "Rab44", "Upp1", 
                  "Fth1", "Ltb4r1", "Lrg1", "Ptpn1", "Retnlg", "Cd33", "Srgn", "App", 
                  "Cd300lf", "Slfn4", "Crispld2", "Csf2rb", "Cxcl2", "S100a9", "Lilr4b", 
                  "Cxcr2", "S100a8")

# ---- 11. Matrix Generation (Choose your list here) ----
# To generate the heatmap for a specific pathway, replace 'ifn_type1' with 
# 'ifn_type2', 'tnf_pathway', or 'il1b_pathway'.
mat_data <- AverageExpression(workingseur, features = ifn_type1, group.by = 'orig.ident')

# Extract counts from RNA assay
mat <- as.matrix(mat_data$RNA)

# Perform row scaling (z-score calculation)
scaled_mat <- t(scale(t(mat)))


# ---- 11. Complex Heatmaps ----
col_fun = colorRamp2(c(min(scaled_mat), 0, max(scaled_mat)), c("blue", "white", "red"))

# Female Plot
ComplexHeatmap::pheatmap(scaled_mat, column_split=c(4,1,5,2,6,3), column_gap=unit(2,"mm"), 
                         column_title_gp=grid::gpar(fill=c("#56B4E9","#56B4E9","#56B4E9", "red", "red", "red")), 
                         heatmap_legend_param=list(title='Average Score (row scaled)'), 
                         cellwidth = 35, col = col_fun, cluster_cols = F, fontsize_row = 14)

# ---- 12. DEG and EnhancedVolcano ----
Idents(workingseur) <- "Condition"
DEGtest <- FindMarkers(workingseur, ident.1 = "Sham", ident.2 = "ORX", test.use = 'MAST')

# Volcano with complete list of gene labels as requested
EnhancedVolcano(DEGtest, lab = rownames(DEGtest), 
                selectLab = c("Junb", "Fosb", "Il1b", "Nr4a1", "Ptgs2", "Il1b", "Acod1", "Tnfaip3", "Tnfaip2", "H2-Q7", "H2-K1", "Cd80", "Wfdc17", "Clec4d", "Camp", "Ngp", "Ifi27la2", "Ltf", "Slfn5", "Oasl2", "Slfn4", "Rsad2", "Slfn1", "Ifitm3", "Trim30a", "Slfn2", "Ly6e", "Oasl1", "Irf2", "Isg15"), 
                x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.1, FCcutoff = .25, 
                xlim = c(-1.5,1.5), ylim = c(0, 65), max.overlaps = Inf, 
                drawConnectors = T, labSize = 8)




