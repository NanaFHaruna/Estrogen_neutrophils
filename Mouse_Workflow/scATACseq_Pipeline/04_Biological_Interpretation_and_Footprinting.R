# --- 1. LOAD LIBRARIES ---
library(Signac)
library(Seurat)
library(GenomicRanges)
library(presto)
library(dplyr)
library(AnnotationHub)
library(motifmatchr)
#BiocManager::install("JASPAR2022")
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
library(BSgenome.Mmusculus.UCSC.mm39)
library(EnsDb.Mmusculus.v79) 
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

# --- 2. DATA LOADING & METADATA ---
combined <- readRDS('merge_with_scGLUEannotations2.rds')
DefaultAssay(combined) <- 'peaks'

# Grouping replicates into experimental (E) and vehicle (V) categories
combined$new_dataset <- recode(as.character(combined$dataset),
                               "Erep1" = "E", "Erep2" = "E",
                               "Vehrep1" = "V", "Vehrep2" = "V")

# --- 3. PEAK HARMONIZATION (Chromosome Filtering) ---
# Ensuring peak names match the BSgenome (chr1, chr2...) and removing scaffolds
scatac <- combined
orig_peaks <- granges(scatac[["peaks"]])
seqlevelsStyle(orig_peaks) <- "UCSC"

# Keep standard chromosomes only
peaks_std <- keepStandardChromosomes(orig_peaks, pruning.mode = "coarse")
cat("Kept", length(peaks_std), "peaks out of", length(orig_peaks), "\n")

# Logic to subset the object while maintaining consistency
gr_to_string <- function(gr) paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
common_strings <- intersect(gr_to_string(orig_peaks), gr_to_string(peaks_std))

# Subset the assay to the common peaks
counts_orig <- GetAssayData(scatac, assay = "peaks", slot = "counts")
idx_orig <- match(common_strings, gr_to_string(orig_peaks))
counts_common <- counts_orig[idx_orig, ]
peaks_common <- peaks_std[match(common_strings, gr_to_string(peaks_std))]

# Replace the assay with the clean version
scatac[["peaks"]] <- CreateChromatinAssay(
  counts = counts_common,
  ranges = peaks_common,
  fragments = Fragments(scatac)
)
combined <- scatac # Synchronize back to combined object

# --- 4. ANNOTATION & MOTIF SETUP ---
ah <- AnnotationHub()
ensdb_v110 <- ah[["AH113713"]]
gene.ranges <- GetGRangesFromEnsDb(ensdb = ensdb_v110)
Annotation(combined) <- gene.ranges

# Add Motifs to the object
pfm_list_mm <- getMatrixSet(JASPAR2020, list(species=10090))
combined <- AddMotifs(combined, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm_list_mm)

# --- 5. DIFFERENTIAL ACCESSIBILITY (MPP3 Subset) ---
Idents(combined) <- "new_ID"
MPP3notdiv <- subset(combined, idents = "MPP3-notdiv")
Idents(MPP3notdiv) <- "new_dataset"

da_peaks <- FindMarkers(
  object = MPP3notdiv, ident.1 = 'E', ident.2 = 'V',
  only.pos = FALSE, test.use = 'LR', min.pct = 0.05, latent.vars = 'nCount_peaks'
)

# Annotate the DA peaks
da_peaks_annotated <- ClosestFeature(combined, regions = rownames(da_peaks))
da_peaks_merged <- merge(da_peaks, da_peaks_annotated, by.x = "row.names", by.y = "query_region")
write.csv(da_peaks_merged, 'da_peaks_MPP3notdiv_rep2.csv')

# --- 6. GENE ONTOLOGY (GO) ENRICHMENT ---
# Defining specific gene lists for enrichment based on logFC
open_E <- rownames(da_peaks[da_peaks$avg_log2FC > 2, ])
closest_genes_open_E_HSC <- ClosestFeature(combined, regions = open_E)

E_HSC_ego <- enrichGO(
  gene = closest_genes_open_E_HSC$gene_id,
  keyType = "ENSEMBL", OrgDb = org.Mm.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
)
if(!is.null(E_HSC_ego)) barplot(E_HSC_ego, showCategory = 20)

# --- 7. MOTIF OVERREPRESENTATION ---
# Subset objects as needed for specific analysis steps
MPP3 <- MPP3notdiv 
MPP3 <- RegionStats(MPP3, genome = BSgenome.Mmusculus.UCSC.mm10)

# Match background peaks for GC content
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
enriched.motifs <- FindMotifs(object = MPP3, features = top.da.peak)
write.csv(enriched.motifs, 'MPP3all_enrichedmotif_rep2.csv')
MotifPlot(object = MPP3, motifs = head(rownames(enriched.motifs)))

# --- 8. HSC SUBSET & FOOTPRINTING ---
HSC <- subset(combined, idents = "HSC")

# Visualizing Myl10 region
open_E_HSC <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ]) # Example filter
regions_highlight <- subsetByOverlaps(StringToGRanges(open_E_HSC), LookupGeneCoords(HSC, "Myl10"))

CoveragePlot(
  object = HSC, region = "Myl10", region.highlight = regions_highlight,
  extend.upstream = 1000, extend.downstream = 1000
)

# TF Footprinting
combined <- Footprint(
  object = combined, motif.name = "ESR1",
  genome = BSgenome.Mmusculus.UCSC.mm10
)

p2 <- PlotFootprint(combined, features = c("GATA2", "CEBPA", "EBF1"))
p2 + patchwork::plot_layout(ncol = 1)

# --- 9. CHROMVAR & EXPORT ---
MPP3 <- RunChromVAR(MPP3, genome = BSgenome.Mmusculus.UCSC.mm10)
DimPlot(combined, split.by = "dataset", label = TRUE, ncol = 2)
