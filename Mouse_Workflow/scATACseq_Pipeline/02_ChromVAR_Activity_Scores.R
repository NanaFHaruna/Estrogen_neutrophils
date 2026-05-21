# --- 1. SETUP & LIBRARIES ---
library(Signac)
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


# Load the merged object
combined <- readRDS("merge_with_scGLUEannotations2.rds")
DefaultAssay(combined) <- "peaks"

# --- 2. CHROMOSOME FILTERING ---
# Ensure only standard chromosomes (chr1-19, X, Y, M) are present for ChromVAR
scatac <- combined
orig_peaks <- granges(scatac)
seqlevelsStyle(orig_peaks) <- "UCSC"
peaks_std <- keepStandardChromosomes(orig_peaks, pruning.mode = "coarse")

# Subset object to standard peaks
common_peaks <- intersect(gr_to_string(orig_peaks), gr_to_string(peaks_std))
scatac <- subset(scatac, features = common_peaks)

# Update Chromatin Assay ranges to match standard style
new_assay <- CreateChromatinAssay(
  counts = GetAssayData(scatac, assay = "peaks", slot = "counts"),
  ranges = peaks_std[match(common_peaks, gr_to_string(peaks_std))],
  fragments = Fragments(scatac)
)
scatac[["peaks"]] <- new_assay

# --- 3. RUN CHROMVAR ---
# Fetch JASPAR 2022 motifs
pfm_list_mm <- getMatrixSet(
  x = JASPAR2022, 
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motifs and calculate per-cell deviation scores
scatac <- AddMotifs(scatac, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm_list_mm)
scatac <- RunChromVAR(scatac, genome = BSgenome.Mmusculus.UCSC.mm10, assay = "peaks")

# Save the ChromVAR-enabled object
saveRDS(scatac, 'chromvarall_1.rds')

# --- 4. DIFFERENTIAL MOTIF ACTIVITY (MPP3 Example) ---
# Switch to ChromVAR assay for differential testing
DefaultAssay(scatac) <- 'chromvar'
Idents(scatac) <- "new_ID"

# Subset for specific cell types (Merging all MPP3 states)
MPP3 <- subset(scatac, idents = c("MPP3-notdiv", "MPP3-div", "MPP3-div2"))
Idents(MPP3) <- "dataset"

# Find motifs with differential activity between E and Veh
# Note: ChromVAR markers use 'avg_diff' instead of 'avg_log2FC'
da_motifs <- FindMarkers(
  object = MPP3,
  ident.1 = 'Erep2',
  ident.2 = 'Vehrep2',
  test.use = 'wilcox',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

# Visualize the top results
MotifPlot(
  object = MPP3,
  motifs = head(rownames(da_motifs)),
  assay = 'peaks'
)