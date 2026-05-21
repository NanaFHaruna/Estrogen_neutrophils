# --- 1. LOAD LIBRARIES ---
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

# Optimization: Increase memory limit for large FeatureMatrix calculations
options(future.globals.maxSize = 50 * 1024^3) # 50GB

# --- 2. LOAD PEAK SETS & CREATE UNIFIED PEAKS ---
# Define paths to simplify code
base_path <- "/data/Amblerwg_Schaughency/rawdata/will/250904_VH01191_389_22273LYNX/"

peaks_E1 <- read.table(paste0(base_path, "scATACm_HPSC-E-rep1_outs/peaks.bed"), col.names = c("chr", "start", "end"))
peaks_E2 <- read.table(paste0(base_path, "scATACm_HPSC-E-rep2_outs/peaks.bed"), col.names = c("chr", "start", "end"))
peaks_V1 <- read.table(paste0(base_path, "scATACm_HPSC-Veh-rep1_outs/peaks.bed"), col.names = c("chr", "start", "end"))
peaks_V2 <- read.table(paste0(base_path, "scATACm_HPSC-Veh-rep2_outs/peaks.bed"), col.names = c("chr", "start", "end"))

# Convert to GRanges and reduce to create a single master peak set
gr.Erep1 <- makeGRangesFromDataFrame(peaks_E1)
gr.Erep2 <- makeGRangesFromDataFrame(peaks_E2)
gr.Vehrep1 <- makeGRangesFromDataFrame(peaks_V1)
gr.Vehrep2 <- makeGRangesFromDataFrame(peaks_V2) # FIXED: previously used V1 dataframe

combined.peaks <- reduce(x = c(gr.Erep1, gr.Erep2, gr.Vehrep1, gr.Vehrep2))

# Filter peaks by width (removes outliers/noise)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

# --- 3. LOAD METADATA & FRAGMENTS ---
# Helper function to load and filter metadata
load_filtered_md <- function(path) {
  md <- read.table(file = path, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
  return(md[md$passed_filters > 500, ])
}

md.Erep1 <- load_filtered_md(paste0(base_path, "scATACm_HPSC-E-rep1_outs/singlecell.csv"))
md.Erep2 <- load_filtered_md(paste0(base_path, "scATACm_HPSC-E-rep2_outs/singlecell.csv"))
md.Vehrep1 <- load_filtered_md(paste0(base_path, "scATACm_HPSC-Veh-rep1_outs/singlecell.csv"))
md.Vehrep2 <- load_filtered_md(paste0(base_path, "scATACm_HPSC-Veh-rep2_outs/singlecell.csv"))

# Create Fragment Objects (links the large .tsv.gz files)
frags.Erep1 <- CreateFragmentObject(path = paste0(base_path, "scATACm_HPSC-E-rep1_outs/fragments.tsv.gz"), cells = rownames(md.Erep1))
frags.Erep2 <- CreateFragmentObject(path = paste0(base_path, "scATACm_HPSC-E-rep2_outs/fragments.tsv.gz"), cells = rownames(md.Erep2))
frags.Vehrep1 <- CreateFragmentObject(path = paste0(base_path, "scATACm_HPSC-Veh-rep1_outs/fragments.tsv.gz"), cells = rownames(md.Vehrep1))
frags.Vehrep2 <- CreateFragmentObject(path = paste0(base_path, "scATACm_HPSC-Veh-rep2_outs/fragments.tsv.gz"), cells = rownames(md.Vehrep2))

# --- 4. QUANTIFY & CREATE OBJECTS ---
# Generate the count matrices based on unified peaks
Erep1.counts <- FeatureMatrix(fragments = frags.Erep1, features = combined.peaks, cells = rownames(md.Erep1))
Erep2.counts <- FeatureMatrix(fragments = frags.Erep2, features = combined.peaks, cells = rownames(md.Erep2))
Vehrep1.counts <- FeatureMatrix(fragments = frags.Vehrep1, features = combined.peaks, cells = rownames(md.Vehrep1))
Vehrep2.counts <- FeatureMatrix(fragments = frags.Vehrep2, features = combined.peaks, cells = rownames(md.Vehrep2))

# Create Seurat objects
Erep1 <- CreateSeuratObject(CreateChromatinAssay(Erep1.counts, fragments = frags.Erep1), assay = "ATAC", meta.data = md.Erep1)
Erep2 <- CreateSeuratObject(CreateChromatinAssay(Erep2.counts, fragments = frags.Erep2), assay = "ATAC", meta.data = md.Erep2)
Vehrep1 <- CreateSeuratObject(CreateChromatinAssay(Vehrep1.counts, fragments = frags.Vehrep1), assay = "ATAC", meta.data = md.Vehrep1)
Vehrep2 <- CreateSeuratObject(CreateChromatinAssay(Vehrep2.counts, fragments = frags.Vehrep2), assay = "ATAC", meta.data = md.Vehrep2)

# --- 5. ATTACH scGLUE ANNOTATIONS & MERGE ---
glue_path <- "/data/Amblerwg_Schaughency/projects/1105_ATAC_Mouse/scGlue/"

# FIXED: Corrected the AddMetaData call for Erep2 (previously referred to allData)
Erep1 <- AddMetaData(Erep1, read.csv(paste0(glue_path, 'HPSC-E_rep1_ATAC_obs.csv'), row.names=1))
Erep2 <- AddMetaData(Erep2, read.csv(paste0(glue_path, 'HPSC-E_rep2_ATAC_obs.csv'), row.names=1))
Vehrep1 <- AddMetaData(Vehrep1, read.csv(paste0(glue_path, 'HPSC-Veh_rep1_ATAC_obs.csv'), row.names=1))
Vehrep2 <- AddMetaData(Vehrep2, read.csv(paste0(glue_path, 'HPSC-Veh_rep2_ATAC_obs.csv'), row.names=1))

# Set dataset identity and merge
Erep1$dataset <- 'Erep1'; Erep2$dataset <- 'Erep2'; Vehrep1$dataset <- 'Vehrep1'; Vehrep2$dataset <- 'Vehrep2'

combined <- merge(x = Erep1, y = list(Erep2, Vehrep1, Vehrep2), add.cell.ids = c("Erep1", "Erep2", "Vehrep1", "Vehrep2"))

# --- 6. DIMENSIONALITY REDUCTION & SAVING ---
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi') # Skip 1st dim (depth bias)

saveRDS(combined, 'merge_with_scGLUEannotations2.rds')
