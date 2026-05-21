# --- 1. LIBRARIES & SETUP ---library(Seurat)
library(magrittr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(SingleR)
library(scRNAseq)
library(scater)
library(magrittr)
library(knitr)
library(stringr)
library(reshape2)
library(celldex)
library(ggrepel)
library(RColorBrewer)
library(scCustomize)

# Global variables
project <- "NIAMS-34"
sample_list <- c("WA_036_merge_outs")
sample_list2 <- c("WA_036")
qcDir <- "/data/Amblerwg_Schaughency/will/Gonadectomy_final/Sample_QC" 
dataDir <- "/data/Amblerwg_Schaughency/rawdata/will"
hto_thresh <- 10000

# Helper Function: Calculate 3-MAD Thresholds in log-space
calc_thresh <- function(x) {
  m <- median(log1p(x))
  d <- mad(log1p(x))
  return(c(low = exp(m - 3*d) - 1, high = exp(m + 3*d) - 1))
}

# --- 2. SAMPLE PROCESSING LOOP ---
for (i in 1:length(sample_list)) {
  sample <- sample_list[i] 
  sample2 <- sample_list2[i]
  print(paste("Processing Sample:", sample2))
  
  # Load Data
  data_path <- paste0(dataDir, "/", sample, "/per_sample_outs/", sample2, "/count/sample_filtered_feature_bc_matrix/")
  data <- Read10X(data_path)
  
  # Create Seurat Object
  if (inherits(data, "dgCMatrix")) {
    seur <- CreateSeuratObject(counts = data, project = project) 
  } else {
    seur <- CreateSeuratObject(counts = data$'Gene Expression', project = project) 
    # Add HTO Assay
    hto_counts <- data$`Multiplexing Capture`
    hto_keep <- hto_counts[rowSums(hto_counts) > hto_thresh, , drop = FALSE]
    seur[['HTO']] <- CreateAssayObject(counts = hto_keep)
  }
  
  # Directories
  workdir <- file.path(qcDir, sample)
  if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
  setwd(workdir)
  
  # QC Metrics (Mouse pattern)
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^mt-")
  
  # Statistical Thresholding
  t_feat  <- calc_thresh(seur$nFeature_RNA)
  t_count <- calc_thresh(seur$nCount_RNA)
  t_mito  <- calc_thresh(seur$percent.mito)
  
  # Identify cells to remove
  cells_remove <- unique(c(
    colnames(seur)[which(seur$nFeature_RNA < t_feat['low'] | seur$nFeature_RNA > t_feat['high'])],
    colnames(seur)[which(seur$nCount_RNA < t_count['low'] | seur$nCount_RNA > t_count['high'])],
    colnames(seur)[which(seur$percent.mito > t_mito['high'])]
  ))
  
  # Pre-Filter Plots
  p1 <- FeatureScatter(seur, "nCount_RNA", "percent.mito")
  p2 <- FeatureScatter(seur, "nCount_RNA", "nFeature_RNA")
  ggsave("PreFilter_Scatter.png", p1 + p2, width = 10, height = 5)
  
  # Apply Filter
  seur <- subset(seur, cells = cells_remove, invert = TRUE)
  
  # --- 3. DEMULTIPLEXING (If HTO present) ---
  if ("HTO" %in% Assays(seur)) {
    seur <- NormalizeData(seur, assay = "HTO", normalization.method = "CLR")
    seur <- ScaleData(seur, assay = "HTO")
    seur <- HTODemux(seur, assay = "HTO", positive.quantile = 0.99)
    seur <- MULTIseqDemux(seur, assay = "HTO", quantile = 0.99)
    
    # Remove Doublets and Negatives
    seur <- subset(seur, subset = HTO_classification.global == "Singlet")
    write.csv(table(seur$HTO_classification.global), 'doublet_count.csv')
  }
  
  # --- 4. NORMALIZATION & CLUSTERING ---
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur, nfeatures = 2000)
  seur <- ScaleData(seur, features = rownames(seur))
  seur <- RunPCA(seur, npcs = 50)
  
  png("ElbowPlot.png", width=8, height=6, units='in', res=300)
  print(ElbowPlot(seur, ndims = 50))
  dev.off()
  
  seur <- FindNeighbors(seur, dims = 1:30)
  seur <- FindClusters(seur, resolution = 0.8, algorithm = 3)
  seur <- RunUMAP(seur, dims = 1:30)
  
  # Visualization
  png("UMAP_RNA.png", width=1800, height=1600, res = 300)
  print(DimPlot(seur, label = TRUE) + ggtitle(paste("RNA clusters:", sample2)))
  dev.off()
  
  # --- 5. AUTOMATED ANNOTATION (SingleR) ---
  sce <- as.SingleCellExperiment(seur)
  
  # References
  refs <- list(Blueprint = BlueprintEncodeData(), Nova = NovershternHematopoieticData())
  
  for (ref_name in names(refs)) {
    for (level in c("label.main", "label.fine")) {
      task_name <- paste0(ref_name, "_", gsub("label.", "", level))
      message("Running SingleR: ", task_name)
      
      pred <- SingleR(test = sce, ref = refs[[ref_name]], labels = refs[[ref_name]][[level]])
      seur[[task_name]] <- pred$pruned.labels
      
      # Save individual prediction objects and plots
      saveRDS(pred, paste0("pred_", task_name, ".rds"))
      pdf(paste0("UMAP_", task_name, ".pdf"))
      print(DimPlot(seur, group.by = task_name, label = TRUE) + ggtitle(task_name))
      dev.off()
    }
  }
  
  # Final Save per sample
  saveRDS(seur, paste0(sample2, "_Processed_Annotated.rds"))
  writeLines(capture.output(sessionInfo()), paste0(sample2, "_sessionInfo.txt"))
}
