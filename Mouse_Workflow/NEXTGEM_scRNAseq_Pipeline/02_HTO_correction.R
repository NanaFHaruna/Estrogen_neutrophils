# --- 1. Load Libraries ---
library(Seurat)
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
library(celldex, lib.loc = "/data/Amblerwg_Schaughency/R_libs")
library(ggrepel)
library(RColorBrewer)
library(scCustomize, lib.loc = "/data/Amblerwg_Schaughency/R_libs")
library(stringr)
library(tibble)

# --- 2. SETUP MASTER LOOKUP TABLE ---
# Load the master metadata file once
master_metadata <- read.csv("mouse_metadata_11_27_2023.csv")
master_metadata <- master_metadata[, -3] # Remove column 3 as per original script

# Standardize the join key to "Sample.HTO"
# Note: Original code had mismatching separators. Standardizing to "." ensures success.
master_metadata$correctedSample.CMO <- str_c(master_metadata$Sample, master_metadata$HTO, sep = ".")

# --- 3. DEFINE SAMPLES AND PATHS ---
# Add all sample IDs to this list
samples <- c("WA_011", "WA_013", "WA_017", "WA_018", "WA_033", "WA_034", "WA_035", "WA_036", "WA_060", "WA_061")
base_path <- "/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Sample_QC_Gonadectomy/"

# --- 4. PROCESSING LOOP ---
for (s_id in samples) {
  message("--- Processing Sample: ", s_id, " ---")
  
  # Handle the different folder naming conventions (_merge vs _merge_outs)
  folder_suffix <- if (s_id %in% c("WA_060", "WA_061")) "_merge" else "_merge_outs"
  file_path <- paste0(base_path, s_id, folder_suffix, "/seur_cluster.rds")
  
  if (!file.exists(file_path)) {
    warning("File not found for ", s_id, ". Skipping."); next
  }
  
  # Load Seurat object
  seur <- readRDS(file_path)
  seur$orig.ident <- s_id
  
  # Create the Seurat-side join key (e.g., WA_011.HTO1)
  seur$correctedSample.CMO <- str_c(seur$orig.ident, seur$hash.ID, sep = ".")
  
  # Join metadata using tidyverse logic
  # Using 'rownames_to_column' ensures barcodes and metadata stay perfectly aligned
  updated_meta <- seur@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(master_metadata, by = "correctedSample.CMO") %>%
    column_to_rownames("barcode")
  
  # Verification: Check for successful matches
  matches <- sum(!is.na(updated_meta$Condition))
  message("Matched ", matches, " cells out of ", ncol(seur))
  
  # Select only the specific columns requested for the final object
  final_meta_columns <- updated_meta[, c("correctedSample.CMO", "CellType", "Batch", "Sex", "Condition")]
  
  # Add updated metadata back to Seurat object
  seur <- AddMetaData(seur, metadata = final_meta_columns)
  
  # Save individual corrected object
  saveRDS(seur, paste0('correct_HTO_', s_id, '.rds'))
  
  # Free up memory for next iteration
  rm(seur, updated_meta, final_meta_columns)
  gc() 
}

message("HTO Correction Complete.")
