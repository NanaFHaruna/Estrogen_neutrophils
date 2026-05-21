# --- 1. SETUP & GLOBAL ANNOTATIONS ---
library(Signac)
library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load object and set default assay
atac <- combined
DefaultAssay(atac) <- "peaks"

# Re-group replicates for differential testing
atac$new_dataset <- recode(as.character(atac$dataset),
                           "Erep1" = "E", "Erep2" = "E",
                           "Vehrep1" = "V", "Vehrep2" = "V")

# --- 2. DEFINE REUSABLE FUNCTIONS ---

# Function A: DA Analysis & Gene Mapping
analyze_subset_da <- function(seurat_obj, subset_name, idents_list) {
  # Subset and find markers
  Idents(seurat_obj) <- "new_ID"
  sub_obj <- subset(seurat_obj, idents = idents_list)
  Idents(sub_obj) <- "new_dataset"
  
  da_peaks <- FindMarkers(sub_obj, ident.1 = 'E', ident.2 = 'V', 
                          test.use = 'wilcox', min.pct = 0.1)
  
  da_peaks_df <- rownames_to_column(da_peaks, var = "query_region")
  
  # Split into Open in E vs Open in Veh and annotate
  for (cond in c("E", "Veh")) {
    log_filter <- if(cond == "E") da_peaks_df$avg_log2FC > 2 else da_peaks_df$avg_log2FC < -2
    subset_peaks <- da_peaks_df[log_filter, ]
    
    if(nrow(subset_peaks) > 0) {
      ann <- ClosestFeature(seurat_obj, regions = subset_peaks$query_region, annotation = gene.ranges)
      merged <- inner_join(subset_peaks, ann, by = "query_region")
      # Flatten for CSV
      flat <- data.frame(lapply(merged, as.character), stringsAsFactors = FALSE)
      write.csv(flat, paste0("20260312_", subset_name, cond, "_da.csv"), row.names = FALSE)
    }
  }
  return(list(obj = sub_obj, da = da_peaks_df))
}

# Function B: Motif Enrichment
run_motif_enrichment <- function(seurat_subset, da_df, subset_name) {
  for (cond in c("E2", "Veh")) {
    log_val <- if(cond == "E2") 2 else -2
    # Filter for top high-confidence peaks
    top_peaks <- da_df %>% 
      filter(if(cond == "E2") avg_log2FC > 2 else avg_log2FC < -2) %>%
      filter(p_val < 0.005 & pct.1 > 0.2) %>% 
      pull(query_region)
    
    if(length(top_peaks) >= 10) {
      motifs <- FindMotifs(object = seurat_subset, features = top_peaks)
      write.csv(motifs, paste0(subset_name, "_", cond, "_enrichedmotifs.csv"))
    }
  }
}

# --- 3. EXECUTE PIPELINE FOR ALL CELL TYPES ---

# Run DA Analysis for all subsets
hsc_res         <- analyze_subset_da(atac, "HSC", "HSC")
mpp_res         <- analyze_subset_da(atac, "MPP", "MPP")
mep_res         <- analyze_subset_da(atac, "MEP", "MEP")
mast_baso_res   <- analyze_subset_da(atac, "Mast_baso", "Mast/baso")
mpp3notdiv_res  <- analyze_subset_da(atac, "MPP3notdiv", "MPP3-notdiv")
mpp3div_res     <- analyze_subset_da(atac, "MPP3div", c("MPP3-div", "MPP3-div2"))
mpp4notdiv_res  <- analyze_subset_da(atac, "MPP4notdiv", "MPP4-notdiv")
mpp4div_res     <- analyze_subset_da(atac, "MPP4div", "MPP4-div")

# Run Motif Enrichment for all subsets using the results from above
run_motif_enrichment(hsc_res$obj, hsc_res$da, "HSC")
run_motif_enrichment(mpp_res$obj, mpp_res$da, "MPP")
run_motif_enrichment(mep_res$obj, mep_res$da, "MEP")
run_motif_enrichment(mast_baso_res$obj, mast_baso_res$da, "Mast_baso")
run_motif_enrichment(mpp3notdiv_res$obj, mpp3notdiv_res$da, "MPP3notdiv")
run_motif_enrichment(mpp3div_res$obj, mpp3div_res$da, "MPP3div")
run_motif_enrichment(mpp4notdiv_res$obj, mpp4notdiv_res$da, "MPP4notdiv")
run_motif_enrichment(mpp4div_res$obj, mpp4div_res$da, "MPP4div")

# --- 4. TF FOOTPRINTING (MPP3 focus) ---
DefaultAssay(mpp3notdiv_res$obj) <- "peaks"
mpp3notdiv_res$obj <- Footprint(
  object = mpp3notdiv_res$obj,
  motif.name = c("ESR1", "STAT1", "IRF9"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# Plot Footprints
p <- PlotFootprint(mpp3notdiv_res$obj, features = c("ESR1", "STAT1", "IRF9"))
p + patchwork::plot_layout(ncol = 1)
