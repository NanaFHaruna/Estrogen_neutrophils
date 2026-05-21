# AIS: Acute Ischemic Stroke scRNA-seq Integration and Analysis Pipeline

## Overview

This folder contains the scRNA-seq workflow used for preprocessing, quality control, integration, clustering, visualization, and downstream analysis of Acute Ischemic Stroke (AIS) datasets.

The workflow is implemented in **R** using the **Seurat** library. The scripts detail the stages of raw sample quality control and dataset merging, integration, dimensionality reduction, clustering, visualization, and downstream transcriptional analysis.


---

## Workflow Summary

The pipeline performs the following tasks:

- Quality control filtering of raw AIS scRNA-seq samples
- Preprocess individual Seurat objects
- Merge AIS datasets across samples and experimental conditions
- Perform downstream merge-level exploratory analysis and quality assessment
- Integrate datasets using Seurat integration workflows to reduce batch effects
- Generate dimensionality reduction embeddings, including PCA and UMAP
- Perform clustering and subclustering analyses
- Generate visualization outputs for AIS immune populations
- Downstream transcriptional and comparative analyses

---
