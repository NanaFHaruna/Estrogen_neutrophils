# AIS: Acute Ischemic Stroke scRNA-seq Integration and Analysis Pipeline

## Overview
This folder contains the complete single-cell RNA sequencing (scRNA-seq) workflow used for preprocessing, hashtag oligonucleotide (HTO) correction, merging, batch analysis, integration, clustering, and downstream analysis of mouse NextGEM datasets.

The workflow is implemented **R** using the **Seurat** library.  The scripts detail the stages of sample quality control and HTO correction, followed by dataset merging, batch assessment, integration, clustering, and downstream transcriptional analysis.



---

## Workflow Summary

The pipeline performs the following tasks:

- Perform quality control filtering of raw NextGEM scRNA-seq datasets
- Correct and harmonize HTO-based sample assignments
- Merge datasets across samples and experimental conditions
- Assess batch structure and technical variation
- Integrate datasets using Seurat integration workflows
- Generate dimensionality reduction embeddings, including PCA and UMAP
- Perform clustering and subclustering analyses
- Conduct downstream transcriptional and comparative analyses
- Generate publication-quality visualizations and DEG outputs

---
