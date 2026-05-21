# Mouse scATAC-seq Regulatory Analysis Pipeline

## Overview

This folder contains the complete single-cell ATAC sequencing (scATAC-seq) workflow used for chromatin accessibility integration, regulatory activity analysis, differential accessibility testing, motif enrichment analysis, transcription factor footprinting, and downstream visualization of mouse scATAC-seq datasets.

The workflow is implemented in **R** using **Seurat/Signac-based** analysis frameworks and is organized into modular analysis stages beginning with dataset integration and merging, followed by chromatin accessibility quantification, motif activity analysis, differential accessibility testing, biological interpretation, and visualization.


---

## Workflow Summary

The pipeline performs the following tasks:

- Integrate and merge scATAC-seq datasets across samples and conditions
- Generate unified chromatin accessibility objects
- Compute chromVAR transcription factor activity scores
- Identify differentially accessible chromatin regions
- Perform motif enrichment and regulatory analysis
- Conduct transcription factor footprinting analysis
- Interpret biological regulatory programs across conditions
- Generate publication-quality visualizations and volcano plots

---
