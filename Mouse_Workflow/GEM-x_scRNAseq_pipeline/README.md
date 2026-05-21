# Mouse GEM-X Neutrophil scRNA-seq Pipeline

## Overview

This folder contains the scRNA-seq workflow used for preprocessing, hashtag oligonucleotide (HTO) demultiplexing, quality control, integration, clustering, and downstream analysis of mouse GEM-X datasets.

The workflow is implemented in  **R** using the **Seurat** library. The scripts begin with HTO-based sample splitting and quality control, followed by dataset merging, integration, neutrophil compartment separation, and downstream blood and bone marrow (BM) analyses.

---

## Workflow Summary

The pipeline performs the following tasks:

- Perform HTO-based sample demultiplexing and condition separation
- Generate filtered Seurat objects from multiplexed GEM-X experiments
- Perform quality control filtering of scRNA-seq datasets
- Merge GEM-X samples across experimental conditions
- Conduct merge-level exploratory analysis and quality assessment
- Integrate datasets using Seurat integration workflows
- Perform dimensionality reduction and clustering analyses
- Separate neutrophil populations by biological compartment and condition
- Generate blood-specific and BM/cKit-specific downstream analyses
- Perform differential gene expression and comparative analysis
- Generate Venn diagrams and overlap analyses for DEG comparisons

---
