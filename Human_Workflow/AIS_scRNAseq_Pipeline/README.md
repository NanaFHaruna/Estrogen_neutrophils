# SLE: Neutrophil and Low-Density Granulocyte scRNA-seq Pipeline

## Overview

This folder contains the complete single-cell RNA sequencing (scRNA-seq) workflow used for integration, preprocessing, clustering, and downstream analysis of SLE neutrophil and PBMC datasets.

The workflow is implemented in **R** using the **Seurat** library. The scripts are organized into sequential analysis stages beginning with sample merging and object preparation, followed by integration, LDG extraction, metadata annotation, clustering, visualization, and downstream biological analysis.

---

## Workflow Summary

The pipeline performs the following major tasks:

- Merge and preprocess raw Seurat objects from SLE cohorts
- Independently merge neutrophil and PBMC datasets
- Integrate datasets across conditions and batches using Seurat integration workflows
- Identify and subset low-density granulocyte (LDG) populations from PBMC-derived cells
- Combine LDGs with neutrophil datasets for integrated analysis
- Add and harmonize metadata annotations across samples and conditions
- Perform dimensionality reduction, clustering, visualization, and downstream transcriptional analysis of neutrophil and LDG populations

---
