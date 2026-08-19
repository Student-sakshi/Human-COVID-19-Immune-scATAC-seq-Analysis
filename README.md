# Human COVID-19 Immune scATAC-seq Analysis

Single-cell chromatin accessibility analysis of immune-cell populations from COVID-19 and healthy samples.

## Overview

This repository contains an ArchR-based analysis of human single-cell ATAC sequencing (scATAC-seq) data from immune-cell samples from COVID-19 and healthy individuals.

The analysis covers:

- Quality control and preprocessing
- Doublet detection and filtering
- Peak calling
- LSI and UMAP dimensionality reduction
- Batch correction and clustering
- Cluster-specific accessible regions
- Gene activity analysis
- Transcription-factor motif activity
- Integration with scRNA-seq data
- Peak-to-gene linkage
- COVID-19 versus healthy differential accessibility
- TF motif enrichment and footprinting
- Co-accessibility and regulatory element analysis

The complete analysis and results are documented in the project report.

## Dataset and Sample Metadata

The dataset contains four immune-cell scATAC-seq samples:

| Sample | Condition |
|---|---|
| `ATAC_555_2_fragments_fragments` | COVID-19 |
| `ATAC_557_fragments_fragments` | COVID-19 |
| `ATAC_EV08_fragments_fragments` | Healthy |
| `ATAC_HIP02_frozen_fragments_fragments` | Healthy |

The original input files are not included in this repository.

## Analysis Workflow

The analysis followed this workflow:

**scATAC-seq Fragments → Quality Control → Doublet Detection → TSS/Fragment Filtering → Peak Calling → Iterative LSI → UMAP → Harmony Batch Correction → Clustering → Peak/Gene Activity Analysis → TF Motif Analysis → scRNA-seq Integration → Peak-to-Gene Linkage → Differential Accessibility → Motif Enrichment → TF Footprinting → Co-accessibility**

## Quality Control and Preprocessing

Initial quality control was performed using fragment counts and TSS enrichment.

Cells were retained using the following thresholds:

- Minimum TSS enrichment: `5`
- Minimum fragments: `550`
- Doublet enrichment: `< 2`
- Final TSS enrichment: `> 8`
- Final fragments: `> 3000`

Fragment-size distributions, TSS enrichment, and fragment count versus TSS enrichment were examined across samples.

## Peak Calling

Reproducible accessible chromatin regions were identified using **MACS2** through ArchR.

Pseudo-bulk sample coverages were generated before reproducible peak calling, followed by construction of the PeakMatrix for downstream analysis.

## Dimensionality Reduction and Clustering

Iterative Latent Semantic Indexing (LSI) was performed on the TileMatrix using 30 dimensions.

UMAP was used for visualization, followed by **Harmony** batch correction using sample identity.

Cell populations were clustered using Seurat-based graph clustering at a resolution of `0.8`.

## Gene Activity and Motif Analysis

ArchR GeneScoreMatrix was used to represent gene activity from chromatin accessibility.

Cluster-specific marker genes were identified and visualized using UMAPs, including MAGIC-style imputation.

TF motif activity was assessed using **chromVAR** through ArchR with the **CIS-BP** motif set.

Variable TF motifs were visualized across the UMAP and across clusters.

## scRNA-seq Integration

The scATAC-seq dataset was integrated with an scRNA-seq reference using ArchR gene integration.

RNA-derived cell-type labels were projected onto the ATAC dataset, and a confusion matrix was used to compare ATAC clusters with RNA-derived labels.

Gene activity and gene expression were also compared using matrix correlation.

## Peak-to-Gene Linkage and Regulatory Analysis

Peak-to-gene links were identified to associate accessible regulatory regions with genes.

Co-accessibility analysis was performed to identify correlated accessible regions and potential regulatory interactions around selected marker genes.

Regulatory regions were visualized using ArchR browser tracks.

## COVID-19 vs Healthy Differential Accessibility

Differential accessibility analysis was performed between COVID-19 and healthy samples using the PeakMatrix.

COVID-19- and healthy-associated accessible peaks were identified using FDR and log2 fold-change thresholds.

MA and volcano plots were generated to visualize differential chromatin accessibility.

## TF Motif Enrichment and Footprinting

TF motif enrichment was performed for COVID-19- and healthy-associated accessible peaks using matched background peaks.

Condition-associated motifs were visualized using motif activity UMAPs and cluster-level heatmaps.

TF footprinting was performed for significant condition-associated motifs to examine local accessibility patterns around motif positions.

## Key Findings

- Quality control and doublet filtering were used to obtain high-quality scATAC-seq profiles.
- Reproducible accessible chromatin regions were identified using MACS2.
- LSI, UMAP, Harmony, and clustering were used to characterize chromatin accessibility patterns.
- Cluster-specific accessible regions and gene activity profiles were identified.
- TF motif activity was assessed using chromVAR-based deviation scores.
- scRNA-seq integration provided RNA-derived cell-type labels for the ATAC clusters.
- Peak-to-gene linkage and co-accessibility analysis were used to investigate potential regulatory relationships.
- COVID-19 versus healthy analysis identified condition-associated differences in chromatin accessibility.
- TF motif enrichment and footprinting were used to investigate potential transcription-factor regulatory activity.

## Tools and Methods

- **Programming:** R
- **scATAC-seq analysis:** ArchR
- **Peak calling:** MACS2
- **Dimensionality reduction:** Iterative LSI, UMAP
- **Batch correction:** Harmony
- **Clustering:** Seurat graph clustering
- **Gene activity:** ArchR GeneScoreMatrix
- **TF motif activity:** chromVAR / CIS-BP
- **scRNA-seq integration:** ArchR gene integration
- **Differential accessibility:** ArchR
- **Peak-to-gene linkage:** ArchR
- **Co-accessibility:** ArchR
- **TF footprinting:** ArchR

## Repository Contents

The repository contains:

- Analysis script(s)
- Project report
- Selected analysis figures
- Analysis results

The complete methodology and figures are described in the project report.

## Report

**Human COVID-19 Immune scATAC-seq Analysis — Project Report**

## Author

**Sakshi Parate**  
M.Sc. Bioinformatics, Saarland University
