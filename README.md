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

After strict filtering and doublet removal, **5,089 cells** were retained for downstream analysis.

## Peak Calling

Reproducible accessible chromatin regions were identified using **MACS2** through ArchR.

Pseudo-bulk sample coverages were generated before reproducible peak calling, followed by construction of the PeakMatrix for downstream analysis.

## Dimensionality Reduction and Clustering

Iterative Latent Semantic Indexing (LSI) was performed on the TileMatrix using 30 dimensions.

UMAP was used for visualization, followed by **Harmony** batch correction using sample identity.

Cell populations were clustered using **Louvain graph clustering** at a resolution of `0.8`.

Seven clusters were identified in the final scATAC-seq dataset.

## Gene Activity and Motif Analysis

ArchR GeneScoreMatrix was used to represent gene activity from chromatin accessibility.

Cluster-specific marker genes were identified and visualized using UMAPs, including MAGIC-style smoothing.

TF motif activity was assessed using **chromVAR** through ArchR with the **CIS-BP** motif set.

Variable TF motifs were visualized across the UMAP and across clusters. SPIB_336 and SPI1_322 were among the most variable motifs examined.

## scRNA-seq Integration

The scATAC-seq dataset was integrated with an scRNA-seq reference using ArchR gene integration.

RNA-derived cell-type labels were projected onto the ATAC dataset, and a confusion matrix was used to compare ATAC clusters with RNA-derived labels.

Gene activity and gene expression were also compared using matrix correlation.

The integration showed good but imperfect agreement between ATAC-derived clusters and RNA-derived cell types.

## Peak-to-Gene Linkage and Regulatory Analysis

Peak-to-gene links were identified to associate accessible regulatory regions with genes.

Co-accessibility analysis was performed to identify correlated accessible regions and potential regulatory interactions around selected marker genes.

Regulatory regions were visualized using ArchR browser tracks.

Potential distal regulatory elements were identified for selected marker genes based on co-accessibility with their transcription start sites.

## COVID-19 vs Healthy Differential Accessibility

Differential accessibility analysis was performed between COVID-19 and healthy samples using the PeakMatrix.

COVID-19- and healthy-associated accessible peaks were evaluated using FDR and log2 fold-change thresholds.

MA and volcano plots were generated to visualize differential chromatin accessibility.

Most peaks were centred around log2 fold change ≈ 0, with only a very small number of statistically significant differential peaks. The significant changes detected were down-regulated in COVID-19.

## TF Motif Enrichment and Footprinting

TF motif enrichment was assessed for COVID-19- and Healthy-associated accessible peaks using matched background peaks.

No robust set of Healthy-specific or COVID-19-specific peaks passed the significance thresholds required for confident condition-specific TF motif enrichment.

Condition-associated motif activity was visualized using motif activity UMAPs and cluster-level summaries.

TF footprinting was **attempted** for condition-associated motifs. However, because no significantly enriched condition-specific TF motifs were identified, a meaningful condition-specific top-three-motif footprinting result could not be obtained.

Tn5 sequence bias correction was applied during the footprinting analysis to reduce sequence-dependent insertion bias.

## Key Findings

- Quality control metrics were consistent across the four samples, with no clear low-quality outlier sample identified.
- After strict filtering and doublet removal, **5,089 cells** were retained for downstream analysis.
- Reproducible accessible chromatin regions were identified using MACS2.
- Iterative LSI, UMAP, Harmony, and Louvain clustering were used to characterize chromatin accessibility patterns.
- Seven scATAC-seq clusters were identified in the final dataset.
- Cluster-specific gene activity patterns were identified, with MAGIC used to visualize smoothed gene-activity patterns.
- TF motif activity was assessed using chromVAR-based deviation scores, with SPIB_336 and SPI1_322 among the most variable motifs examined.
- scRNA-seq integration enabled comparison of chromatin accessibility with gene expression and showed good but imperfect agreement between ATAC clusters and RNA-derived cell types.
- Differential accessibility between COVID-19 and Healthy samples was limited, with most peaks near log2FC ≈ 0 and only a small number of significant changes.
- No robust condition-specific TF motif enrichment was detected.
- Because significant condition-specific TF motifs were not identified, meaningful condition-specific TF footprinting could not be obtained.
- Co-accessibility analysis identified potential distal regulatory elements for selected genes, with candidate enhancers detected for genes including **CD44, CD69, PAX5, CD14, IRF8, and CCR7**.

## Tools and Methods

- **Programming:** R
- **scATAC-seq analysis:** ArchR
- **Peak calling:** MACS2
- **Dimensionality reduction:** Iterative LSI, UMAP
- **Batch correction:** Harmony
- **Clustering:** Louvain graph clustering
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
