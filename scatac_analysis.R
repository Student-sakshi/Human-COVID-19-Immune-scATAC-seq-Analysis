# ============================================================
# Single-cell ATAC-seq Analysis
# Chromatin accessibility analysis of immune-cell populations from COVID-19 and healthy samples
# Author: Sakshi Parate - M.Sc. Bioinformatics, Saarland University
# Reference genome: hg38
# ============================================================

# 1. Load Packages
required_packages <- c("ArchR","ggplot2","dplyr","Seurat","hexbin")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_packages)>0) stop(paste("The following packages are missing:",paste(missing_packages,collapse=", "),"\nInstall them before running this analysis."))
library(ArchR); library(ggplot2); library(dplyr); library(Seurat); library(grid)
set.seed(42)

# 2. Analysis Configuration
data_dir <- "data"; results_dir <- "results"
dir.create(results_dir,showWarnings=FALSE,recursive=TRUE)
addArchRGenome("hg38"); addArchRThreads(threads=1)
lsi_dims <- 1:30; downstream_lsi_dims <- 2:30
genes_of_interest <- c("CD8A","CD14","GATA1","PAX5","TBX21")
integration_genes <- c("PAX5","EGR1","OLIG2")
marker_genes <- c("CD34","GATA1","PAX5","MS4A1","MME","CD14","MPO","IRF8","CD3D","CD8A","CD4","TBX21","CD3G","NCAM1","FCGR3A","FOXP3","GATA3","RORC","PDCD1","HLA-DRA","CD28","IL2RA","CD69","CD44","CCR7")

# 3. Data Import and Quality Control
input_files <- list.files(data_dir,pattern="\\.tsv\\.gz$",full.names=TRUE)
if(length(input_files)==0) stop("No .tsv.gz input files were found in the data directory.")
sample_names <- sub("\\.tsv\\.gz$","",basename(input_files)); print(sample_names)

ArrowFiles <- createArrowFiles(inputFiles=input_files,sampleNames=sample_names,minTSS=5,minFrags=550,addTileMat=TRUE,addGeneScoreMat=TRUE,force=TRUE)
print(ArrowFiles)

proj <- ArchRProject(ArrowFiles=ArrowFiles,outputDirectory=file.path(results_dir,"ArchRProject"),copyArrows=TRUE)
print(proj)

# 4. Doublet Detection and Filtering
proj <- addDoubletScores(input=proj,k=10,knnMethod="UMAP",LSIMethod=1,force=TRUE)
doublet_enrichment <- proj$DoubletEnrichment
doublet_filter <- doublet_enrichment < 2
filtered_proj <- proj[doublet_filter,]; print(filtered_proj)

# 5. Initial Quality-Control Statistics
initial_qc <- data.frame(metric=c("Number of cells","Median TSS enrichment","Median fragments"),value=c(nCells(proj),median(proj$TSSEnrichment),median(proj$nFrags)))
print(initial_qc)
initial_tile_dimensions <- dim(getMatrixFromProject(proj,useMatrix="TileMatrix",binarize=TRUE)); print(initial_tile_dimensions)
initial_cells_per_sample <- table(proj$Sample); print(initial_cells_per_sample)

fragment_plot <- plotFragmentSizes(proj)
plotPDF(fragment_plot,name="Fragment_Length_Distribution",ArchRProj=proj,addDOC=FALSE)

tss_plot <- plotTSSEnrichment(ArchRProj=proj,groupBy="Sample")
plotPDF(tss_plot,name="TSS_Enrichment_by_Sample",ArchRProj=proj,addDOC=FALSE)

qc_metadata <- as.data.frame(getCellColData(proj,select=c("Sample","nFrags","TSSEnrichment")))
fragment_tss_plots <- lapply(unique(qc_metadata$Sample),function(sample_name) {
  sample_data <- qc_metadata[qc_metadata$Sample==sample_name,]
  ggplot(sample_data,aes(x=nFrags,y=TSSEnrichment))+geom_point(alpha=0.35,size=0.6)+scale_x_log10()+theme_minimal()+labs(title=sample_name,x="Number of fragments",y="TSS enrichment")
})
plotPDF(plotList=fragment_tss_plots,name="Fragments_vs_TSS_by_Sample",ArchRProj=proj,addDOC=FALSE,width=6,height=4)

# 6. Final Quality Filtering
strict_filtered_proj <- filtered_proj[filtered_proj$TSSEnrichment>8 & filtered_proj$nFrags>3000,]
print(strict_filtered_proj)

final_qc <- data.frame(metric=c("Number of cells","Median TSS enrichment","Median fragments"),value=c(nCells(strict_filtered_proj),median(strict_filtered_proj$TSSEnrichment),median(strict_filtered_proj$nFrags)))
print(final_qc)
write.csv(final_qc,file=file.path(results_dir,"final_QC_summary.csv"),row.names=FALSE)

final_cells_per_sample <- table(strict_filtered_proj$Sample); print(final_cells_per_sample)
write.csv(as.data.frame(final_cells_per_sample),file=file.path(results_dir,"final_cells_per_sample.csv"),row.names=FALSE)
final_tile_dimensions <- dim(getMatrixFromProject(strict_filtered_proj,useMatrix="TileMatrix",binarize=TRUE)); print(final_tile_dimensions)

# 7. Peak Calling
proj_peaks <- strict_filtered_proj

pathToMacs2 <- tryCatch(findMacs2(),error=function(e) "")
if(!nzchar(pathToMacs2)||!file.exists(pathToMacs2)) stop("MACS2 was not found. Please install MACS2 and make it available to ArchR.")
print(pathToMacs2)

proj_peaks <- addGroupCoverages(ArchRProj=proj_peaks,groupBy="Sample",force=TRUE)
proj_peaks <- addReproduciblePeakSet(ArchRProj=proj_peaks,groupBy="Sample",peakMethod="Macs2",pathToMacs2=pathToMacs2,force=TRUE)
proj_peaks <- addPeakMatrix(proj_peaks,force=TRUE)
print(getAvailableMatrices(proj_peaks))
peak_matrix <- getMatrixFromProject(proj_peaks,useMatrix="PeakMatrix"); print(peak_matrix)
print(length(getPeakSet(proj_peaks)))

# 8. Dimensionality Reduction
proj_peaks <- addIterativeLSI(ArchRProj=proj_peaks,useMatrix="TileMatrix",name="IterativeLSI",iterations=2,clusterParams=list(resolution=0.2,sampleCells=10000,n.start=10),varFeatures=25000,dimsToUse=lsi_dims,force=TRUE)

proj_peaks <- addUMAP(ArchRProj=proj_peaks,reducedDims="IterativeLSI",name="UMAP",dimsToUse=downstream_lsi_dims,nNeighbors=30,minDist=0.5,metric="cosine",force=TRUE)

umap_sample <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="Sample",embedding="UMAP")
umap_tss <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="TSSEnrichment",embedding="UMAP")
umap_fragments <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="nFrags",embedding="UMAP")
plotPDF(plotList=list(umap_sample,umap_tss,umap_fragments),name="UMAP_QC",ArchRProj=proj_peaks,addDOC=FALSE)

# 9. Batch Correction and Clustering
proj_peaks <- addHarmony(ArchRProj=proj_peaks,reducedDims="IterativeLSI",dimsToUse=downstream_lsi_dims,name="Harmony",groupBy="Sample",force=TRUE)
proj_peaks <- addUMAP(ArchRProj=proj_peaks,reducedDims="Harmony",name="UMAP_Harmony",nNeighbors=30,minDist=0.5,metric="cosine",force=TRUE)

umap_harmony_sample <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="Sample",embedding="UMAP_Harmony")
umap_harmony_tss <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="TSSEnrichment",embedding="UMAP_Harmony")
umap_harmony_fragments <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="nFrags",embedding="UMAP_Harmony")
plotPDF(plotList=list(umap_harmony_sample,umap_harmony_tss,umap_harmony_fragments),name="UMAP_Harmony_QC",ArchRProj=proj_peaks,addDOC=FALSE)

proj_peaks <- addClusters(input=proj_peaks,reducedDims="Harmony",dimsToUse=downstream_lsi_dims,method="Seurat",name="Clusters",resolution=0.8,seed=42,force=TRUE)
cluster_umap <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="Clusters",embedding="UMAP_Harmony")
plotPDF(cluster_umap,name="UMAP_Clusters",ArchRProj=proj_peaks,addDOC=FALSE)

cluster_cell_counts <- table(proj_peaks$Clusters); print(cluster_cell_counts)
cluster_sample_composition <- table(proj_peaks$Clusters,proj_peaks$Sample); print(cluster_sample_composition)
write.csv(as.data.frame(cluster_cell_counts),file=file.path(results_dir,"cluster_cell_counts.csv"),row.names=FALSE)
write.csv(as.data.frame(cluster_sample_composition),file=file.path(results_dir,"cluster_sample_composition.csv"),row.names=FALSE)

# 10. Cluster-Specific Accessible Regions
marker_peaks <- getMarkerFeatures(ArchRProj=proj_peaks,useMatrix="PeakMatrix",groupBy="Clusters",bias=c("TSSEnrichment","nFrags"),testMethod="wilcoxon")
marker_peak_list <- getMarkers(marker_peaks,cutOff="FDR <= 0.01 & Log2FC >= 1"); print(marker_peak_list)

marker_peak_heatmap <- plotMarkerHeatmap(seMarker=marker_peaks,cutOff="FDR <= 0.01 & Log2FC >= 1",transpose=TRUE)
plotPDF(marker_peak_heatmap,name="Cluster_Marker_Peak_Heatmap",ArchRProj=proj_peaks,addDOC=FALSE,width=8,height=8)

gene_tracks <- plotBrowserTrack(ArchRProj=proj_peaks,groupBy="Clusters",geneSymbol=genes_of_interest,upstream=50000,downstream=50000,loops=NULL)
plotPDF(plotList=gene_tracks,name="Marker_Gene_Accessibility_Tracks",ArchRProj=proj_peaks,addDOC=FALSE,width=8,height=8)

# 11. Gene Activity and Marker Genes
if(!"GeneScoreMatrix" %in% getAvailableMatrices(proj_peaks)) proj_peaks <- addGeneScoreMatrix(ArchRProj=proj_peaks,matrixName="GeneScoreMatrix",force=TRUE)
gene_score_matrix <- getMatrixFromProject(ArchRProj=proj_peaks,useMatrix="GeneScoreMatrix"); print(gene_score_matrix)

marker_genes_results <- getMarkerFeatures(ArchRProj=proj_peaks,useMatrix="GeneScoreMatrix",groupBy="Clusters",bias=c("TSSEnrichment","log10(nFrags)"),testMethod="wilcoxon")
marker_gene_list <- getMarkers(seMarker=marker_genes_results,cutOff="FDR <= 0.05 & Log2FC >= 1"); print(marker_gene_list)

# 12. Gene Activity Visualization with MAGIC
top_marker_genes <- unique(unlist(lapply(marker_gene_list,function(x) head(x$name,5))))
top_marker_genes <- head(top_marker_genes,5); print(top_marker_genes)

umap_gene_activity <- plotEmbedding(ArchRProj=proj_peaks,colorBy="GeneScoreMatrix",name=top_marker_genes,embedding="UMAP_Harmony")
plotPDF(plotList=umap_gene_activity,name="Marker_Gene_Activity_UMAP",ArchRProj=proj_peaks,addDOC=FALSE)

proj_peaks <- addImputeWeights(ArchRProj=proj_peaks,reducedDims="Harmony",dimsToUse=downstream_lsi_dims,seed=42)
impute_weights <- getImputeWeights(proj_peaks)

umap_gene_activity_magic <- plotEmbedding(ArchRProj=proj_peaks,colorBy="GeneScoreMatrix",name=top_marker_genes,embedding="UMAP_Harmony",imputeWeights=impute_weights)
plotPDF(plotList=umap_gene_activity_magic,name="Marker_Gene_Activity_UMAP_MAGIC",ArchRProj=proj_peaks,addDOC=FALSE)

# 13. Transcription-Factor Motif Activity
proj_peaks <- addMotifAnnotations(ArchRProj=proj_peaks,motifSet="cisbp",name="Motif",force=TRUE)
print(getPeakAnnotation(proj_peaks,"Motif"))
proj_peaks <- addBgdPeaks(ArchRProj=proj_peaks,force=TRUE)
proj_peaks <- addDeviationsMatrix(ArchRProj=proj_peaks,peakAnnotation="Motif",force=TRUE)

variable_motifs <- getVarDeviations(ArchRProj=proj_peaks,name="MotifMatrix",plot=FALSE)
top_variable_motifs <- head(variable_motifs,2); top_variable_motif_names <- top_variable_motifs$name; print(top_variable_motif_names)

motif_features <- getFeatures(proj_peaks,select=paste(top_variable_motif_names,collapse="|"),useMatrix="MotifMatrix")
top_variable_z_motifs <- grep("^z:",motif_features,value=TRUE); print(top_variable_z_motifs)

motif_umap <- plotEmbedding(ArchRProj=proj_peaks,colorBy="MotifMatrix",name=top_variable_z_motifs,embedding="UMAP_Harmony")
plotPDF(plotList=motif_umap,name="Top_Variable_TF_Motif_Activity",ArchRProj=proj_peaks,addDOC=FALSE)

motif_activity_plots <- lapply(top_variable_z_motifs,function(motif) plotGroups(ArchRProj=proj_peaks,groupBy="Clusters",colorBy="MotifMatrix",name=motif,plotAs="violin"))
plotPDF(plotList=motif_activity_plots,name="Top_Variable_Motif_Activity_by_Cluster",ArchRProj=proj_peaks,addDOC=FALSE)

# 14. Integration with scRNA-seq
seurat_data <- readRDS("data/blish_awilk_seu_subset.rds")
print(seurat_data); print(colnames(seurat_data@meta.data))

proj_peaks <- addGeneIntegrationMatrix(ArchRProj=proj_peaks,useMatrix="GeneScoreMatrix",matrixName="GeneIntegrationMatrix",reducedDims="Harmony",dimsToUse=downstream_lsi_dims,seRNA=seurat_data,groupRNA="cell.type",nameCell="predictedCell",nameGroup="predictedGroup",nameScore="predictedScore",addToArrow=TRUE,force=TRUE)

integration_umap <- plotEmbedding(ArchRProj=proj_peaks,colorBy="GeneIntegrationMatrix",name=integration_genes,embedding="UMAP_Harmony")
plotPDF(plotList=integration_umap,name="RNA_Expression_on_ATAC_UMAP",ArchRProj=proj_peaks,addDOC=FALSE)

# 15. Gene Activity versus Gene Expression
cor_results <- correlateMatrices(ArchRProj=proj_peaks,useMatrix1="GeneIntegrationMatrix",useMatrix2="GeneScoreMatrix",reducedDims="Harmony",dimsToUse=downstream_lsi_dims)
print(head(cor_results))

highest_agreement <- cor_results[which.max(cor_results$cor),]
lowest_agreement <- cor_results[which.min(cor_results$cor),]
print("Highest agreement:"); print(highest_agreement)
print("Lowest agreement:"); print(lowest_agreement)

write.csv(as.data.frame(cor_results),file=file.path(results_dir,"ATAC_RNA_gene_correlations.csv"),row.names=FALSE)

# 16. RNA-Derived Cell-Type Labels
rna_label_umap <- plotEmbedding(ArchRProj=proj_peaks,colorBy="cellColData",name="predictedGroup",embedding="UMAP_Harmony")
plotPDF(rna_label_umap,name="RNA_Inferred_Cell_Types",ArchRProj=proj_peaks,addDOC=FALSE)

confusion_matrix <- table(ATAC_Clusters=proj_peaks$Clusters,RNA_Labels=proj_peaks$predictedGroup)
print(confusion_matrix)
write.csv(as.data.frame(confusion_matrix),file=file.path(results_dir,"ATAC_RNA_confusion_matrix.csv"),row.names=FALSE)

cluster_label_map <- apply(confusion_matrix,1,function(x) { if(sum(x)==0) return(NA_character_); names(x)[which.max(x)] })
proj_peaks$CellType <- unname(cluster_label_map[as.character(proj_peaks$Clusters)])
print(table(proj_peaks$Clusters,proj_peaks$CellType))

# 17. Peak-to-Gene Linkage
proj_peaks <- addPeak2GeneLinks(ArchRProj=proj_peaks,reducedDims="Harmony",dimsToUse=downstream_lsi_dims,useMatrix="GeneIntegrationMatrix",corCutOff=0.75,force=TRUE)
peak_gene_links <- getPeak2GeneLinks(ArchRProj=proj_peaks,corCutOff=0.5,resolution=10000); print(peak_gene_links)

peak_gene_heatmap <- plotPeak2GeneHeatmap(ArchRProj=proj_peaks)
plotPDF(peak_gene_heatmap,name="Peak_to_Gene_Linkage",ArchRProj=proj_peaks,addDOC=FALSE,width=8,height=8)

# 18. COVID versus Healthy Differential Accessibility
condition_map <- c(
  "ATAC_555_2_fragments_fragments"="COVID",
  "ATAC_557_fragments_fragments"="COVID",
  "ATAC_EV08_fragments_fragments"="Healthy",
  "ATAC_HIP02_frozen_fragments_fragments"="Healthy"
)

proj_peaks$Condition <- unname(condition_map[as.character(proj_peaks$Sample)])
if(any(is.na(proj_peaks$Condition))) stop("Some samples were not assigned to COVID/Healthy. Check condition_map.")
print(table(proj_peaks$Condition)); print(table(proj_peaks$Condition,proj_peaks$Sample))

differential_peaks <- getMarkerFeatures(ArchRProj=proj_peaks,useMatrix="PeakMatrix",groupBy="Condition",bias=c("TSSEnrichment","log10(nFrags)"),testMethod="wilcoxon")

covid_peaks <- getMarkers(differential_peaks,cutOff="FDR <= 0.1 & Log2FC >= 0.5")$COVID
healthy_peaks <- getMarkers(differential_peaks,cutOff="FDR <= 0.1 & Log2FC >= 0.5")$Healthy
print(covid_peaks); print(healthy_peaks)

ma_plot <- plotMarkers(seMarker=differential_peaks,name="COVID",cutOff="FDR <= 0.05 & abs(Log2FC) >= 0.5",plotAs="MA")
volcano_plot <- plotMarkers(seMarker=differential_peaks,name="COVID",cutOff="FDR <= 0.05 & abs(Log2FC) >= 0.5",plotAs="Volcano")
plotPDF(plotList=list(ma_plot,volcano_plot),name="COVID_vs_Healthy_Differential_Accessibility",ArchRProj=proj_peaks,addDOC=FALSE)

# 19. TF Motif Enrichment
healthy_motif_enrichment <- peakAnnoEnrichment(seMarker=differential_peaks,ArchRProj=proj_peaks,peakAnnotation="Motif",cutOff="FDR <= 0.1 & Log2FC >= 0.5",background="bgdPeaks")
covid_motif_enrichment <- peakAnnoEnrichment(seMarker=differential_peaks,ArchRProj=proj_peaks,peakAnnotation="Motif",cutOff="FDR <= 0.1 & Log2FC >= 0.5",background="bgdPeaks")
print(healthy_motif_enrichment); print(covid_motif_enrichment)

get_top_enriched_motifs <- function(enrichment_result,n=5) {
  if(nrow(enrichment_result)==0) return(character(0))
  enrichment_scores <- assay(enrichment_result,"mlog10Padj")
  if(is.null(dim(enrichment_scores))) enrichment_scores <- as.numeric(enrichment_scores) else enrichment_scores <- enrichment_scores[,1]
  valid <- !is.na(enrichment_scores)
  motif_names <- rownames(enrichment_result)[valid]
  enrichment_scores <- enrichment_scores[valid]
  motif_names[order(enrichment_scores,decreasing=TRUE)][seq_len(min(n,length(motif_names)))]
}

top_healthy_motifs <- get_top_enriched_motifs(healthy_motif_enrichment,n=5)
top_covid_motifs <- get_top_enriched_motifs(covid_motif_enrichment,n=5)
print(top_healthy_motifs); print(top_covid_motifs)

motifs_for_visualisation <- unique(c(top_healthy_motifs,top_covid_motifs))

if(length(motifs_for_visualisation)>0) {
  motif_activity_umap <- plotEmbedding(ArchRProj=proj_peaks,colorBy="MotifMatrix",name=motifs_for_visualisation,embedding="UMAP_Harmony")
  plotPDF(plotList=motif_activity_umap,name="Condition_Associated_Motif_Activity",ArchRProj=proj_peaks,addDOC=FALSE)

  motif_activity_heatmap <- plotGroups(ArchRProj=proj_peaks,groupBy="Clusters",colorBy="MotifMatrix",name=motifs_for_visualisation,plotAs="heatmap")
  plotPDF(plotList=list(motif_activity_heatmap),name="Condition_Associated_Motif_Activity_Heatmap",ArchRProj=proj_peaks,addDOC=FALSE)
} else message("No condition-associated motifs passed the enrichment criteria.")

# 20. TF Footprinting
if(length(motifs_for_visualisation)>0) {
  motif_positions <- getPositions(ArchRProj=proj_peaks,name="Motif")
  footprint_motifs <- unlist(lapply(motifs_for_visualisation,function(motif) grep(motif,names(motif_positions),value=TRUE)))
  footprint_motifs <- unique(footprint_motifs)

  if(length(footprint_motifs)>0) {
    proj_peaks <- addGroupCoverages(ArchRProj=proj_peaks,groupBy="Condition",force=TRUE)
    footprint_positions <- motif_positions[footprint_motifs]
    footprint_results <- getFootprints(ArchRProj=proj_peaks,positions=footprint_positions,groupBy="Condition",useGroups=c("COVID","Healthy"),minCells=25)
    footprint_plot <- plotFootprints(seFoot=footprint_results,ArchRProj=proj_peaks,normMethod="Subtract",smoothWindow=5,addDOC=FALSE,plot=TRUE)
  } else message("No matching motif positions were available for footprinting.")
} else message("Condition-specific footprinting was not performed because no significant condition-associated motifs were identified.")

# 21. Co-accessibility
proj_peaks <- addCoAccessibility(ArchRProj=proj_peaks,reducedDims="IterativeLSI",dimsToUse=downstream_lsi_dims,force=TRUE)
coaccessibility_links <- getCoAccessibility(ArchRProj=proj_peaks,corCutOff=0.5,resolution=10000); print(coaccessibility_links)

coaccessibility_marker_genes <- c("CD14","CD8A")
coaccessibility_tracks <- plotBrowserTrack(ArchRProj=proj_peaks,groupBy="Condition",geneSymbol=coaccessibility_marker_genes,upstream=50000,downstream=50000,loops=coaccessibility_links)
plotPDF(plotList=coaccessibility_tracks,name="CoAccessibility_Marker_Genes",ArchRProj=proj_peaks,addDOC=FALSE,width=8,height=8)

# 22. Potential Regulatory Elements
regulatory_tracks <- plotBrowserTrack(ArchRProj=proj_peaks,groupBy="Condition",geneSymbol=marker_genes,upstream=50000,downstream=50000,loops=coaccessibility_links)
plotPDF(plotList=regulatory_tracks,name="Marker_Gene_Regulatory_Elements",ArchRProj=proj_peaks,addDOC=FALSE,width=8,height=8)

# 23. Final Analysis Summary
final_summary <- data.frame(
  metric=c("Final cells","Final median TSS enrichment","Final median fragments","Number of clusters","Number of peak features","Number of gene features"),
  value=c(nCells(proj_peaks),median(proj_peaks$TSSEnrichment),median(proj_peaks$nFrags),length(unique(proj_peaks$Clusters)),length(getPeakSet(proj_peaks)),nrow(getMatrixFromProject(proj_peaks,useMatrix="GeneScoreMatrix")))
)
print(final_summary)
write.csv(final_summary,file=file.path(results_dir,"final_analysis_summary.csv"),row.names=FALSE)

# 24. Save Final ArchR Project
# Optional: save the complete final ArchR project locally.
# Keep commented for GitHub because the ArchR project can be very large.
# proj_peaks <- saveArchRProject(ArchRProj=proj_peaks,outputDirectory=file.path(results_dir,"ArchRProject_Final"),load=TRUE)

# 25. Session Information
sessionInfo()
