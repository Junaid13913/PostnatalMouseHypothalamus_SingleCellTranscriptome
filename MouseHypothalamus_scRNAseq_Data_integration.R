# Loading necessary libraries
# Load required libraries
# Seurat: For single-cell RNA sequencing (scRNA-seq) data analysis and visualization
# ggplot2: For creating elegant data visualizations
# dplyr: For data manipulation and transformation
library(Seurat)
library(ggplot2)
library(dplyr)

# Merge all datasets into a single Seurat object
# This step consolidates all individual datasets into one combined dataset for downstream analysis.
hyp <- merge(GSE74672_seu, 
             y = c(GSE87544_seu, GSE146692_seu, GSE132355_se,  GSE125065_seu, GSE126836_seu, GSE93374_seu, GSE139923_seu, GSE113576_seu, Anderson_seu), 
             add.cell.ids = c("GSE74672", "GSE87544", "GSE146692", "GSE132355", "GSE125065", "GSE126836", "GSE93374", "GSE139923", "GSE113576", "Anderson"), 
             project = "Hypothalamus", 
             merge.data = TRUE) 

# Add metadata column 'new' to the Seurat object based on the 'Study' column
hyp <- AddMetaData(hyp, hyp$Study, col.name = "new")

# Renaming values in the 'new' column for better understanding and clarity
hyp$new[hyp$new==c("GSE139923")] <- "Connect_seq"
hyp$new[hyp$new==c("Anderson")] <- "Retro_seq"
# Remaining datasets are grouped as 'Rest'
hyp$new[hyp$new==c("GSE87544")] <- "Rest"
hyp$new[hyp$new==c("GSE146692")] <- "Rest"
hyp$new[hyp$new==c("GSE132355")] <- "Rest"
hyp$new[hyp$new==c("GSE125065")] <- "Rest"
hyp$new[hyp$new==c("GSE126836")] <- "Rest"
hyp$new[hyp$new==c("GSE93374")] <- "Rest"
hyp$new[hyp$new==c("GSE113576")] <- "Rest"
hyp$new[hyp$new==c("GSE74672")] <- "Rest"
# Quality control and cell filtering steps
# This section aims to filter out cells that might be of low quality or artifacts.
Idents(object = hyp) <- "hyp"
hyp[["percent.mt"]] <- PercentageFeatureSet(hyp, pattern = "^mt-")
# Visualizing feature distributions across cells
VlnPlot(hyp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0) + NoLegend()
# Filtering cells based on certain thresholds to retain high-quality cells
hyp <- subset(hyp, subset = nFeature_RNA < 2800 & nFeature_RNA > 100 & percent.mt < 15)

# Standard pre-processing steps to visualize batch effects
# The following steps normalize the data, identify variable features, scale the data, and reduce its dimensionality.
hyp <- NormalizeData(object = hyp)
hyp <- FindVariableFeatures(object = hyp)
hyp <- ScaleData(object = hyp)
hyp <- RunPCA(object = hyp)
# Visualization of PCA results to inspect potential batch effects
VizDimLoadings(valid_data_seu, dims = 1:2, reduction = "pca")
DimHeatmap(hyp, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(hyp)
# Further pre-processing steps for clustering and visualization
hyp <- FindNeighbors(object = hyp, dims = 1:20)
hyp <- FindClusters(object = hyp)
hyp <- RunUMAP(object = hyp, dims = 1:20)

# Integration process to correct for batch effects
# This section works on integrating data from different studies to mitigate any batch effects.
hyp.list <-  SplitObject(hyp, split.by = "Study")

# Standard preprocessing on each object in the list
for (i in 1:length(hyp.list)) {
  hyp.list[[i]] <- NormalizeData(hyp.list[[i]], verbose = FALSE)
  hyp.list[[i]] <- subset(hyp.list[[i]])
  hyp.list[[i]] <- FindVariableFeatures(hyp.list[[i]], selection.method = "vst", verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = hyp.list)
hyp.list <- lapply(X = hyp.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Finding and applying integration anchors
anchors <- FindIntegrationAnchors(object.list = hyp.list)
hyp <- IntegrateData(anchorset = anchors)

# Post-integration processing and visualization
# After integration, the data undergoes further preprocessing, dimensionality reduction, and visualization.
hyp <- ScaleData(hyp, verbose = FALSE)
hyp <- FindVariableFeatures(hyp, selection.method = "vst", verbose = FALSE)
hyp <- RunPCA(hyp, verbose = FALSE)
DimPlot(hyp, reduction = "pca", group.by = "Study")
VizDimLoadings(hyp, dims = 1:2, reduction = "pca")
hyp <- RunUMAP(hyp, reduction = "pca", dims = 1:20)
hyp <- FindNeighbors(hyp, reduction = "pca", dims = 1:20)
hyp <- FindClusters(hyp, resolution = c(0.25))

# Visualizing UMAP embeddings grouped by different metadata categories
# This allows for a visual assessment of the data based on various attributes.
# Multiple plots are created to visualize the data grouped by different metadata categories.
DimPlot(hyp, reduction = 'umap', group.by = 'Study', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Sex', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Age', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 10))
DimPlot(hyp, reduction = 'umap', group.by = 'Strain', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Technology', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Region', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 9))
DimPlot(hyp, reduction = 'umap', group.by = 'new', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))

# Save the hyp object to an RDS file
saveRDS(hyp, file= "exmpl15.rds") 

# Load the saved 'exmpl15.rds' file
exmpl15 <- readRDS("~/workbench/Hypothalamus/exmpl15.rds")

# Feature plots for canonical markers in non-neural and neural clusters
FeaturePlot(exmpl15, features = c("Sox9", "Olig1", "C1qa","Mustn1", "Cldn5"), cols = c("lightgrey", "#f40909", "lightgrey"))
FeaturePlot(exmpl15, features = c("Syp", "Elavl2", "Slc17a6", "Slc32a1", "Snap25"), cols =  c("lightgrey", "#203354", "#ff0000"))

# Violin plots for expression of canonical markers in neural and non-neural clusters
#### Cluster2 OPC
VlnPlot(exmpl15, features = c("Pdgfra", "Gpr17"), pt.size = 0)
#### Cluster3 Astro
VlnPlot(exmpl15, features = c("Agt", "Sox2"), pt.size = 0)
###### Cluster Endo1/Endo2
VlnPlot(exmpl15, features = c("Itm2a", "Cldn5", "Pecam1"), pt.size = 0)
###### Cluster Micro
VlnPlot(exmpl15, features = c("Cx3cr1", "Tmem119"), pt.size = 0)
#### Cluster  Tany
VlnPlot(exmpl15, features = c("Rax", "Col23a1", "Slc16a2", "Lhx2"), pt.size = 0)
##### Cluster MO
VlnPlot(exmpl15, features = c("Mobp", "Apod", "Mog", "Plp1", "Mbp", "Pmp22", "Klk6", "Tspan2", "Mog"), pt.size = 0)
### Cluster  NFO
VlnPlot(exmpl15, features = c("Fyn", "Bmp4", "Sirt2", "Gpr17"), pt.size = 0)
####Cluster  Mural 
VlnPlot(exmpl15, features = c("Mustn1", "Vim"), pt.size = 0)
### Cluster Tany-like-ependy
VlnPlot(exmpl15, features = c("Ccdc153", "Six3"), pt.size = 0)
#### Cluster  Macr
VlnPlot(exmpl15, features = c("Mrc1"), pt.size = 0)
#### cluster  Eryth
VlnPlot(exmpl15, features = c("Alas2", "Ube2l6"), pt.size = 0)
#### cluster  IPC
VlnPlot(exmpl15, features = c("Mki67", "Hmgb2", "Ascl1"), pt.size = 0)
####### VLMCS
VlnPlot(exmpl15, features = c("Dcn"), pt.size = 0)
#### Neuronal clusters
VlnPlot(exmpl15, features = c("Slc32a1", "Slc17a6", "Cck", "Elavl2", "Syp", "Gad1", "Gad2"), pt.size = 0)

# Rename cluster identities for better interpretability
new.cluster.ids <- c("GABA1", "GABA2", "MO", "Astro2", "GLUT1", "OPC", "GABA3", "Endo1", "Tany", "Micro", "GLUT2", "Mural", "NFO", "Mac", "Tany-like-Ependy", "Eryth", "Agt+/Olig1+", "VLMCs", "IPCs", "Astro1", "Endo2", "Astro3", "Cck-N")
names(new.cluster.ids) <- levels(exmpl15)
exmpl15 <- RenameIdents(exmpl15, new.cluster.ids)

# Plot UMAP with cluster labels
DimPlot(exmpl15, reduction = "umap", label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

# Save the modified object to an RDS file
saveRDS(exmpl15, file= "exmpl15.rds")


# Subset Analysis
##########

# Subsetting the dataset to include only neurons based on specific cluster identities
neuron <- subset(exmpl15, idents = c("GABA1", "GABA2", "GLUT1", "GLUT2", "Cck-N", "GABA3"))

# Visualizing the UMAP plot for the subset of neurons
DimPlot(neuron, reduction = "umap", label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

# Displaying the dimensions (number of cells and features) of the subset
dim(neuron)

# Scaling the data for all genes
all.genes <- rownames(neuron)
neuron <- ScaleData(neuron, features = all.genes)

# Running PCA on the scaled data using variable features
neuron <- RunPCA(neuron, features = VariableFeatures(object = neuron))

# Displaying the elbow plot to determine the number of principal components
ElbowPlot(neuron)

# Running UMAP dimensional reduction on the first 20 principal components
neuron <- RunUMAP(neuron, dims = 1:20)

# Finding neighbors based on the PCA
neuron <- FindNeighbors(neuron, reduction = "pca", dims = 1:20)

# Identifying clusters within the subset
neuron <- FindClusters(neuron, resolution = 0.3)

# Updating the metadata for cell types
exmpl15$Cell_type <- as.character(Idents(exmpl15))

# Adding a new metadata column for neuron clusters
neuron <- AddMetaData(neuron,
                      metadata=as.character(neuron$seurat_clusters),
                      col.name="neuron_cluster")

# Setting the default assay to "RNA" for subsequent analysis
DefaultAssay(neuron) <- "RNA"

# Renaming neuron cluster identities based on the seurat clusters
neuron$neuron_cluster[neuron$seurat_clusters == 0] <- "GABA1"
neuron$neuron_cluster[neuron$seurat_clusters == 1] <- "GABA2"
neuron$neuron_cluster[neuron$seurat_clusters == 12] <- "GLUT4"
neuron$neuron_cluster[neuron$seurat_clusters == 13] <- "GABA2"
neuron$neuron_cluster[neuron$seurat_clusters == 2] <- "GLUT1"
neuron$neuron_cluster[neuron$seurat_clusters == 3] <- "GABA3"
neuron$neuron_cluster[neuron$seurat_clusters == 4] <- "GABA4"
neuron$neuron_cluster[neuron$seurat_clusters == 5] <- "GABA5"
neuron$neuron_cluster[neuron$seurat_clusters == 6] <- "GLUT2"
neuron$neuron_cluster[neuron$seurat_clusters == 7] <- "GABA6"
neuron$neuron_cluster[neuron$seurat_clusters == 8] <- "GLUT3"
neuron$neuron_cluster[neuron$seurat_clusters == 9] <- "GABA7"
neuron$neuron_cluster[neuron$seurat_clusters == 10] <- "GLUT4"
neuron$neuron_cluster[neuron$seurat_clusters == 11] <- "GABA8"
neuron$neuron_cluster[neuron$seurat_clusters == 14] <- "Cck-N"

# Visualizing the UMAP plot colored by the newly defined neuron clusters
DimPlot(neuron, reduction = "umap", group.by = "neuron_cluster",  label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

# Generating violin plots for neurotransmitter-specific markers and discriminatory markers for neurons

VlnPlot(neuron, features = c("Slc32a1"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Slc17a6"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Gad2"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Gad1"), pt.size = 0, group.by = "neuron_cluster") 
VlnPlot(neuron, features = c("Klhl1"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Tubb3"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Meis2"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Six3"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Htr2c"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Npy2r"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Avp"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Lhx8"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Crabp1"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Agrp"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Pomc"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Adcyap1"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Lhx1"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Tmem163"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Slc1a3"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Cck"), pt.size = 0, group.by = "neuron_cluster")
VlnPlot(neuron, features = c("Nr4a2"), pt.size = 0, group.by = "neuron_cluster")

# Saving the modified neuron object for future use
saveRDS(neuron, file= "neuro_meta.rds")

# Mapping cell barcodes of the neuronal subset back to the original Seurat object based on their cluster assignments

# Update the "Cell_type" column in the exmpl15 object for cells that match those in the "GABA1" cluster of the neuron object
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA1"))), collapse = "|"))] <- "GABA1"

# Update the "Cell_type" column for cells that match those in the "GABA2" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA2"))), collapse = "|"))] <- "GABA2"

# Update the "Cell_type" column for cells that match those in the "GABA3" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA3"))), collapse = "|"))] <- "GABA3"

# Update the "Cell_type" column for cells that match those in the "GABA4" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA4"))), collapse = "|"))] <- "GABA4"

# Update the "Cell_type" column for cells that match those in the "GABA5" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA5"))), collapse = "|"))] <- "GABA5"

# Update the "Cell_type" column for cells that match those in the "GABA6" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA6"))), collapse = "|"))] <- "GABA6"

# Update the "Cell_type" column for cells that match those in the "GABA7" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA7"))), collapse = "|"))] <- "GABA7"

# Update the "Cell_type" column for cells that match those in the "GABA8" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA8"))), collapse = "|"))] <- "GABA8"

# Update the "Cell_type" column for cells that match those in the "GLUT1" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT1"))), collapse = "|"))] <- "GLUT1"

# Update the "Cell_type" column for cells that match those in the "GLUT2" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT2"))), collapse = "|"))] <- "GLUT2"

# Update the "Cell_type" column for cells that match those in the "GLUT3" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT3"))), collapse = "|"))] <- "GLUT3"

# Update the "Cell_type" column for cells that match those in the "GLUT4" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT4"))), collapse = "|"))] <- "GLUT4"

# Update the "Cell_type" column for cells that match those in the "Cck-N" cluster
exmpl15$Cell_type[mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("Cck-N"))), collapse = "|"))] <- "Cck-N"

#SAME FOR TANYCYTES, ASTROCYTES AND MATURE OLIGODENTROCYTES

# Saving the modified final seurat object after cell types and subtypes analyses for future use
saveRDS(exmpl15, file= "final_hypo_ann.rds")

# -----------------------------------
# UMAP Visualization with Final Annotations
# -----------------------------------

# Generate UMAP plot with custom color mapping for each cell type, showing the distribution of various cell types in the dataset.
DimPlot(object = final_hypo_ann, reduction = 'umap', group.by = "Cell_type", label = T, raster = FALSE, cols = c("IPCs" = "#000000","OPC" = "#CC0000", ...)) + theme_bw()

# -----------------------------------
# Expression Analysis for Proliferation Markers
# -----------------------------------

# Assign the Seurat object 'final_hypo_ann' to 'Seu_test' for downstream analysis.
Seu_test <- final_hypo_ann
DefaultAssay(Seu_test) <- "RNA"  # Set the default assay to RNA for subsequent operations.

# Peek into the metadata of the Seu_test object to ensure data integrity.
view(Seu_test@meta.data)

# Annotate cells based on 'Mki67' expression, which is a marker for proliferating cells.
Seu_test <- AddMetaData(Seu_test, metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Mki67"),], col.name="Mki67_STATUS")
Seu_test$Mki67_STATUS[Seu_test$Mki67_STATUS==0] <- "Mki67_negative"
Seu_test$Mki67_STATUS[Seu_test$Mki67_STATUS!=c("Mki67_negative")] <- "Mki67_positive"

# Annotate cells based on 'Sox2' expression, which is a marker for neural progenitors.
Seu_test <- AddMetaData(Seu_test, metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Sox2"),], col.name="Sox2_STATUS")
Seu_test$Sox2_STATUS[Seu_test$Sox2_STATUS==0] <- "Sox2_negative"
Seu_test$Sox2_STATUS[Seu_test$Sox2_STATUS!=c("Sox2_negative")] <- "Sox2_positive"

# Annotate cells based on 'Hmgb2' expression. HMGB2 is involved in various DNA-related processes.
Seu_test <- AddMetaData(Seu_test, metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Hmgb2"),], col.name="Hmgb2_STATUS")
Seu_test$Hmgb2_STATUS[Seu_test$Hmgb2_STATUS==0] <- "Hmgb2_negative"
Seu_test$Hmgb2_STATUS[Seu_test$Hmgb2_STATUS!=c("Hmgb2_negative")] <- "Hmgb2_positive"

# -----------------------------------
# Visualize Proliferation Marker Expression on UMAP
# -----------------------------------

# Plot UMAP visualizations of the three proliferation markers (Mki67, Sox2, and Hmgb2) to understand their distribution across cells.
DimPlot(Seu_test, reduction = "umap", group.by = "Mki67_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw()
DimPlot(Seu_test, reduction = "umap", group.by = "Sox2_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw()
DimPlot(Seu_test, reduction = "umap", group.by = "Hmgb2_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw()


# Load essential libraries for the analysis.
library(Seurat)
library(tidyverse)
library(ggplot2)

# Specify the directory containing validation data.
data_dirv <- "/mnt/S7/data2/workbench/Users/sblim/Hypothalamus/validation data"
# List the files in the specified directory.
list.files(data_dirv)

# Read the 10X Genomics data from the directory.
valid_data <- Read10X(data.dir = data_dirv)

# Create a Seurat object using the count matrix.
valid_data_seu   <- CreateSeuratObject(count = valid_data)

# Add metadata to the Seurat object.
# Set study name for all cells.
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("nuConnect-seq"), 1533), col.name = "Study")
# Set age information for all cells.
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("8 weeks"), 1533), col.name = "Age")
# Set brain region for all cells.
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("Whole"), 1533), col.name = "Region")
# Set sequencing technology for all cells.
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("FACS"), 1533), col.name = "Technology")
# Set mouse strain for all cells.
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("CRH-Cre"), 1533), col.name = "Strain")
# Set sex for all cells.
valid_data_seu <- AddMetaData(valid_data_seu , rep(c("Male"), 1533), col.name = "Sex")

# Normalize the data.
valid_data_seu <- NormalizeData(object = valid_data_seu)
# Identify variable features in the data.
valid_data_seu <- FindVariableFeatures(object = valid_data_seu)
# Scale the data.
valid_data_seu <- ScaleData(object = valid_data_seu)
# Perform principal component analysis.
valid_data_seu <- RunPCA(object = valid_data_seu)

# Display heatmap of the top dimensions.
DimHeatmap(hyp, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(valid_data_seu, dims = 1:10, cells = 500, balanced = TRUE)

# Plot PCA results.
DimPlot(valid_data_seu, reduction = "pca", raster = FALSE)

# Display heatmap of the top dimension.
DimHeatmap(valid_data_seu, dims = 1, cells = 500, balanced = TRUE)
# Plot elbow plot to help choose the number of dimensions.
ElbowPlot(valid_data_seu)

# Identify neighboring cells.
valid_data_seu <- FindNeighbors(object = valid_data_seu, dims = 1:20)
# Perform clustering on the data.
valid_data_seu <- FindClusters(object = valid_data_seu, resolution = 0.25)
# Reduce dimensions using UMAP.
valid_data_seu <- RunUMAP(object = valid_data_seu, dims = 1:20)
# Plot the UMAP results.
DimPlot(valid_data_seu, reduction = 'umap', raster = FALSE, label = F)

# Assign cell types based on cluster results.
valid_data_seu$Cell_type <- as.character(Idents(valid_data_seu))
# Manually assign cell type names based on cluster IDs.
# ...
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 0] <- "Mac"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 1] <- "Neurons"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 2] <- "Mac"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 3] <- "Astro"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 4] <- "MO"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 5] <- "Endo1"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 6] <- "Neurons"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 7] <- "Neurons"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 8] <- "IPCs"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 9] <- "Endo2"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 10] <- "OPC"
valid_data_seu$Cell_type[valid_data_seu$seurat_clusters == 11] <- "Neurons"

# Plot violin plots for various genes across identified cell types.
# ...
VlnPlot(valid_data_seu, features = c("Syp"), pt.size = 0, group.by = "Cell_typ")
VlnPlot(valid_data_seu, features = c("Mog"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Itm2a"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Pecam1"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Agt"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Mki67"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Pdgfra"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("C1qa"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Elavl2"), pt.size = 0, group.by = "Cell_type")

# Display a UMAP plot with cell type labels.
DimPlot(object = valid_data_seu, reduction = "umap", group.by = "Cell_type", label = T, repel = TRUE, raster = F) + theme_bw() + NoLegend()

# Save the processed Seurat object as an RDS file.
saveRDS(valid_data_seu, file= "Validation.rds")

# Label Transfer Analysis

# Load integrated dataset and validation dataset
final_hypo_ann <- readRDS("~/workbench/Hypothalamus/final_hypo_ann.rds")
Validation <- readRDS("~/workbench/Hypothalamus/Validation.rds")

# Set reference as the integrated dataset
reference <- final_hypo_ann

# Find anchors between the reference and query datasets using Seurat's integration methods
anchors <- FindTransferAnchors(
  reference = reference,
  query = Validation,
  normalization.method = "LogNormalize", # Specify normalization method
  reduction = "pcaproject",              # Specify reduction method
  reference.reduction = "pca")

# Set default assay to RNA
DefaultAssay(reference) <- "RNA"

# Find variable features for the reference dataset
reference <- FindVariableFeatures(reference)

# Extract variable features from the validation dataset
varfeatures <- VariableFeatures(object = Validation)

# Set default assay to RNA for validation set
DefaultAssay(Validation) <- "RNA"

# Scale and perform PCA on the validation dataset
Validation <- ScaleData(Validation, features = varfeatures2, verbose = F)
Validation <- RunPCA(Validation, features = varfeatures2, verbose = FALSE)

# Find neighbors for the validation set
Validation <- FindNeighbors(Validation, dims = 1:20)

# Compute UMAP embeddings for the reference dataset
reference <- RunUMAP(reference, reduction = "pca", dims = 1:20, return.model = T)

# Create a new UMAP dimensionality reduction object
reference[["umap.new"]] <- CreateDimReducObject(embeddings = reference[["umap"]]@cell.embeddings, key = "UMAPnew", assay = "RNA")

# Transfer the cell type labels from the reference dataset to the validation dataset using the previously computed anchors
Validation <- MapQuery(
  anchorset = anchors,
  query = Validation,
  reference = reference,
  refdata = list(
    celltype = reference$Cell_type),
  reference.reduction = "pca", 
  reduction.model = "umap")

# Visualize the UMAP plots colored by actual and predicted cell types
DimPlot(Validation, reduction = "ref.umap", group.by = "Cell_type", label = T)
DimPlot(Validation, reduction = "ref.umap", group.by = "predicted.celltype", label = T)

# Evaluate the prediction quality by visualizing prediction scores on UMAP plots for various features
DefaultAssay(Validation) <- "prediction.score.celltype"
FeaturePlot(Validation, features = c("GABA4", "GLUT1", "GLUT2", "Mac", "IPCs", "OPC", "Tany", "VLMCs", "Myl-OL"), reduction = "ref.umap", cols = c("lightgrey", "#ff4433"), ncol = 3) & theme(plot.title = element_text(size = 17))



