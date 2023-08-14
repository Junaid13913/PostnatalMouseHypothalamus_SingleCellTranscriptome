library(Seurat)
library(ggplot2)
library(dplyr)

#Merge all datasets
hyp <- merge(GSE74672_seu, y = c(GSE87544_seu, GSE146692_seu, GSE132355_se,  GSE125065_seu, GSE126836_seu, GSE93374_seu, GSE139923_seu, GSE113576_seu, Anderson_seu), 
             add.cell.ids = c("GSE74672", "GSE87544", "GSE146692", "GSE132355", "GSE125065", "GSE126836", "GSE93374", "GSE139923", "GSE113576", "Anderson"), 
             project = "Hypothalamus", 
             merge.data = TRUE) 

####adding colum in reference to study
hyp <- AddMetaData(hyp, hyp$Study, col.name = "new")
#####renaming values
hyp$new[hyp$new==c("GSE139923")] <- "Connect_seq"
hyp$new[hyp$new==c("Anderson")] <- "Retro_seq"
hyp$new[hyp$new==c("GSE87544")] <- "Rest"
hyp$new[hyp$new==c("GSE146692")] <- "Rest"
hyp$new[hyp$new==c("GSE132355")] <- "Rest"
hyp$new[hyp$new==c("GSE125065")] <- "Rest"
hyp$new[hyp$new==c("GSE126836")] <- "Rest"
hyp$new[hyp$new==c("GSE93374")] <- "Rest"
hyp$new[hyp$new==c("GSE113576")] <- "Rest"
hyp$new[hyp$new==c("GSE74672")] <- "Rest"

# Qc and filtering
Idents(object = hyp) <- "hyp"

hyp[["percent.mt"]] <- PercentageFeatureSet(hyp, pattern = "^mt-")
VlnPlot(hyp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0) + NoLegend()
hyp <- subset(hyp, subset = nFeature_RNA < 2800 & nFeature_RNA > 100 & percent.mt < 15)


# perform standard workflow steps to figure out if we see any batch effects --------
hyp <- NormalizeData(object = hyp)
hyp <- FindVariableFeatures(object = hyp)
hyp <- ScaleData(object = hyp)
hyp <- RunPCA(object = hyp)
VizDimLoadings(valid_data_seu, dims = 1:2, reduction = "pca")
DimHeatmap(hyp, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(hyp)
hyp <- FindNeighbors(object = hyp, dims = 1:20)
hyp <- FindClusters(object = hyp)
hyp <- RunUMAP(object = hyp, dims = 1:20)


# perform integration to correct for batch effects ------
hyp.list <-  SplitObject(hyp, split.by = "Study")

# perform standard pre-processing on each object
for (i in 1:length(hyp.list)) {
  hyp.list[[i]] <- NormalizeData(hyp.list[[i]], verbose = FALSE)
  hyp.list[[i]] <- subset(hyp.list[[i]])
  hyp.list[[i]] <- FindVariableFeatures(
    hyp.list[[i]], selection.method = "vst",
    verbose = FALSE
  )
}

features <- SelectIntegrationFeatures(object.list = hyp.list)
hyp.list <- lapply(X = hyp.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# find anchors
anchors <- FindIntegrationAnchors(object.list = hyp.list)

# integrate data
hyp <- IntegrateData(anchorset = anchors)

# Run the standard workflow for visualization and clustering
hyp <- ScaleData(hyp, verbose = FALSE)
hyp <- FindVariableFeatures(hyp, 
                            selection.method = "vst",
                            verbose = FALSE)
...........................
hyp <- RunPCA(hyp, verbose = FALSE)

DimPlot(hyp, reduction = "pca", group.by = "Study")
VizDimLoadings(hyp, dims = 1:2, reduction = "pca")
hyp <- RunUMAP(hyp, reduction = "pca", dims = 1:20)
hyp <- FindNeighbors(hyp, reduction = "pca", dims = 1:20)
hyp <- FindClusters(hyp, resolution = c(0.25))

DimPlot(hyp, reduction = 'umap', group.by = 'Study', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Sex', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Age', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 10))
DimPlot(hyp, reduction = 'umap', group.by = 'Strain', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Technology', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(hyp, reduction = 'umap', group.by = 'Region', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 9))
DimPlot(hyp, reduction = 'umap', group.by = 'new', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))


saveRDS(hyp, file= "exmpl15.rds")

head(final_hypo_ann)
DimPlot(final_hypo_ann, reduction = 'umap', group.by = 'Cell_type', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))

##### Load "exmpl15.rds"
exmpl15 <- readRDS("~/workbench/Hypothalamus/exmpl15.rds")

## canonical MARKERs for non-neural cluster
FeaturePlot(exmpl15, features = c("Sox9", "Olig1", "C1qa","Mustn1", "Cldn5"), cols = c("lightgrey", "#f40909", "lightgrey"))

## canonical MARKERs for neural cluster
FeaturePlot(exmpl15, features = c("Syp", "Elavl2", "Slc17a6", "Slc32a1", "Snap25"), cols =  c("lightgrey", "#203354", "#ff0000"))

### canonical markers for both neural and non-neural clusters
VlnPlot(Object, features = c("Oxt", "Elavl2", "Syp", "Snap25", "Cck", "Dcn", "Slc32a1", "Slc17a6", "Sox9", "Olig1", "C1qa", "Cldn5","Mustn1", "Ube2l6"), pt.size = 0)
VlnPlot(exmpl15, features = c("Syp"), pt.size = 0)
VlnPlot(exmpl15, features = c("Elavl2"), pt.size = 0)
VlnPlot(exmpl15, features = c("Slc17a6"), pt.size = 0)
VlnPlot(exmpl15, features = c("Gad2"), pt.size = 0)
VlnPlot(exmpl15, features = c("Olig1"), pt.size = 0)
VlnPlot(exmpl15, features = c("Sox9"), pt.size = 0)
VlnPlot(exmpl15, features = c("C1qa"), pt.size = 0)
VlnPlot(exmpl15, features = c("Cldn5"), pt.size = 0)
VlnPlot(exmpl15, features = c("Mustn1"), pt.size = 0)
VlnPlot(exmpl15, features = c("Cck"), pt.size = 0)
VlnPlot(exmpl15, features = c("Alas2"), pt.size = 0)


###### Cluster specific markers

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


new.cluster.ids <- c("GABA1", "GABA2", "MO", "Astro2", "GLUT1", "OPC", "GABA3", "Endo1", "Tany", "Micro", "GLUT2", "Mural", "NFO", "Mac", "Tany-like-Ependy", "Eryth", "Agt+/Olig1+", "VLMCs", "IPCs", "Astro1", "Endo2", "Astro3", "Cck-N")
names(new.cluster.ids) <- levels(exmpl15)
exmpl15 <- RenameIdents(exmpl15, new.cluster.ids)

DimPlot(exmpl15, reduction = "umap", label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()
saveRDS(exmpl15, file= "exmpl15.rds")


########subset analysis
####subset neurons
neuron <- subset(exmpl15, idents = c("GABA1", "GABA2", "GLUT1", "GLUT2", "Cck-N", "GABA3"))
DimPlot(neuron, reduction = "umap", label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

dim(neuron)

# perform standard workflow steps
all.genes <- rownames(neuron)
neuron <- ScaleData(neuron, features = all.genes)
neuron <- RunPCA(neuron, features = VariableFeatures(object = neuron))
ElbowPlot(neuron)
neuron <- RunUMAP(neuron, dims = 1:20)
neuron <- FindNeighbors(neuron, reduction = "pca", dims = 1:20)
neuron <- FindClusters(neuron, resolution = 0.3)


exmpl15$Cell_type <- as.character(Idents(exmpl15))
neuron <- AddMetaData(neuron,
                      metadata=as.character(neuron$seurat_clusters),
                      col.name="neuron_cluster")
DefaultAssay(neuron) <- "RNA"
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

DimPlot(neuron, reduction = "umap", group.by = "neuron_cluster",  label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

##### Neurotransmitters specific markers along with discriminatory markers
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

saveRDS(neuron, file= "neuro_meta.rds")

##### mapping cell barcodes of neuronal subset to original seurat object
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA1"))), collapse = "|"))] <- "GABA1"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA2"))), collapse = "|"))] <- "GABA2"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA3"))), collapse = "|"))] <- "GABA3"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA4"))), collapse = "|"))] <- "GABA4"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA5"))), collapse = "|"))] <- "GABA5"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA6"))), collapse = "|"))] <- "GABA6"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA7"))), collapse = "|"))] <- "GABA7"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GABA8"))), collapse = "|"))] <- "GABA8"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT1"))), collapse = "|"))] <- "GLUT1"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT2"))), collapse = "|"))] <- "GLUT2"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT3"))), collapse = "|"))] <- "GLUT3"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("GLUT4"))), collapse = "|"))] <- "GLUT4"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(neuron, subset = neuron_cluster == c("Cck-N"))), collapse = "|"))] <- "Cck-N"

#####subset Tany
tany <- subset(exmpl15, idents = c("Tany"))

all.genes <- rownames(tany)
tany <- ScaleData(tany, features = all.genes)
tany <- RunPCA(tany, features = VariableFeatures(object = tany))
ElbowPlot(tany)
tany<- RunUMAP(tany, dims = 1:20)
tany <- FindNeighbors(tany, reduction = "pca", dims = 1:20)
tany <- FindClusters(tany, resolution = 0.2)

###
tany <- AddMetaData(tany,
                    metadata=as.character(tany$seurat_clusters),
                    col.name="Alphabeta")

DefaultAssay(tany) <- "RNA"

tany$Alphabeta[tany$seurat_clusters == 0] <- "α2"
tany$Alphabeta[tany$seurat_clusters == 1] <- "β2"
tany$Alphabeta[tany$seurat_clusters == 2] <- "α2"
tany$Alphabeta[tany$seurat_clusters == 3] <- "β1"
tany$Alphabeta[tany$seurat_clusters == 4] <- "α1"
tany$Alphabeta[tany$seurat_clusters == 5] <- "Mig-astro"
tany$Alphabeta[tany$seurat_clusters == 6] <- "Mig-astro"
tany$Alphabeta[tany$seurat_clusters == 7] <- "α2"
tany$Alphabeta[tany$seurat_clusters == 8] <- "Mig-astro"

DimPlot(tany, reduction = "umap", group.by = "Alphabeta",  label = TRUE, label.size = 4.5, pt.size = 0.6) + NoLegend()

### plotting
VlnPlot(tany, features = c("S100b"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Scn7a"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Pdzph1"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Penk"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Mef2c"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Col23a1"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Slc16a2"), group.by = "Alphabeta", pt.size = 0)
VlnPlot(tany, features = c("Trpm3"), group.by = "Alphabeta", pt.size = 0)

saveRDS(tany, file= "tany_meta.rds")

##### mapping cell barcodes of neuronal subset to original seurat object

exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(tany, subset = Alphabeta == c("α1"))), collapse = "|"))] <- "α1"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(tany, subset = Alphabeta == c("α2"))), collapse = "|"))] <- "α2"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(tany, subset = Alphabeta == c("β1"))), collapse = "|"))] <- "β1"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(tany, subset = Alphabeta == c("β2"))), collapse = "|"))] <- "β2"
exmpl15$Cell_type [mapply(grepl, colnames(exmpl15), paste(colnames(subset(tany, subset = Alphabeta == c("Mig-astro"))), collapse = "|"))] <- "Mig-astro"

#######subset MO cluster
mo_cluster <- subset(exmpl15, idents = c("MO"))
all.genes <- rownames(mo_cluster)
mo_cluster <- ScaleData(mo_cluster, features = all.genes)
mo_cluster <- RunPCA(mo_cluster, features = VariableFeatures(object = mo_cluster))
ElbowPlot(mo_cluster)
mo_cluster<- RunUMAP(mo_cluster, dims = 1:20)
mo_cluster <- FindNeighbors(mo_cluster, reduction = "pca", dims = 1:20)
mo_cluster <- FindClusters(mo_cluster, resolution = 0.2)

DimPlot(mo_cluster, reduction = 'umap', raster = FALSE) + theme_bw() + theme(legend.text = element_text(size = 14))

######
mo_cluster <- AddMetaData(mo_cluster,
                          metadata=as.character(mo_cluster$seurat_clusters),
                          col.name="Oligo")
DefaultAssay(mo_cluster) <- "RNA"

mo_cluster$Oligo[mo_cluster$seurat_clusters == 0] <- "Myl-OL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 1] <- "MOL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 2] <- "Myl-OL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 3] <- "Myl-OL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 4] <- "Myl-OL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 5] <- "Myl-OL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 6] <- "EMOL"
mo_cluster$Oligo[mo_cluster$seurat_clusters == 7] <- "Myl-OL"

saveRDS(mo_cluster, file= "MO_meta.rds")

### plotting
VlnPlot(mo_cluster, features = c("Klk6"), pt.size = 0, group.by = "Oligo")
VlnPlot(mo_cluster, features = c("Anxa5"), pt.size = 0,  group.by = "Oligo")
VlnPlot(mo_cluster, features = c("Tspan2"), pt.size = 0, group.by = "Oligo")
VlnPlot(mo_cluster, features = c("Plp1"), pt.size = 0, group.by = "Oligo")
VlnPlot(mo_cluster, features = c("Grm3"), pt.size = 0, group.by = "Oligo")
VlnPlot(mo_cluster, features = c("Egr2"), pt.size = 0, group.by = "Oligo")

##### mapping cell barcodes of neuronal subset to original seurat object

data$Cell_type [mapply(grepl, colnames(data), paste(colnames(subset(mo_cluster, subset = Oligo == c("MOL"))), collapse = "|"))] <- "MOL"
data$Cell_type [mapply(grepl, colnames(data), paste(colnames(subset(mo_cluster, subset = Oligo == c("Myl-OL"))), collapse = "|"))] <- "Myl-OL"
data$Cell_type [mapply(grepl, colnames(data), paste(colnames(subset(mo_cluster, subset = Oligo == c("EMOL"))), collapse = "|"))] <- "EMOL"

saveRDS(exmpl15, file= "final_hypo_ann.rds")

##### Astrocytes subtypes specific markers
Astro <- subset(exmpl15, idents = c( "Astro1","Astro2", "Astro3"))
DimPlot(Astro, reduction = "umap", group.by = "Region")
VlnPlot(Astro, features = c("Agt"), pt.size = 0)
VlnPlot(Astro, features = c("Unc13c"), pt.size = 0)
VlnPlot(Astro, features = c("Gfap"), pt.size = 0)
VlnPlot(Astro, features = c("Fam107a"), pt.size = 0)
VlnPlot(Astro, features = c("C1qa"), pt.size = 0)
VlnPlot(Astro, features = c("Rgcc"), pt.size = 0)
VlnPlot(Astro, features = c("Olig2"), pt.size = 0)

### UMAP plot showing final annotation
DimPlot(object = final_hypo_ann, reduction = 'umap', group.by = "Cell_type", label = T, raster = FALSE, cols = c("IPCs" = "#000000","OPC" = "#CC0000", "NFO" = "#3D85C6", "MOL" = "#28b348", "EMOL" = "#CDA34F", "Myl-OL" = "#1C4B6C" , "Astro1" = "#d3ee7d", "Astro2" = "#e4a921", "Astro3" = "#fc7a08", "VLMCs" = "#9212FF",
                                                                                                                 "Tany-like-Ependy" = "#6aa84f", "Agt+/Olig1+"= "#B500D9",
                                                                                                                 "Mig-astro" = "#e7574e", "α1" = '#416649', "α2" = '#AC9179',
                                                                                                                 "β1" = '#6E4F70', "β2" = '#5B85A5',
                                                                                                                 "Mural" = "#595F37","Mac" ="#6B8E4E", "Endo1" = "#7F0000","Endo2"= "#965D79","Eryth"= "#B097CD",
                                                                                                                 "Micro"= "#7752D4" , "GABA1" = "#28b348", "GABA2" = "#CDA34F", "GABA3" = "#1C4B6C", "GABA4" = '#6E4F70', "GABA5" = '#5B85A5', "GABA6"= "#016DFF", "GABA7" = "#B8336A", "GABA8"= "#D39E30", "GLUT1"= "#9F4D36", "GLUT2"="#00204C","GLUT3" = "#AC9179", "GLUT4" = '#416649', "GLUT3" = '#AC9179',  "Cck-N"= "#FF4400")) + theme_bw()

##### expression of proliferating specific markers 
Seu_test <- final_hypo_ann
DefaultAssay(Seu_test) <- "RNA"
view(Seu_test@meta.data)

Seu_test <- AddMetaData(Seu_test,
                        metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Mki67"),],
                        col.name="Mki67_STATUS")

Seu_test$Mki67_STATUS[Seu_test$Mki67_STATUS==0] <- "Mki67_negative"
Seu_test$Mki67_STATUS[Seu_test$Mki67_STATUS!=c("Mki67_negative")] <- "Mki67_positive"

Seu_test <- AddMetaData(Seu_test,
                        metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Sox2"),],
                        col.name="Sox2_STATUS")

Seu_test$Sox2_STATUS[Seu_test$Sox2_STATUS==0] <- "Sox2_negative"
Seu_test$Sox2_STATUS[Seu_test$Sox2_STATUS!=c("Sox2_negative")] <- "Sox2_positive"


Seu_test <- AddMetaData(Seu_test,
                        metadata=as.matrix(GetAssayData(object = Seu_test, slot = "counts"))[c("Hmgb2"),],
                        col.name="Hmgb2_STATUS")

Seu_test$Hmgb2_STATUS[Seu_test$Hmgb2_STATUS==0] <- "Hmgb2_negative"
Seu_test$Hmgb2_STATUS[Seu_test$Hmgb2_STATUS!=c("Hmgb2_negative")] <- "Hmgb2_positive"

##### plotting
DimPlot(Seu_test, reduction = "umap", group.by = "Mki67_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw() #+ NoLegend()
DimPlot(Seu_test, reduction = "umap", group.by = "Sox2_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw() #+ NoLegend()
DimPlot(Seu_test, reduction = "umap", group.by = "Hmgb2_STATUS", label = F, label.size = 3, repel = TRUE, raster = F) + theme_bw() #+ NoLegend()


##### Trajectory analysis 
GLs <- subset(exmpl15, idents = c("Tany-like-Ependy", "Tany", "Astro3", "Astro2", "Astro1"))

library(monocle3)
library(SeuratWrappers)
library(tidyverse)

# ...1 Convert to cell_data_set object 
cds <- as.cell_data_set(GLs)
colData(cds)
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)
reacreate.partitions <- c(rep(1, length(cds@colData@rownames)))

names(reacreate.partitions) <- cds@colData@rownames

reacreate.partitions <- as.factor(reacreate.partitions)
cds@clusters$UMAP$partitions <- reacreate.partitions

list_cluster <- GLs@active.ident
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- GLs@reductions$umap@cell.embeddings

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = "cluster",
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right") + scale_color_manual(values = c ("red", "blue", "green"))
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds,
           color_cells_by = "ident",
           label_cell_groups = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == "Tany-like-Ependy"]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.psedu <- as.data.frame(colData(cds))
ggplot(data.psedu, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) + geom_boxplot()

########slingshot using UMAP reduction

sce <- as.SingleCellExperiment(GLs)
sce  <- slingshot(sce,
                  reducedDim = "UMAP",
                  clusterLabels = GLs$Cell_type,
                  start.clus = "Tany-like-Ependy")
sds <- SlingshotDataSet(sce)

# Add into Seurat object
GLs$slingshot_1 <- sce$slingPseudotime_1
FeaturePlot(GLs, c("slingshot_1"), label = TRUE, 
            cols = rainbow(12))

rm(validat)
######## Loads validation data

library(Seurat)
library(tidyverse)
library(ggplot2)

data_dirv <- "/mnt/S7/data2/workbench/Users/sblim/Hypothalamus/validation data"
list.files(data_dirv)

valid_data <- Read10X(data.dir = data_dirv)

valid_data_seu   <- CreateSeuratObject(count = valid_data)

####metadata standardization
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("Valid_data"), 1533), col.name = "Study")
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("6~14 weeks"), 1533), col.name = "Age")
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("Whole"), 1533), col.name = "Region")
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("FACS"), 1533), col.name = "Technology")
valid_data_seu <- AddMetaData(valid_data_seu, rep(c("CRH-Cre"), 1533), col.name = "Strain")
valid_data_seu <- AddMetaData(valid_data_seu , rep(c("Both"), 1533), col.name = "Sex")
# perform standard workflow steps--------
valid_data_seu <- NormalizeData(object = valid_data_seu)
valid_data_seu <- FindVariableFeatures(object = valid_data_seu)
valid_data_seu <- ScaleData(object = valid_data_seu)
valid_data_seu <- RunPCA(object = valid_data_seu)
DimHeatmap(hyp, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(valid_data_seu, dims = 1:10, cells = 500, balanced = TRUE)
DimPlot(valid_data_seu, reduction = "pca", raster = FALSE)
DimHeatmap(valid_data_seu, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(valid_data_seu)
valid_data_seu <- FindNeighbors(object = valid_data_seu, dims = 1:20)
valid_data_seu <- FindClusters(object = valid_data_seu, resolution = 0.25)
valid_data_seu <- RunUMAP(object = valid_data_seu, dims = 1:20)
DimPlot(valid_data_seu, reduction = 'umap', raster = FALSE, label = F)

valid_data_seu$Cell_type <- as.character(Idents(valid_data_seu))

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

####cluster specific markers
VlnPlot(valid_data_seu, features = c("Syp"), pt.size = 0, group.by = "Cell_typ")
VlnPlot(valid_data_seu, features = c("Mog"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Itm2a"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Pecam1"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Agt"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Mki67"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Pdgfra"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("C1qa"), pt.size = 0, group.by = "Cell_type")
VlnPlot(valid_data_seu, features = c("Elavl2"), pt.size = 0, group.by = "Cell_type")

DimPlot(object = valid_data_seu, reduction = "umap", group.by = "Cell_type", label = T, repel = TRUE, raster = F) +theme_bw() + NoLegend()

saveRDS(valid_data_seu, file= "Validation.rds")

###### Label transfer analysis
#####Load integrated dataset and validation dataset
final_hypo_ann <- readRDS("~/workbench/Hypothalamus/final_hypo_ann.rds")
Validation <- readRDS("~/workbench/Hypothalamus/Validation.rds")
reference <- final_hypo_ann

####### finding transfer anchor between reference and query
anchors <- FindTransferAnchors(
  reference = reference,
  query = Validation,
  normalization.method = "LogNormalize",
  reduction = "pcaproject",
  reference.reduction = "pca")

DefaultAssay(reference) <- "RNA"
reference <- FindVariableFeatures(reference)
varfeatures <- VariableFeatures(object = Validation)
DefaultAssay(Validation) <- "RNA"
Validation <- ScaleData(Validation, features = varfeatures2, verbose = F)
Validation <- RunPCA(Validation, features = varfeatures2, verbose = FALSE)
Validation <- FindNeighbors(Validation, dims = 1:30)

#set umap cell embeddings
reference <- RunUMAP(reference, reduction = "pca", dims = 1:30, return.model = T)
reference[["umap.new"]] <- CreateDimReducObject(embeddings = reference[["umap"]]@cell.embeddings, key = "UMAPnew", assay = "RNA")
#transferring cell type labels from reference to the query
Validation <- MapQuery(
  anchorset = anchors,
  query = Validation,
  reference = reference,
  refdata = list(
    celltype = reference$Cell_type),
  reference.reduction = "pca", 
  reduction.model = "umap")

DimPlot(Validation, reduction = "ref.umap", group.by = "Cell_type", label = T)
DimPlot(Validation, reduction = "ref.umap", group.by = "predicted.celltype", label = T)

#evaluate prediction of cell types
DefaultAssay(Validation) <- "prediction.score.celltype"
FeaturePlot(Validation, features = c("GABA4",
                                     "GLUT1", "GLUT2", "Mac", 
                                     "IPCs", "OPC", "Tany", "VLMCs", "Myl-OL"),
            reduction = "ref.umap", cols = c("lightgrey", "#ff4433"), ncol = 3) & theme(plot.title = element_text(size = 17))
#merge reference and query
reference$id <- 'reference'
Validation$id <- 'query'

table(Validation@meta.data$predicted.celltype.score > 0.50)
DefaultAssay(Validation) <- "RNA"
querytest <- subset(Validation, subset = predicted.celltype.score > 0.50) 

refquerynew <- merge(reference, querytest)
refquerynew[["pca"]] <- merge(reference[["pca"]], querytest[["ref.pca"]])
refquerynew[["umap.new"]] <- merge(reference[["umap"]], querytest[["ref.umap"]])

##### plotting 
DimPlot(refquerynew, group.by = "id", reduction = "umap.new",,  shuffle = TRUE, raster=FALSE, label = F) + theme_bw() 
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Study', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) #+ NoLegend()
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Sex', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) + NoLegend()
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Age', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) 
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Strain', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) + NoLegend()
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Technology', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Region', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) + NoLegend()
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'Cell_type', raster = FALSE, shuffle = TRUE, label = T) + theme_bw() + theme(legend.text = element_text(size = 14))
DimPlot(refquerynew, reduction = 'umap.new', group.by = 'NN_cluster', raster = FALSE, shuffle = TRUE, label = F) + theme_bw() + theme(legend.text = element_text(size = 14)) #+ NoLegend()

saveRDS(refquerynew, file= "refquerynew.rds")
