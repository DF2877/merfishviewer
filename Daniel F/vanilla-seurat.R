library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
cea.data <- Read10X(data.dir = "D:/df/nextseq_211004/2/s2/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
cea <- CreateSeuratObject(counts = cea.data, project = "ceah2bgfp", min.cells = 3, min.features = 200)
cea

targetIDs <- read.csv("D:/df/targets_and_families.csv", skip=1)
rownames(targetIDs) <- make.names(targetIDs$MGI.symbol, unique = TRUE)

cea[["RNA"]] <- AddMetaData(cea[["RNA"]], targetIDs)

plot1 <- VlnPlot(cea, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
plot2 <- FeatureScatter(cea, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Stop and subset ##

cea <- subset(cea, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)

cea <- SCTransform(cea, ncells = 7500, n_genes = 4500, variable.features.n = 4500, conserve.memory = TRUE, verbose = TRUE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cea), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cea)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

cea <- RunPCA(cea, features = VariableFeatures(object = cea))

DimPlot(cea, reduction = "pca") + NoLegend()

DimHeatmap(cea, dims = 1:15, cells = 500, balanced = TRUE)

cea <- IntegrateLayers(object = cea, method = CCAIntegration, normalization.method = "SCT")

ElbowPlot(cea)

## Stop and set dims ##
## Can skip clustering if using metadata ##

dimm <- 1:14

cea <- FindNeighbors(cea, dims = dimm)
cea <- FindClusters(cea, resolution = 0.5)

cea <- RunUMAP(cea, dims = dimm)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(cea, reduction = "umap", label = TRUE)

## Stop and set resolution ##

FeaturePlot(cea, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"))
VlnPlot(cea, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"), stack = TRUE, same.y.lims = TRUE)

## Resume here if clustering skipped ##

cea <- PrepSCTFindMarkers(cea)
# cea.markers <- FindAllMarkers(cea)

saveRDS(cea, file = "D:/df/nextseq_211004/basic-processed.rds")
# cea <- readRDS("D:/df/nextseq_211004/basic-processed.rds")

## If using a metadata condition, set ident.1 equal to the ##
## metadata condition, and set group.by equal to the metadata column name ##
markers <- FindMarkers(cea, ident.1 = 15)

gc()
targets <- targetIDs$Type[match(rownames(markers),targetIDs$MGI.symbol)]
merfish <- !is.na(merfishSymbols[match(rownames(markers),merfishSymbols)])
# mgiSymbols <- d@assays[["RNA"]]@meta.data$gene_symbol[match(rownames(markers),d@assays[["RNA"]]@meta.data$gene_identifier)]

markers <- cbind(markers, targets, merfish)
# rownames(markers) <- make.names(mgiSymbols, unique = TRUE)
gc()

interest <- markers[markers$targets %in% c("gpcr","vgic","lgic","other_ic"),]
interest <- interest[interest$p_val_adj < 0.05,]

write.csv(markers, file = "D:/df/nextseq_211004/markers.csv")
write.csv(interest, file = "D:/df/nextseq_211004/gpcr-ic.csv")
