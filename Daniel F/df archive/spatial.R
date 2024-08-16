library(Seurat)
library(schard)

#Load data
sce <- schard::h5ad2sce("./abc_download_root/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw.h5ad")
sce

assayNames(sce)[1] <- "counts"

cell_metadata <- read.csv("./abc_download_root/metadata/MERFISH-C57BL6J-638850/20231215/cell_metadata.csv")
# gene <- read.csv("./abc_download_root/metadata/MERFISH-C57BL6J-638850/20231215/gene.csv")
rownames(cell_metadata) <- cell_metadata$cell_label
sce@colData$cluster_alias <- cell_metadata$cluster_alias[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$average_correlation_score <- cell_metadata$average_correlation_score[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$feature_matrix_label <- cell_metadata$feature_matrix_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_label <- cell_metadata$donor_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_genotype <- cell_metadata$donor_genotype[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_sex <- cell_metadata$donor_sex[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$x <- cell_metadata$x[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$y <- cell_metadata$y[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$z <- cell_metadata$z[match(sce@colData$cell_label,cell_metadata$cell_label)]
#Setup the Seurat Object
d <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
#d <- SetAssayData(object = d, layer = "data", new.data = logcounts(sce))
# AddMetaData(d,cell_metadata)
# Set feature metadata, AKA rowData. Super intuitive, right?
d[["RNA"]][[]] <- as.data.frame(rowData(sce))
d

#QC and selecting cells for further analysis
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(d[,sample(colnames(d),ncol(d)/20)], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# plot1 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(d[,sample(colnames(d),ncol(d)/20)], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

d <- subset(d, subset = nFeature_RNA > 24 & nFeature_RNA < 300 & percent.mt < 5)

# run sctransform
#replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
gc()
d <- SCTransform(d, ncells = 5500, n_genes = 550, variable.features.n = 550, conserve.memory = TRUE, verbose = TRUE)
gc()

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(d)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(d)

#Perform linear dimensional reduction
d <- RunPCA(d, features = VariableFeatures(object = d))

# Examine and visualize PCA results a few different ways
print(d[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(d, dims = 1:2, reduction = "pca")

DimPlot(d, reduction = "pca") + NoLegend()

DimHeatmap(d, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(d, dims = 1:15, cells = 500, balanced = TRUE)

saveRDS(d, file = "./merfish.rds")

d <- RunUMAP(d, dims = 1:30)
d <- FindNeighbors(d, reduction = "pca", dims = 1:30)
d <- FindClusters(d, resolution = 0.3)
