library(Seurat)
library(schard)
library(scrabbitr)
library(scattermore)

#Load data
filedir <- "D:/df/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-" #set the file location
firstFile <- "P-log2.h5ad"
files <- c("CTXsp-log2.h5ad","HPF-log2.h5ad","HY-log2.h5ad",
           "Isocortex-1-log2.h5ad","Isocortex-2-log2.h5ad",
           "MB-log2.h5ad","MY-log2.h5ad","OLF-log2.h5ad",
           "PAL-log2.h5ad","CB-log2.h5ad","STR-log2.h5ad",
           "TH-log2.h5ad")

sce <- schard::h5ad2sce(paste0(filedir, firstFile))
message(firstFile)
sce <- downsampleSCE(sce, round(ncol(sce)/8))
for (aFile in files) {
  message(aFile)
  temp <- schard::h5ad2sce(paste0(filedir, aFile))
  temp <- downsampleSCE(temp, round(ncol(temp)/8))
  sce <- cbind(sce, temp, deparse.level=1)
}
# sce <- schard::h5ad2sce("D:/df/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-STR-log2.h5ad")
sce

assayNames(sce)[1] <- "counts"

gc()

targetIDs <- read.csv("D:/df/targets_and_families.csv", skip=1)


# Load metadata again
cell_metadata <- read.csv("D:/df/abc_download_root/metadata/WMB-neighborhoods/20231215/views/10x_cell_metadata_with_group_membership.csv")
rownames(cell_metadata) <- cell_metadata$cell_label

# Attach metadata to the SCe, while also matching cells to put the rows in the same order
# Missing values will be filled with NA
sce@colData$cell_barcode <- cell_metadata$cell_barcode[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$barcoded_cell_sample_label <- cell_metadata$barcoded_cell_sample_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$library_label <- cell_metadata$library_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$feature_matrix_label <- cell_metadata$feature_matrix_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$entity <- cell_metadata$entity[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$brain_section_label <- cell_metadata$brain_section_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$library_method <- cell_metadata$library_method[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$region_of_interest_acronym <- cell_metadata$region_of_interest_acronym[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_label <- cell_metadata$donor_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_genotype <- cell_metadata$donor_genotype[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_sex <- cell_metadata$donor_sex[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$x <- cell_metadata$x[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$y <- cell_metadata$y[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$dataset_label <- cell_metadata$dataset_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$cluster_alias <- cell_metadata$cluster_alias[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$neurotransmitter <- cell_metadata$neurotransmitter[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$class <- cell_metadata$class[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$subclass <- cell_metadata$subclass[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$supertype <- cell_metadata$supertype[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$cluster <- cell_metadata$cluster[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$neurotransmitter_color <- cell_metadata$neurotransmitter_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$class_color <- cell_metadata$class_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$subclass_color <- cell_metadata$subclass_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$supertype_color <- cell_metadata$supertype_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$region_of_interest_order <- cell_metadata$region_of_interest_order[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$region_of_interest_color <- cell_metadata$region_of_interest_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$HY.EA.Glut.GABA <- cell_metadata$HY.EA.Glut.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$MB.HB.CB.GABA <- cell_metadata$MB.HB.CB.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$MB.HB.Glut.Sero.Dopa <- cell_metadata$MB.HB.Glut.Sero.Dopa[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$NN.IMN.GC <- cell_metadata$NN.IMN.GC[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$Pallium.Glut <- cell_metadata$Pallium.Glut[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$Subpallium.GABA <- cell_metadata$Subpallium.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$TH.EPI.Glut <- cell_metadata$TH.EPI.Glut[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$WholeBrain <- cell_metadata$WholeBrain[match(sce@colData$cell_label,cell_metadata$cell_label)]

sce@rowRanges@elementMetadata$Type <- targetIDs$Type[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Family.id <- targetIDs$Family.id[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Family.name <- targetIDs$Family.name[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Target.id <- targetIDs$Target.id[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Target.name <- targetIDs$Target.name[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Subunit.id <- targetIDs$Subunit.id[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Subunit.name <- targetIDs$Subunit.name[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Target.systematic.name <- targetIDs$Target.systematic.name[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$Target.abbreviated.name <- targetIDs$Target.abbreviated.name[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]
sce@rowRanges@elementMetadata$synonyms <- targetIDs$synonyms[match(sce@rowRanges@elementMetadata$gene_identifier, targetIDs$Mouse.Ensembl.Gene)]

gc()

#Setup the Seurat Object
d <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
d <- SetAssayData(object = d, layer = "data", new.data = counts(sce))

# Set feature metadata, AKA rowData. Super intuitive, right?
d[["RNA"]][[]] <- as.data.frame(rowData(sce))
d
gc()
#QC and selecting cells for further analysis
# d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

table(d@meta.data[["supertype"]])["1010 DMX VII Tbx20 Chol_1"]

# Specify subclass or supertype
sbc <- d@meta.data$subclass

# Choose subclasses/supertypes of interest
io <- sbc %in% c("079 CEA-BST Six3 Cyp26b1 Gaba","080 CEA-AAA-BST Six3 Sp9 Gaba","082 CEA-BST Ebf1 Pdyn Gaba","083 CEA-BST Rai14 Pdyn Crh Gaba")

d <- AddMetaData(d, io, col.name = "io")

scattermoreplot(na.omit(d@meta.data$x),na.omit(d@meta.data$y),pch=".")

# cea <- subset(d, subset = subclass %in% c("079 CEA-BST Six3 Cyp26b1 Gaba","080 CEA-AAA-BST Six3 Sp9 Gaba","082 CEA-BST Ebf1 Pdyn Gaba","083 CEA-BST Rai14 Pdyn Crh Gaba"))
# unique(cea@meta.data[["subclass"]])
gc()
saveRDS(d, file = "D:/df/abc_sc_P.rds")
gc()

## Move on to gene-idents.R if doing DE on metadata conditions ##

## Continue only if clustering needed ##

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA"), 
                 ncol = 3, pt.size = 0.1, alpha = 0.05)
plot2 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Stop and look ##

lcut <- 800 #set lower cutoff for nFeature_RNA
ucut <- 10000 #set upper cutoff for nFeature_RNA

d <- subset(d, subset = nFeature_RNA > lcut & nFeature_RNA < ucut)

# run sctransform
#replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
gc()
d <- SCTransform(d, ncells = 7500, n_genes = 4500, variable.features.n = 4500, conserve.memory = TRUE, verbose = TRUE)
gc()

#Normalizing the data
# d <- NormalizeData(d)

#Identification of highly variable features (feature selection)
# d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(d)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
# all.genes <- rownames(d)
# gc()
# d <- ScaleData(d)
# gc()

#Perform linear dimensional reduction
d <- RunPCA(d, features = VariableFeatures(object = d))

# Examine and visualize PCA results a few different ways
print(d[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(d, dims = 1:2, reduction = "pca")

DimPlot(d, reduction = "pca") + NoLegend()

# DimHeatmap(d, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(d, dims = 1:15, cells = 500, balanced = TRUE)

saveRDS(d, file = "D:/df/abc_sc_P_pca.rds")

## Move on to Wilcoxon.R ##
