## Should come after seurat-pre.R ##
## You may also want to do spatial-seurat-prep.R first ##

# Load scRNAseq data
d <- readRDS("D:/df/abc_sc_300k.rds")
# Load table with common protein types and associated gene identifiers
targetIDs <- read.csv("D:/df/targets_and_families.csv", skip=1)

# Load MERFISH data, if you want to know which genes are shared
vizgen.obj <- readRDS("D:/df/merfish_whole.rds")

# Get all merfish genes
if (exists("vizgen.obj")) {
  merfishSymbols <- unique(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]])
}

test <- read.csv("D:/df/abc_download_root/metadata/Zhuang-ABCA-1/20231215/gene.csv")
merfishSymbols <- test$gene_symbol

"079 CEA-BST Six3 Cyp26b1 Gaba"
"080 CEA-AAA-BST Six3 Sp9 Gaba"
"082 CEA-BST Ebf1 Pdyn Gaba"
"083 CEA-BST Rai14 Pdyn Crh Gaba"

# knownMarkers <- c("Six3","Cyp26b1","Sp9","Ebf1","Pdyn","Rai14","Crh")
# knownMarkerTypes <- targetIDs$Type[match(knownMarkers,targetIDs$MGI.symbol)]
# 
# knownMarkers <- cbind(knownMarkers, knownMarkerTypes)

# cea <- subset(d, subset = subclass %in% c("079 CEA-BST Six3 Cyp26b1 Gaba","080 CEA-AAA-BST Six3 Sp9 Gaba","082 CEA-BST Ebf1 Pdyn Gaba","083 CEA-BST Rai14 Pdyn Crh Gaba"))

# sbc <- d@meta.data$subclass
# io <- sbc %in% c("079 CEA-BST Six3 Cyp26b1 Gaba","080 CEA-AAA-BST Six3 Sp9 Gaba","082 CEA-BST Ebf1 Pdyn Gaba","083 CEA-BST Rai14 Pdyn Crh Gaba")
# 
# d <- AddMetaData(d, io, col.name = "io")

# cea <- subset(d, subset = io)
# unique(cea@meta.data[["subclass"]])

# cea.hy 28439x117
# cea.ctx 32285x538
# cea.pal 32285x4238
# cea.str 32285x15538

gc()
# Use Wilcoxon rank-sum test to identify differentially expressed genes
markers <- FindMarkers(d, ident.1 = "0386 CEA-BST Ebf1 Pdyn Gaba_3", group.by = "supertype", verbose = TRUE)
gc()
# Label protein types
targets <- targetIDs$Type[match(rownames(markers),targetIDs$Mouse.Ensembl.Gene)]
# Get MGI symbols for easier reading
mgiSymbols <- d@assays[["RNA"]]@meta.data$gene_symbol[match(rownames(markers),d@assays[["RNA"]]@meta.data$gene_identifier)]

if (exists("merfishSymbols")) {
  # Identify genes in both scRNAseq and MERFISH
  merfish <- !is.na(merfishSymbols[match(mgiSymbols,merfishSymbols)])
  markers <- cbind(markers, targets, merfish)
} else if (exists("vizgen.obj")) {
  merfishSymbols <- vizgen.obj[["Vizgen"]]@meta.data$gene_symbol
  # Identify genes in both scRNAseq and MERFISH
  merfish <- !is.na(merfishSymbols[match(mgiSymbols,merfishSymbols)])
  markers <- cbind(markers, targets, merfish)
} else {
  markers <- cbind(markers, targets)
}

rownames(markers) <- make.names(mgiSymbols, unique = TRUE)
gc()

markers <- markers[markers$p_val_adj < 0.005,]

interest <- markers[markers$targets %in% c("gpcr","vgic","lgic","other_ic"),]

write.csv(markers, file = "D:/df/WMB-scRNAseq/m2-0386-markers.csv")
write.csv(interest, file = "D:/df/WMB-scRNAseq/m2-0386-gpcr-ic.csv")

## Identify which clusters a gene is most expressed in

a <- AggregateExpression(d, group.by = "supertype", verbose = TRUE)
amgi <- d@assays[["RNA"]]@meta.data$gene_symbol[match(a[["RNA"]]@Dimnames[[1]],d@assays[["RNA"]]@meta.data$gene_identifier)]
a[["RNA"]]@Dimnames[[1]] <- make.names(amgi, unique = TRUE)

geneOfInterest <- "Gpr150"
nClusters <- 15

a$RNA[geneOfInterest,order(a$RNA[geneOfInterest,], decreasing=TRUE)[1:nClusters]]
