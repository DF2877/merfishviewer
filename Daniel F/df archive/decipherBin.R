GFP21 <- readRDS("D:/df/2nd_Batch/GFP21.rds")
GFP21_Inh <- readRDS("D:/df/2nd_Batch/GFP21_Inh.rds")
GFP21_dim10reso0.5_15clusters_4510 <- readRDS("D:/df/2nd_Batch/GFP21_dim10reso0.5_15clusters_4510.rds")
GFP2s2 <- readRDS("D:/df/2nd_Batch/GFP2s2.rds")
targetIDs <- read.csv("D:/df/targets_and_families.csv", skip=1)

d <- GFP21_dim10reso0.5_15clusters_4510
DimPlot(d, reduction = "umap", label = TRUE)
FeaturePlot(d, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"))
VlnPlot(d, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"), stack = TRUE)

markers <- FindMarkers(d, ident.1 = 8)
gc()
targets <- targetIDs$Type[match(rownames(markers),targetIDs$MGI.symbol)]
merfish <- !is.na(merfishSymbols[match(rownames(markers),merfishSymbols)])
# mgiSymbols <- d@assays[["RNA"]]@meta.data$gene_symbol[match(rownames(markers),d@assays[["RNA"]]@meta.data$gene_identifier)]

markers <- cbind(markers, targets, merfish)
# rownames(markers) <- make.names(mgiSymbols, unique = TRUE)
gc()

interest <- markers[markers$targets %in% c("gpcr","vgic","lgic","other_ic"),]
interest <- interest[interest$p_val_adj < 0.05,]

write.csv(markers, file = "D:/df/2nd_Batch/m-round2-markers.csv")
write.csv(interest, file = "D:/df/2nd_Batch/m-round2-gpcr-ic.csv")