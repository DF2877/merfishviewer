## Should come after seurat-pre.R ##

library(Seurat)

#Load data

d <- readRDS(file = "D:/df/abc_sc_300k_pca.rds") #Seurat object, processed, ready for clustering

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(d)

## Stop and set dims ##

#Cluster the cells
dims = 1:12 #Set dims according to elbow plot

d <- FindNeighbors(d, dims = dims)

## Stop and set resolution ##

res <- 0.0008 #Set resolution by trial and error (look at UMAP)

gc()
d <- FindClusters(d, resolution = res) #Set resolution by trial and error (look at UMAP)
#Run non-linear dimensional reduction (UMAP/tSNE)
d <- RunUMAP(d, dims = dims)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(d, reduction = "umap", label = TRUE)
gc()

## Stop and set resolution, then plot again, repeat until good resolution found ##


#Find marker genes for each cluster
#Default is Wilcoxon rank sum test
# markers <- FindAllMarkers(d, min.pct = 0.01, logfc.threshold = 0.3, min.diff.pct = 0.01, max.cells.per.ident = 3000)
# markers <- FindAllMarkers(d)
markers <- FindMarkers(d, min.pct = 0.01, logfc.threshold = 0.3, min.diff.pct = 0.01, max.cells.per.ident = 3000,
                       ident.1 = 0)
gc()
markers <- markers[markers$p_val_adj < 0.05,]

#ENSMUSG00000030500 vglut2  glutamatergic 
#ENSMUSG00000028645 Slc2a1  glutamatergic 
#ENSMUSG00000020178 Adora2a gabaergic 
#ENSMUSG00000061718 Ppp1r1b gabaergic 
#ENSMUSG00000070880 Gad1    gabaergic 
#ENSMUSG00000026787 Gad2    gabaergic
#ENSMUSG00000029819 Npy     gabaergic
#ENSMUSG00000045573 Penk    gabaergic
#ENSMUSG00000053279 Aldh1a1 dopaminergic
#ENSMUSG00000037025 Foxa2   dopamingeric
#ENSMUSG00000024858 Girk2   dopaminergic
#ENSMUSG00000000214 Th      dopaminergic
#ENSMUSG00000055197 Fev     seretonergic
#ENSMUSG00000020838 Slc6a4  seretonergic
#ENSMUSG00000006764 tph2    seretonergic
#ENSMUSG00000021919 Chat    cholinergic
#ENSMUSG00000037771 vgat    anti-target 

"glutamatergic"
markers[c("ENSMUSG00000030500","ENSMUSG00000028645"),]

"gabaergic"
markers[c("ENSMUSG00000020178","ENSMUSG00000061718","ENSMUSG00000070880","ENSMUSG00000026787","ENSMUSG00000029819","ENSMUSG00000045573"),]

"dopaminergic"
markers[c("ENSMUSG00000053279","ENSMUSG00000037025","ENSMUSG00000024858","ENSMUSG00000000214"),]

"seretonergic"
markers[c("ENSMUSG00000055197","ENSMUSG00000020838","ENSMUSG00000006764"),]

"cholinergic"
markers[c("ENSMUSG00000021919"),]

"motor"
markers[c("ENSMUSG00000021919","ENSMUSG00000042258","ENSMUSG00000027967"),]

#Find marker genes for each cluster
#Default is Wilcoxon rank sum test
ids<-data.frame(gl=integer(),ga=integer(),do=integer(),se=integer(),ch=integer(),mo=integer())
posid<-data.frame(gl=integer(),ga=integer(),do=integer(),se=integer(),ch=integer(),mo=integer())
for (idx in as.numeric(levels(Idents(d)))) {
  message("Identifying cluster ",idx)
  check <- FindMarkers(d, min.pct = 0.01, logfc.threshold = 0.3, min.diff.pct = 0.01, max.cells.per.ident = 3000,
                       ident.1 = idx)
  check <- check[check$p_val_adj < 0.05,]
  ids[idx+1,"gl"] <- sum(sign(check[c("ENSMUSG00000030500","ENSMUSG00000028645"),"avg_log2FC"]),na.rm=T)
  ids[idx+1,"ga"] <- sum(sign(check[c("ENSMUSG00000020178","ENSMUSG00000061718","ENSMUSG00000070880","ENSMUSG00000026787","ENSMUSG00000029819","ENSMUSG00000045573"),"avg_log2FC"]),na.rm=T)
  ids[idx+1,"do"] <- sum(sign(check[c("ENSMUSG00000053279","ENSMUSG00000037025","ENSMUSG00000024858","ENSMUSG00000000214"),"avg_log2FC"]),na.rm=T)
  ids[idx+1,"se"] <- sum(sign(check[c("ENSMUSG00000055197","ENSMUSG00000020838","ENSMUSG00000006764"),"avg_log2FC"]),na.rm=T)
  ids[idx+1,"ch"] <- sum(sign(check[c("ENSMUSG00000021919"),"avg_log2FC"]),na.rm=T)
  ids[idx+1,"mo"] <- sum(sign(check[c("ENSMUSG00000021919","ENSMUSG00000042258","ENSMUSG00000027967","ENSMUSG00000023484"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"gl"] <- max(0,sign(check[c("ENSMUSG00000030500","ENSMUSG00000028645"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"ga"] <- max(0,sign(check[c("ENSMUSG00000020178","ENSMUSG00000061718","ENSMUSG00000070880","ENSMUSG00000026787","ENSMUSG00000029819","ENSMUSG00000045573"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"do"] <- max(0,sign(check[c("ENSMUSG00000053279","ENSMUSG00000037025","ENSMUSG00000024858","ENSMUSG00000000214"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"se"] <- max(0,sign(check[c("ENSMUSG00000055197","ENSMUSG00000020838","ENSMUSG00000006764"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"ch"] <- max(0,sign(check[c("ENSMUSG00000021919"),"avg_log2FC"]),na.rm=T)
  posid[idx+1,"mo"] <- max(0,sign(check[c("ENSMUSG00000021919","ENSMUSG00000042258","ENSMUSG00000027967","ENSMUSG00000023484"),"avg_log2FC"]),na.rm=T)
}
rownames(ids) <- levels(Idents(d))
rownames(posid) <- levels(Idents(d))
as.numeric(rownames(posid[rowSums(posid)>1,]))

DimPlot(d, reduction = "umap", label = TRUE)
FeaturePlot(d, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"))
VlnPlot(d, features = c("Prkcd","Arc","Egr4","Fos","Egr1","Junb","Ntsr1","Tacr1"), stack = TRUE)

saveRDS(d, file = paste("./",filedir,"2.rds",sep=""))

chol.mot <- FindMarkers(d, ident.1 = 0,ident.2 = 4)
chol.mot <- chol.mot[chol.mot$p_val_adj < 0.05,]

write.csv(chol.mot,file = paste0("./chol-mot.csv"))

library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(chol.mot)
df$id <- NA
G_list <- getBM(filters= "ensembl_gene_id", attributes="mgi_symbol",values=genes,mart= mart)
write.table(G_list,file = "./genes.txt",row.names = F,quote=F)
