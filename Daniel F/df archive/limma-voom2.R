library(edgeR)

#Read in the counts table
filedir <- "WMB-10Xv2-OLF-raw"

scedata <- schard::h5ad2sce(paste("./",filedir,"/",filedir,".h5ad",sep = ""))
scedata

table(colData(scedata)$anatomical_division_label)

assayNames(scedata)[1] <- "counts"

counts <- counts(scedata)
head(counts)