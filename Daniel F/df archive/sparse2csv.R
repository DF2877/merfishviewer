library(Matrix)
library(schard)

# Convert a .h5ad to .csv
filedir <- "C57BL6J-638850-log2" #set the file location
outf <- paste("./",filedir,".csv",sep="")
sce <- schard::h5ad2sce(paste0("./",filedir,".h5ad"))
assayNames(sce)[1] <- "counts"
test <- sce@assays@data@listData[["counts"]]
s.sce <- as.matrix(test)
s.sce # shows  (i, j, x)  [columns of a data frame]
write.csv(s.sce, file = outf, row.names=TRUE)

## Move on to spatial-seurat-prep.R ##