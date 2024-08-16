##
#Implmentation of limma-voom as shown in https://subread.sourceforge.net/RNAseqCaseStudy.html
##

#Supposed to be fast, according to https://www.embopress.org/doi/full/10.15252/msb.20188746
#RAM usage scales extremely poorly with dataset size

# load libraries
library(Rsubread)
library(limma)
library(edgeR)
library(HDF5Array)

#Read in the counts table
filedir <- "WMB-10Xv2-OLF-raw"

data <- schard::h5ad2sce(paste("./",filedir,".h5ad",sep = ""))
data

cond <- factor(colData(data)$anatomical_division_label)
design <- model.matrix(~cond)
table(colData(data)$anatomical_division_label)

assayNames(data)[1] <- "counts"

#Create DGEList object
x <- DGEList(counts(data))

# filter out low-count genes
isexpr <- rowSums(edgeR::cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
gc()
y <- voom(x,design,plot=TRUE)

# cluster libraries
gc()
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)