##
#ZINB-WaVE implementation as shown in https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html
##

#Supposed to be very accurate, according to https://www.embopress.org/doi/full/10.15252/msb.20188746
#Too slow to be effective for large dataset

library(zinbwave)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(sparseMatrixStats)
library(edgeR)
library(schard)

# Register BiocParallel Execution
BiocParallel::register(BiocParallel::SnowParam())
options(mc.cores = 4)

#Load data
filedir <- "WMB-10Xv2-OLF-raw"

sce <- schard::h5ad2sce(paste("./",filedir,"/",filedir,".h5ad",sep = ""))
sce

table(colData(sce)$anatomical_division_label)

#Gene filtering
filter <- rowSums(assay(sce)>5)>5
table(filter)

sce <- sce[filter,]

assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
head(vars)

sce <- sce[names(vars)[1:1000],] #Arbitrary number. Recommend >1000

assayNames(sce)[1] <- "counts"

#ZINB-WaVE
#The zinbFit function
gc()
zinb <- zinbFit(sce, verbose = TRUE, K=10, epsilon=1e6)  #Arbitrary number. Recommend 1e6-1e13 (1e12)

sce_zinb <- zinbwave(sce, fitted_model = zinb, K = 10, epsilon=1e6, #Arbitrary number. Recommend 1e6-1e13 (1e12)
                          observationalWeights = TRUE, verbose = TRUE)

W <- reducedDim(sce_zinb)

data.frame(W, bio=colData(sce)$anatomical_division_label
           #,coverage=colData(sce)$Coverage_Type
           ) %>%
  ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

#Differential Expression
weights <- assay(sce_zinb, "weights")

#Differential expression with edgeR
dge <- DGEList(assay(sce_zinb))
dge <- calcNormFactors(dge)

design <- model.matrix(~anatomical_division_label, data = colData(sce))
dge$weights <- weights
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

lrt <- glmWeightedF(fit, coef = 3)
topTags(lrt)
