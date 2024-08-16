##
#limma-voom implementation as shown in https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
##

#Supposed to be fast, according to https://www.embopress.org/doi/full/10.15252/msb.20188746
#RAM usage scales extremely poorly with dataset size

library(edgeR)

#Read in the counts table
filedir <- "WMB-10Xv2-OLF-raw"

scedata <- schard::h5ad2sce(paste("./",filedir,"/",filedir,".h5ad",sep = ""))
scedata

table(colData(scedata)$anatomical_division_label)

assayNames(scedata)[1] <- "counts"

counts <- counts(scedata)
head(counts)

#Create DGEList object
d <- DGEList(counts)

#2. Preprocessing
#Calculate normalization factors
d <- calcNormFactors(d)
d
#Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.

#Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d), 1, max) < cutoff)
d <- d[-drop,] 
dim(d) # number of genes left

#Derive experiment information from the sample names
snames <- colnames(counts) # Sample names
snames

id1 <- substr(snames, nchar(snames) - 6, nchar(snames) - 4) 
id2 <- substr(snames, nchar(snames) - 2, nchar(snames)) 

group <- interaction(id1,id2)
group

#Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

#3. Voom transformation and calculation of variance weights
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

#4. Fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))

cond1 <- group081.A01
cond2 <- group081.B01

contr <- makeContrasts(cond1 - cond2, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)

#DE Output
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

# logFC: log2 fold change of cond1/cond2
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value (based on t) from test that logFC differs from 0
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# B: log-odds that gene is DE (arguably less useful than the other columns)

#Number of DE genes
length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = paste(filedir, ".txt", sep = ""), row.names = F, sep = "\t", quote = F)
