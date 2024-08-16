library(Seurat)
# library(here)
library(ggplot2)
# library(dplyr)
library(schard)
library(Biobase)
# library(dbscan)
library(scCustomize)
library(scattermore)
library(rgl)

## Unused ##
# Path to the directory with Vizgen MERFISH files; requires at least one of the following files present:
# “cell_by_gene.csv”: used for reading count matrix
# “cell_metadata.csv”: used for reading cell spatial coordinate matrices
# “detected_transcripts.csv”: used for reading molecule spatial coordinate matrices
# I only had cell_metadata, since my count matrix was not a csv
# vizgendir <- "D:/df/abc_m/" 

# Path to the count matrix. 
# I had a .h5ad file, which was proccessed using schard.
# If you have a .csv or .txt, singleCellTK can read it
countpath <- "D:/df/abc_download_root/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw.h5ad"  

# Path to metadata.
# AbcProjectCache provided another, more comprehensive, metadata file in a different cache, so I am using that
metadatapath <- "D:/df/abc_download_root/metadata/WMB-neighborhoods/20231215/views/merfish_cell_metadata_with_group_membership.csv"  

# Load metadata
cell_metadata <- read.csv(metadatapath, colClasses = "character")
rownames(cell_metadata) <- cell_metadata$cell_label

# Load spatial data
# vizgen.input <- ReadVizgen(data.dir = vizgendir, type = "centroids")
vizgen.input <- list()
vizgen.input$centroids <- data.frame(x=as.numeric(cell_metadata$x), y=as.numeric(cell_metadata$y), cell=cell_metadata$cell_label)

vizgen.input$centroids %>% head()

# Loading the counts as an SCE
sce <- schard::h5ad2sce(countpath)
assayNames(sce)[1] <- "counts"

# Attach metadata to the SCe, while also matching cells to put the rows in the same order
# Missing values will be filled with NA
# Certain metadata are required for the app to work. Some are supported and used, but not required.

# Required
sce@colData$cell_label <- cell_metadata$cell_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$x <- as.numeric(cell_metadata$x[match(sce@colData$cell_label,cell_metadata$cell_label)])
sce@colData$y <- as.numeric(cell_metadata$y[match(sce@colData$cell_label,cell_metadata$cell_label)])
sce@colData$z <- cell_metadata$z[match(sce@colData$cell_label,cell_metadata$cell_label)]
# Recommended
sce@colData$subclass <- cell_metadata$subclass[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$supertype <- cell_metadata$supertype[match(sce@colData$cell_label,cell_metadata$cell_label)]

# The rest are unnecessary. They will be unused.
# I don't know why you would bother adding these.
# I don't know why I bothered adding these.
# They literally don't do anything
sce@colData$cluster_alias <- cell_metadata$cluster_alias[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$average_correlation_score <- cell_metadata$average_correlation_score[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$feature_matrix_label <- cell_metadata$feature_matrix_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_label <- cell_metadata$donor_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_genotype <- cell_metadata$donor_genotype[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_sex <- cell_metadata$donor_sex[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$brain_section_label <- cell_metadata$brain_section_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$neurotransmitter <- cell_metadata$neurotransmitter[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$class <- cell_metadata$class[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$cluster <- cell_metadata$cluster[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$neurotransmitter_color <- cell_metadata$neurotransmitter_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$class_color <- cell_metadata$class_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$subclass_color <- cell_metadata$subclass_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$supertype_color <- cell_metadata$supertype_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$cluster_color <- cell_metadata$cluster_color[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$HY.EA.Glut.GABA <- cell_metadata$HY.EA.Glut.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$MB.HB.CB.GABA <- cell_metadata$MB.HB.CB.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$MB.HB.Glut.Sero.Dopa <- cell_metadata$MB.HB.Glut.Sero.Dopa[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$NN.IMN.GC <- cell_metadata$NN.IMN.GC[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$Pallium.Glut <- cell_metadata$Pallium.Glut[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$Subpallium.GABA <- cell_metadata$Subpallium.GABA[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$TH.EPI.Glut <- cell_metadata$TH.EPI.Glut[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$WholeBrain <- cell_metadata$WholeBrain[match(sce@colData$cell_label,cell_metadata$cell_label)]

# Create a Seurat object using counts and metadata
vizgen.obj <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)), assay = "Vizgen")
vizgen.obj[["Vizgen"]][[]] <- as.data.frame(rowData(sce))
# vizgen.obj <- SetAssayData(object = vizgen.obj, layer = "data", new.data = counts(sce))

# Convert spatial data to FOV for Seurat to use
cents <- CreateCentroids(vizgen.input$centroids)
coords <- CreateFOV(
  coords = cents,
  type = "centroids",
  molecules = NULL,
  assay = "Vizgen"
)

#Add FOV to the Seurat object
vizgen.obj@images[["abc"]] <- coords

# Preview
GetTissueCoordinates(vizgen.obj[["abc"]][["centroids"]]) %>% head()
ImageDimPlot(subset(vizgen.obj, subset=z==3), fov = "abc", cols = "polychrome", axes = TRUE, flip_xy = FALSE)

# A little extra processing
vizgen.obj <- NormalizeData(vizgen.obj, normalization.method = "LogNormalize", scale.factor = 10000) %>% ScaleData() 
# vizgen.obj <- RunPCA(vizgen.obj, features = rownames(vizgen.obj))
gc()
# Split the Seurat object into a list of Seurat objects, one for each section. This takes a while.
# We want to do this now so that we can load up individual sections more easily in the future, since we wont always need (or want) the whole brain
sections <- SplitObject(vizgen.obj, split.by = "z")

# Save both for future use
saveRDS(vizgen.obj, file = "D:/df/merfish_whole.rds")
saveRDS(sections, file = "D:/df/merfish_split.rds")

## Move on to app.R ##

# The rest of this is just testing and plotting

# Loading them instead of running the first half of the script
# vizgen.obj <- readRDS("D:/df/merfish_whole.rds")
# sections <- readRDS("D:/df/merfish_split.rds")

# Plotting parameters
pointx <- 6
pointy <- 6
eps <- 0.5
eps2 <- 0.25
sectionn <- 3
cr <- colorRamp(c("white", "red"))
maxx <- ceiling(max(vizgen.obj@meta.data[["x"]], na.rm = TRUE))
minx <- floor(min(vizgen.obj@meta.data[["x"]], na.rm = TRUE))
maxy <- ceiling(max(vizgen.obj@meta.data[["y"]], na.rm = TRUE))
miny <- floor(min(vizgen.obj@meta.data[["y"]], na.rm = TRUE))

# Find cells near a point
distances <- list()
dat <- list()
for (z in names(sections)) {
  distances[[z]] <- matchpt(as.matrix(sections[[z]][[c("x","y")]]),
                            as.matrix(data.frame(pointx,pointy)))
  distances[[z]]$distance <- sqrt(((as.numeric(z)-sectionn)**2)+(distances[[z]]$distance**2))
  dat[[z]] <- distances[[z]][rownames(distances[[z]][distances[[z]]$distance<eps,]),]
  dat[[z]]$io <- dat[[z]]$distance<eps2
}
datdat <- do.call(rbind,unname(dat))
datCells <- rownames(datdat)
localVizgen <- subset(vizgen.obj, subset = cell_label %in% datCells)
localVizgen <- AddMetaData(localVizgen, datdat)

# Plot a point and all nearby cells
vizgen.input$centroids %>%
  filter(cell %in% rownames(dat[[as.character(sectionn)]]))%>%
  ggplot(aes(x=x, y = y)) +
  geom_point() +
  ggforce::geom_circle(aes(x0 = pointx , y0 = pointy, r = eps)) +
  geom_point(aes(x=pointx, y=pointy), color = "red", size = 3) +
  coord_fixed()

# Plot a point and all nearby cells
vizgen.input$centroids %>%
  filter(cell %in% rownames(datdat))%>%
  ggplot(aes(x=x, y = y)) +
  geom_point() +
  ggforce::geom_circle(aes(x0 = pointx , y0 = pointy, r = eps)) +
  geom_point(aes(x=pointx, y=pointy), color = "red", size = 3) +
  coord_fixed()

plot3d( 
  x=datdat$x, y=datdat$y, z=datdat$z, 
  # col = data$color, 
  type = 's', 
  radius = .01,
)

# Plot a section
df <- FetchData(sections[[as.character(sectionn)]],c("x","y","ENSMUSG00000021919"),clean = "all")
DF <- data.frame(x=df[1],df[2],feat=df[3])
DF2 <- data.frame(x=round(df[1],2),y=round(df[2],2),feat=round(df[3]),1)
DF2 <- DF[!duplicated(DF2),]
colnames(DF2) <- c("x","y","Chat")
par(bg="black")
scattermoreplot(DF2$x,DF2$y,pch=".", col = rgb(cr(unlist(DF2[3]) / max(DF2[3])),max=255, alpha=200),size = c(720,720),cex=2, ylim = c(maxy,miny), xlim = c(minx,maxx))
