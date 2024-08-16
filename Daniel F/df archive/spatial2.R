library(Seurat)
library(here)
library(ggplot2)
library(dplyr)
library(schard)
library(Biobase)
library(dbscan)

vizgen.input <- ReadVizgen(data.dir = "D:/df/abc_m/", type = "centroids")

vizgen.input$centroids %>% head()

sce <- schard::h5ad2sce("D:/df/abc_download_root/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw.h5ad")
assayNames(sce)[1] <- "counts"
cell_metadata <- read.csv("D:/df/abc_download_root/metadata/MERFISH-C57BL6J-638850/20231215/cell_metadata.csv")
rownames(cell_metadata) <- cell_metadata$cell_label
sce@colData$cluster_alias <- cell_metadata$cluster_alias[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$average_correlation_score <- cell_metadata$average_correlation_score[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$feature_matrix_label <- cell_metadata$feature_matrix_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_label <- cell_metadata$donor_label[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_genotype <- cell_metadata$donor_genotype[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$donor_sex <- cell_metadata$donor_sex[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$x <- cell_metadata$x[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$y <- cell_metadata$y[match(sce@colData$cell_label,cell_metadata$cell_label)]
sce@colData$z <- cell_metadata$z[match(sce@colData$cell_label,cell_metadata$cell_label)]
vizgen.obj <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)), assay = "Vizgen")
vizgen.obj[["Vizgen"]][[]] <- as.data.frame(rowData(sce))

cents <- CreateCentroids(vizgen.input$centroids)

coords <- CreateFOV(
  coords = cents,
  type = "centroids",
  molecules = NULL,
  assay = "Vizgen"
)

vizgen.obj[["abc"]] <- coords

GetTissueCoordinates(vizgen.obj[["abc"]][["centroids"]]) %>%
  head()

ImageDimPlot(subset(vizgen.obj, subset=z==0.4), fov = "abc", cols = "polychrome", axes = TRUE, flip_xy = FALSE)

vizgen.obj <- NormalizeData(vizgen.obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  ScaleData() 
vizgen.obj <- RunPCA(vizgen.obj, features = rownames(vizgen.obj))

vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:13)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:13)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.13)

saveRDS(vizgen.obj, file = "D:/df/merfish2.rds")

vizgen.obj <- readRDS("D:/df/merfish2.rds")



pointx <- 5.5
pointy <- 2.25
eps <- 0.25

d <- matchpt(as.matrix(vizgen.obj[[c("x","y")]]),
             as.matrix(data.frame(pointx,pointy))) 
min_row <- rownames(d[d$distance==min(d$distance),])

dat<- vizgen.obj[[c("x","y")]][rownames(d[d$distance<eps,]),]
# those cells' positions
head(dat)

vizgen.input$centroids %>%
  filter(cell %in% rownames(dat))%>%
  ggplot(aes(x=x, y = y)) +
  geom_point() +
  ggforce::geom_circle(aes(x0 = pointx , y0 = pointy, r = eps)) +
  geom_point(aes(x=pointx, y=pointy), color = "red", size = 3) +
  coord_fixed()
