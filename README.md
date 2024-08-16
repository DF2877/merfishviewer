# merfishviewer

## df archive
Scripts that were abandoned or not intended to be reused and left undocumented. Some may be functional, but most are probably not.
## app.R
Prerequisites: spatial-seurat-prep.R
The app lets you view each brain section individually or all sections at once. You can view cell subclasses/supertypes (clusters) and expression data for any of the genes included in the dataset. Expression data can be inverted for negative expression. Any combination of genes and subclasses can be used. 
Additionally, the app lets users select two concentric spherical areas for more detailed analysis. Users can generate violin plots to compare gene expression in those two areas, or can run a full DE analysis using the Wilcoxon rank-sum test. This lets users find differentially expressed genes between a specified brain region and its surrounding cells. There is also an option for strict target selection if known markers or clusters are being used.

## gene-idents.R
Prerequisites: seurat-pre.R, spatial-seurat-prep.R (optional)
Compares gene expression in two groups of cells labeled in metadata. For DE on clustering results, see wilcoxon.R.
Additionally, can identify metadata labels which are most associated with a gene.

## seurat-pre.R
Prepares scRNAseq data for further analysis. Takes many .h5ad that are split from a single dataset. For unsplit datasets, see vanilla-seurat.R. Expects pre-processed data (qc’d, normalized, batch corrected, and feature selected).
Followed by: gene-idents.R (optional), wilcoxon.R (optional)

## spatial-seurat-prep.R
Prepares MERFISH data for usage in app.R. 
Followed by: app.R, gene-idents.R (optional)

## vanilla-seurat.R
Runs differential expression analysis for scRNAseq data. Takes data from a single .h5ad file. For split files, see seurat-pre.R. Does not expect pre-processed data (not qc’d, normalized, etc.), but the processing steps can be skipped. Uses clustering, but can be modified to use metadata. 
This script is pretty much just the Seurat example vignette, hence the name.

## wilcoxon.R
Prerequisites: seurat-pre.R
Clusters scRNAseq data and then runs differential expression analysis using those clusters.

