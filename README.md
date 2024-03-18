# NPmatch
Latent Batch Effect Correction of Omics data by Nearest-Pair Matching

Batch effects (BEs) are a predominant source of noise in omics data and often mask real biological signals. BEs remain common in existing datasets. Current methods for BE correction mostly rely on specific assumptions or complex models, and may not detect and adjust BEs adequately, impacting downstream analysis and discovery power. To address these challenges we developed NPmatch, a nearest-neighbor matching-based method that adjusts BEs satisfactorily and outperforms current methods in a wide range of datasets.

NPmatch (Nearest-Pair Matching) relies on distance-based matching to deterministically search for nearest neighbors with opposite labels, so-called “nearest-pair”, among samples. NPmatch requires knowledge of the phenotypes but not of the batch assignment. Differently to many other algorithms, NPmatch does not rely on specific models or underlying distribution. It does not require special experimental designs, randomized controlled experiments, control genes or batch information. NPmatch is based on the simple rationale that samples sharing a biological state (e.g., phenotype, condition) should empirically pair based on distance in biological profiles, such as transcriptomics profiles.

# Installation
NPmatch can be installed with the following command, in R:

remotes::install_github("bigomics/NPmatch")

# Usage example
Here we provide a basic example on how to use NPmatch to correct batch effects in a biological dataset (e.g., RNA-seq data)
``` r
## Load NPmatch source code and other needed libraries
library("NPmatch")
library("limma")

## Intra-sample normalization of the raw gene expression matrix with counts-per-million (CPM)
## You can use the normalize.log2CPM.R function provided
nX <- function(X) 

## Inter-sample normalization by quantile normalization
nX <- limma::normalizeQuantiles(nX)

## Batch correction with NPmatch
cX <- NPmatch(X=nX, y=Meta[,"pheno"], dist.method="cor", sdtop=5000)

## Check batch effects by plotting UMAP or t-SNE
nb <- max(1, min(30, round(ncol(cX) / 5)))
# pos <- Rtsne::Rtsne(t(cX), perplexity=nb)$Y ## t-SNE
pos <- uwot::tumap(t(cX), n_neighbors = max(2, nb)) 
pos <- data.frame(Dim1=pos[,1], Dim2=pos[,2], Pheno=Meta[,"pheno"], Batch=Meta[,"batch"])
rownames(pos) <- colnames(cX)
pos[,1:2] <- apply(pos[,1:2], 2, function(x) as.numeric(x))
pos$Col <- as.numeric(factor(pos$Pheno))

## Plot
plot(pos$Dim1,
     pos$Dim2,
     xlab = "Dim1", 
     ylab = "Dim2",
     col = pos$Col,
     pch = 18, 
     cex = 1.8, 
     cex.lab = 2,
     cex.axis = 2,
     las = 1, 
     tcl = -0.1,
     mgp = c(1.5,0.5,0))

mtext("NPmatch corrected data",
      font = 2,
      adj = 0.5,
      cex = 1)

legend("bottomleft",
       unique(pos$Pheno),
       cex = 1.5,
       fill = unique(pos$Col),
       col = unique(pos$Col))

grid(lwd = 1.5)
```
