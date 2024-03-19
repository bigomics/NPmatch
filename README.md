# NPmatch
Latent Batch Effect Correction of Omics data by Nearest-Pair Matching

Batch effects (BEs) are a predominant source of noise in omics data and often mask real biological signals. BEs remain common in existing datasets. Current methods for BE correction mostly rely on specific assumptions or complex models, and may not detect and adjust BEs adequately, impacting downstream analysis and discovery power. To address these challenges we developed NPmatch, a nearest-neighbor matching-based method that adjusts BEs satisfactorily and outperforms current methods in a wide range of datasets.

NPmatch (Nearest-Pair Matching) relies on distance-based matching to deterministically search for nearest neighbors with opposite labels, so-called “nearest-pair”, among samples. NPmatch requires knowledge of the phenotypes but not of the batch assignment. Differently to many other algorithms, NPmatch does not rely on specific models or underlying distribution. It does not require special experimental designs, randomized controlled experiments, control genes or batch information. NPmatch is based on the simple rationale that samples sharing a biological state (e.g., phenotype, condition) should empirically pair based on distance in biological profiles, such as transcriptomics profiles.

NPmatch is freely available on GitHub. It's a main batch correction algorithm in OmicsPlayground, our Bioinformatics platform at BigOmics Analytics. In OmicsPlayground, you can perform NPmatch without coding needs.

## Installation
You can install the NPmatch R package with the following steps:
1. Download NPmatch from https://github.com/bigomics/NPmatch or use "git clone" in the command line;
2. Enter the directory where NPmatch has been downloaded;
3. In your terminal, type: "R CMD INSTALL NPmatch" to install NPmatch.

## Usage example
We provide a basic example on how to use NPmatch to correct batch effects (BEs) in a biological dataset.
We use the GSE10846 dataset (Lenz et al., 2008), which include array gene expression profiling on clinical samples from diffuse large B-cell lymphoma (DLBCL) patients pre-treated with the pharmacological regimens CHOP and Rituximab-CHOP. Dataset includes two biological types of DLBCL: 167 ABC and 183 GCB samples. Treatment was performed prior to expression profiling and samples were split in the two treatment groups for processing. Thus, "treatment" represents a batch variable. We show that BEs in the uncorrected data appear evident with samples clustering by pharmacological treatment. Following NPmatch batch correction, the samples cluster by DLBCL type, reflecting their biological heterogeneity.

``` r
## Load NPmatch and limma
library("NPmatch")
library("limma")

## X is the raw data matrix, with features in rows and samples in columns.
## Meta is a matrix or a dataframe with the full set of metadata associated with X. 
## We also need to make sure that the samples in X and Meta are aligned.
X <- read.table("GSE10846.Expression.txt", sep="\t")
Meta <- read.table("GSE10846.Metadata.txt", sep="\t")
dim(X); class(X)
dim(Meta); class(Meta)
table(rownames(Meta) == colnames(X))

## To correct batch effects, NPmatch requires a vector of phenotype labels per sample.
## To assess  batch effect correction, we also need a vector of batch labels (see below).
## "pheno" is the vector of phenotype labels.
## "batch" is the vector of batch labels.
pheno <- Meta[,"dlbcl.type"]
batch <- Meta[,"Chemotherapy"]

## Intra-sample normalization of the raw gene expression matrix with counts-per-million (CPM)
## You can use the normalize.log2CPM.R function provided
nX <- normalize.log2CPM(X)

## Inter-sample normalization by quantile normalization
nX <- limma::normalizeQuantiles(nX)

## Batch correction with NPmatch
cX <- NPmatch(X=nX, y=pheno, dist.method="cor", sdtop=5000)
table(rownames(Meta) == colnames(cX))

## Check batch effects by plotting UMAP or t-SNE coordinates
nb <- max(1, min(30, round(ncol(cX) / 5)))
# pos <- Rtsne::Rtsne(t(cX), perplexity=nb)$Y ## t-SNE
pos <- uwot::tumap(t(cX), n_neighbors = max(2, nb)) 
pos <- data.frame(Dim1=pos[,1], Dim2=pos[,2], Pheno=pheno, Batch=batch)
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

mtext("NPmatch batch-corrected data",
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

### Support
For support feel free to reach our Bioinformatics Data Science Team at BigOmics Analytics:

Antonino Zito, PhD:  antonino.zito@bigomics.ch

Ivo Kwee, PhD: ivo.kwee@bigomics.ch
