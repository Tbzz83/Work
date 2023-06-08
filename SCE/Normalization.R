library(scRNA.seq.funcs)
library(scater)
library(scran)


set.seed(1234567)

umi <- readRDS("umi.rds")
umi.qc <- umi[! rowData(umi)$discard, ! colData(umi)$discard]

# Again looking at logcounts_raw PCA plot. We see dependence on sequencing depth
umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw")
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")


# Let's look at CPM-normalized data
logcounts(umi.qc) <- log2(calculateCPM(umi.qc) + 1) # Make new logcounts assay
umi.qc <- runPCA(umi.qc) # Don't need to specify exprs_values b/c logcounts used by default for runPCA
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")


## We can use a RLE plot to see how successful our normalization was

# Logcounts_raw
plotRLE(umi.qc, exprs_values = "logcounts_raw",colour_by = "batch") + ggtitle("RLE plot for logcounts_raw")

# log2(CPM)
plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch") + ggtitle("RLE plot for log2(CPM) counts")

# ------------------------------------------------------------------------------

## Scran-normalized clustering. This assumes not all cells contain the same amount of RNA
qclust <- quickCluster(umi.qc, min.size = 30)
table(qclust)

umi.qc <- computeSumFactors(umi.qc, clusters = qclust)

umi.qc <- logNormCounts(umi.qc)

umi.qc <- runPCA(umi.qc)
plotPCA(umi.qc, colour_by = "batch",size_by = "detected", shape_by = "individual")
plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch")

summary(sizeFactors(umi.qc))
# If size factors are all >0 we are good to go and we should use it. If not try increasing the cluster and pool
# sizes until they are positive.

# ------------------------------------------------------------------------------

## PCA with downsampled data

# Reassigning logcounts assay 
logcounts(umi.qc) <- log2(Down_Sample_Matrix(counts(umi.qc)) + 1)
umi.qc <- runPCA(umi.qc)
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch")


# We have successfully normalized for library size, we can use either method above just
# do some investigating to determine the best one


# ------------------------------------------------------------------------------












