library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)

# Whatever pre-processed SingleCellExperiment object that you have
deng <- readRDS("deng-reads.rds")

table(colData(deng)$cell_type2)

deng <- runPCA(deng)

# Basic PCA to start
plotPCA(deng, colour_by = 'cell_type2')


## Using SC3

# estimate k clusters

deng <- sc3_estimate_k(deng)

# print cluster estimate
metadata(deng)$sc3$k_estimation


plotPCA(deng, colour_by = "cell_type1")


# Running SC3
deng <- sc3(deng, ks = 10, biology = TRUE, n_cores = 1)

#----

## tSNE

set.seed(1234)
deng <- runTSNE(deng, exprs_values = "logcounts", perplexity = 30, k = 8)
plotTSNE(deng)

# Still not sure how to colour by k means

#----

## Ultimately you will need to use some form of clustering before moving on to 
## Differencial expression for scRNA data

