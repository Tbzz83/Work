library(SingleCellExperiment)
library(scater)

# Loading saved umi object from QC.R
umi <- readRDS(file = 'umi.rds')

# Well create another object called umi.qc which has trimmed poor genes and lowly expressed cells
umi.qc <- umi[! rowData(umi)$discard, ! colData(umi)$discard]

# ------------------------------------------------------------------------------
## PCA

umi <- runPCA(umi, exprs_values = "counts")
dim(reducedDim(umi, "PCA"))

# May need to change colour_by depending on what you want to visualize
plotPCA(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")

# Plotting raw data doesn't show us much, as we need to normalize and transform the data

umi <- runPCA(umi, exprs_values = "logcounts_raw")
dim(reducedDim(umi, "PCA"))

plotPCA(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")

#### PLEASE NOTE:
# logcounts_raw is not sufficient for downstream analysis. Will need to use logcounts 
# slot of SingleCellExperiment object, which is log-transformed and normalized by library size 
# (eg. CPM normalization) Not curently sure how to do this so logcounts_raw will be used to 
# demonstrate


# Now lets have a look at things after QC using umi.qc

umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw")
dim(reducedDim(umi.qc, "PCA"))
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

# plotPCA by default selects top 500 most variable genes. We can adjust this using 'ntop':

umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw", ntop = nrow(umi.qc))
dim(reducedDim(umi.qc, "PCA"))
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

# Can also look at say the top 50 genes
umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw", ntop = 50)
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")


# ------------------------------------------------------------------------------

## tSNE Map (t-Distributed Stochastic Neighbor Embedding)
# Alternative to PCA. tSNE is stochastic so running the method multiple times on
# the same dataset will create different plots 
# can use the set.seed() option so you always get the same plot just for demonstration
# Perplexity reflects the number of neighbours used to build the nearest-neighbor network.
# High perplexity -> dense nertwork. Low perplexity -> groups of cells to separate from each other


# Below is before QC as you can see by 'umi'
set.seed(123456)
umi <- runTSNE(umi, exprs_values = "logcounts_raw", perplexity = 130)
plotTSNE(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")

# After QC
set.seed(123456)
umi.qc <- runTSNE(umi.qc, exprs_values = "logcounts_raw", perplexity = 200)
plotTSNE(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")


# ------------------------------------------------------------------------------

## Correlation 
# Lets look at our original QC PC plot
umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw")
dim(reducedDim(umi.qc, "PCA"))
plotPCA(umi.qc, colour_by = "batch", size_by = "sum", shape_by = "individual")

# Lets see if any variables correlate with any of the PCs

logcounts(umi.qc) <- assay(umi.qc, "logcounts_raw")

getExplanatoryPCs(umi.qc, variables = "sum")

plotExplanatoryPCs(umi.qc,variables = "sum") 

# 86% of PC1 explained by sum

logcounts(umi.qc) <- NULL


# Explanatory variables

plotExplanatoryVariables(reads,exprs_values = "logcounts_raw",
                         variables = c("detected","sum","batch",
                                       "individual","altexps_ERCC_percent","subsets_Mito_percent"))


# ------------------------------------------------------------------------------




