library(scRNA.seq.funcs)
library(scater)
library(scran)
library(sva)
library(batchelor)
library(kBET)
set.seed(1234567)

## We will use the scran-normalization (logNormCounts)
umi    <- readRDS("umi.rds")
umi.qc <- umi[! rowData(umi)$discard, ! colData(umi)$discard]
qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, clusters = qclust)
umi.qc <- logNormCounts(umi.qc)

## Using ComBat
# If balanced experimental design, can use ComBat to eliminate batch effects while
# preserving biological effects by specifying the biological effects using the mod parameter. 
# If unbalance, will throw error

assay(umi.qc, "combat") <- ComBat(logcounts(umi.qc),batch = umi.qc$replicate)

# can also use ComBat to correct for total features as a co-variate
assay(umi.qc, "combat_tf") <- ComBat(logcounts(umi.qc), batch = umi.qc$detected)


# ------------------------------------------------------------------------------

## mnnCorrect (batchelor)
mnn_out <- fastMNN(umi.qc,batch = umi.qc$replicate)
#assay(umi.qc, "mnn") <- assay(mnn_out, 'reconstructed')


# ------------------------------------------------------------------------------
## Comparing effectiveness


# Below will run a PCA of all our different assays so far
for(n in assayNames(umi.qc)) {
  tmp <- runPCA(umi.qc, exprs_values = n, ncomponents = 20)
  
  print(
    plotPCA(
      tmp,
      colour_by = "batch",
      size_by = "detected",
      shape_by = "individual"
    ) +
      ggtitle(n)
  )
} # mnn for some reason is not working. Something to do with the assay as an object


# Compare the RLE for each assay
res <- list()
for(n in assayNames(umi.qc)) {
  res[[n]] <- suppressWarnings(calc_cell_RLE(assay(umi.qc, n)))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)

# ------------------------------------------------------------------------------

# Will also use kBET although this is complex and works by taking kNN networks around
# random cells and tests the number of cells from each batch against a binomial distribution
compare_kBET_results <- function(sce){
  sce <- umi.qc
  indiv <- unique(as.character(sce$individual))
  norms <- assayNames(sce) # Get all normalizations
  results <- list()
  for (i in indiv){ 
    for (j in norms){
      tmp <- kBET(
        df = t(assay(sce[,sce$individual== i], j)), 
        batch = sce$batch[sce$individual==i], 
        heuristic = TRUE, 
        verbose = FALSE, 
        addTest = FALSE, 
        plot = FALSE)
      results[[i]][[j]] <- tmp$summary$kBET.observed[1]
    }
  }
  return(do.call(rbind.data.frame, results))
}

eff_debatching <- compare_kBET_results(umi.qc)
eff_debatching

## GRAPHING BATCH EFFECTS BASED ON NORMALIZATION

library("reshape2")
library("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Individual", "Normalisation")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Individual, Normalisation, fill=kBET)) +  
  geom_tile() +
  scale_fill_gradient2(
    na.value = "gray",
    low = colorset[2],
    mid=colorset[6],
    high = colorset[10],
    midpoint = 0.5, limit = c(0,1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      size = 12, 
      hjust = 1
    )
  ) + 
  ggtitle("Effect of batch regression methods per individual") 





















