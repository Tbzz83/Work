library(scRNA.seq.funcs)
library(edgeR)
#library(monocle)
library(MAST)
library(ROCR)
set.seed(1)

DE <- read.table("TPs.txt") # From bulk RNA-seq, these are the differentially expressed genes
notDE <- read.table("TNs.txt") # These are the not differentially expressed genes

GroundTruth <- list(
  DE = as.character(unlist(DE)), 
  notDE = as.character(unlist(notDE))
) # compiling both for comparisons between NA19101 and NA19239 only

molecules <- read.table("molecules.txt", sep = "\t")
anno <- read.table("annotation.txt", sep = "\t", header = TRUE)
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0) > 5;
counts <- data[gkeep,]
# Library size normalization
lib_size = colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules

## Kolmogorov-Smirnov test. Non-parametric. Sensitive to large sample size. Allows us to tell the distance 
## between empirical cummulative distributions of the expression of each gene in two populations

pVals <- apply(
  norm, 1, function(x) {
    ks.test(
      x[group == "NA19101"], 
      x[group == "NA19239"]
    )$p.value
  }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")

# what are the significant differentially expressed (DE) genes)
sigDE <- names(pVals)[pVals < 0.05]
length(sigDE)  

# Number of KS-DE genes
sum(GroundTruth$DE %in% sigDE) 
  
# Number of KS-DE genes that are true DE genes
sum(GroundTruth$notDE %in% sigDE)
  
  
tp <- sum(GroundTruth$DE %in% sigDE) # true positive
fp <- sum(GroundTruth$notDE %in% sigDE) # false positive
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05]) # true negative
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05]) # false negative
tpr <- tp/(tp + fn) # true positive rate
fpr <- fp/(fp + tn) # false positive rate
cat(c(tpr, fpr))  

pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                 names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)


## Long story short you can test a bunch of different statistical methods
## and find out based on the groundtruth example data here which are
## going to give you the best results