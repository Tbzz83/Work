# 
# Differential expression analysis with the DESeq2 package.
# 
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#

# Load the library.
suppressPackageStartupMessages(library(DESeq2))

# The name of the file that contains the counts.
counts_file = "counts.csv"

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = "design.csv"

# The final result file.
output_file = "results.csv"

# Read the sample file.
colData <- read.csv(design_file, stringsAsFactors=F)

# Turn conditions into factors.
colData$condition = factor(colData$condition)

# The first level should correspond to the first entry in the file!
# Required later when building a model.
colData$condition = relevel(colData$condition, toString(colData$condition[1]))

# Isolate the sample names.
sample_names <- colData$sample

# Read the data from the standard input.
df = read.csv(counts_file, header=TRUE, row.names=1 )


df[is.na(df)] <- 0

df$ADP3 <- as.integer(df$ADP3)


# Trimming genes that dont have a count of at least 10 in X samples. Where X is 
# the smallest number of samples of your groups
#keep <- rowSums(df >= 10) >= min(table(colData$condition))
#df <- df[keep,]

# Created rounded integers for the count data
countData = round(df[, sample_names])

countData[is.na(countData)] <- 0

# Other columns in the dataframe that are not sample information. 
otherCols = df[!(names(df) %in% sample_names)]

#
# Running DESeq2
#

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

# Trimming genes that dont have a count of at least 10 in X samples. Where X is 
# the smallest number of samples of your groups
#keep <- rowSums(counts(dds) >= 10) >= min(table(colData$condition))
#dds <- dds[keep,]
# Can just do this on the original df as is done above. Doesn't seem to create a big difference in output, one gene maybe

# Run deseq
dse = DESeq(dds)

# ______________________________________________________________________________
# Format the results.

res = results(dse, contrast=c("condition", comparison_group, control_group))
# ______________________________________________________________________________

#
# Some helpful plots
#

plotMA(res, ylim=c(-6,6))
plotCounts(dse, gene=which.min(res$padj), intgroup="condition") # gene= specifies which gene you want to compare.
# here gene=which.min(res$padj) is the gene with the smallest adjusted p-value of all genes


### TRANSFORMATIONS
#
## Variance stabilized transformation
#

if(sum( rowMeans( counts(dse, normalized = TRUE)) > 5)< 1000) {
  vsd <- vst(dse, blind=FALSE, nsub = sum( rowMeans( counts(dse, normalized = TRUE)) > 5 ))
} else {
  vsd <- vst(dse, blind=FALSE)
}


head(assay(vsd), 3)

# and regularized log transformation
rld <- rlog(dds, blind=FALSE)


## Plotting the best transformation
# Flat curves may look promising, but not if there are multiple meaningful differences between groups


# this gives log2(n + 1)
ntd <- normTransform(dse)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Heatmap of count matrix
library("pheatmap")
select <- order(rowMeans(counts(dse,normalized=TRUE)),
                decreasing=TRUE)[1:20] # can change this [1:20] to [1:200] to look at top 200 genes
df <- as.data.frame(colData(dse)[,c("condition")])
rownames(df) <- colnames(dse)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames = TRUE,
         cluster_cols=TRUE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)


## Sample to sample differences

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition)#, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## PCA

plotPCA(vsd, intgroup=c("condition"))

plotPCA(rld, intgroup=c("condition"))

plotPCA(ntd, intgroup=c("condition"))


# Custom PCA through GGplot
#library(ggplot2)
#pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed()


#
# The rest of the code is about formatting the output dataframe.
#

# Turn the DESeq2 results into a data frame.
data = cbind(otherCols, data.frame(res))

# Create the foldChange column.
data$foldChange = 2 ^ data$log2FoldChange

# Rename columns to better reflect reality.
names(data)[names(data)=="pvalue"] <-"PValue"
names(data)[names(data)=="padj"] <- "FDR"

# Create a real adjusted pvalue
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Sort the data by PValue to compute false discovery counts.
data = data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR

# Create the additional columns that we wish to present.
data$baseMeanA = 1
data$baseMeanB = 1

# Get the normalized counts.
normed = counts(dse, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

write.csv(normed, file="norm_counts.csv")

# Save normed
saveRDS(normed, file = "normalized_counts_from_deseq.rds")

# Merge the two datasets by row names.
total <- merge(data, normed, by=0)

# Sort again for output.
total = total[with(total, order(PValue, -foldChange)), ]

## A should be control, B comparison

# Sample names for condition A
col_names_A = colData$sample[colData$condition == control_group]

# Sample names for condition B
col_names_B = colData$sample[colData$condition == comparison_group]

# Create the individual baseMean columns.
if (length(dim(total[, col_names_A])) < 2) {
  # Skip the line of code
} else {
  total$baseMeanA = rowMeans(total[, col_names_A])
}

if (length(dim(total[, col_names_B])) < 2) {
  # Skip the line of code
} else {
  total$baseMeanB = rowMeans(total[, col_names_B])
}

# Bringing some sanity to numbers. Round columns to fewer digits.
total$foldChange = round(total$foldChange, 3)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
total$FDR = round(total$FDR, 4)
total$falsePos = round(total$falsePos, 0)

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)

# Rename the first column.
colnames(total)[1] <- "name"

# Reorganize columns names to make more sense.
new_cols = c("name", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
             "log2FoldChange","lfcSE","stat","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]

# Write the results to the standard output.
write.csv(total, file=output_file, row.names=FALSE, quote=FALSE)

# Inform the user.
print("# Tool: DESeq2")
print(paste("# Design: ", design_file))
print(paste("# Input: ", counts_file))
print(paste("# Output: ", output_file))


