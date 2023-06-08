library(scater)
library(SingleCellExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

# Read the matrix and annotation file 
molecules <- read.delim("molecules.txt",row.names=1)
annotation <- read.delim("annotation.txt",stringsAsFactors = T)

umi <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)

# May not need to do below as we don't use ERCC spike ins
altExp(umi,"ERCC") <- umi[grep("^ERCC-",rownames(umi)), ]
umi <- umi[grep("^ERCC-",rownames(umi),invert = T), ]

# Map ENSEMBL IDs to gene symbols
gene_names <- mapIds(org.Hs.eg.db, keys=rownames(umi), keytype="ENSEMBL", columns="SYMBOL",column="SYMBOL")

# Add symbol column mapping gene names to ENSMBL IDs, view with rowdata(umi)
rowData(umi)$SYMBOL <- gene_names
table(is.na(gene_names))

# Remove genes with no symbols (NA)
umi <- umi[! is.na(rowData(umi)$SYMBOL),]
table(is.na(gene_names))


# Check for mitochondrial proteins. If count is 0, it may be that IDs do not contain 'MT' to indicate mitochondrial
# Can experiment using different database:

ensdb_genes <- genes(EnsDb.Hsapiens.v86)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id
is_mito <- rownames(umi) %in% MT_names
table(is_mito)

grep("^MT-",rowData(umi)$SYMBOL,value = T)

# Check for ribosomal proteins
grep("^RP[LS]",rowData(umi)$SYMBOL,value = T)

# ------------------------------------------------------------------------------

## BASIC QC

# Initially want to look at per cell stats

# Per cell metrics
umi_cell <- perCellQCMetrics(umi,subsets=list(Mito=is_mito))

# Per gene metrics
umi_feature <- perFeatureQCMetrics(umi)

# Adding them all back into the metadata 
umi <- addPerCellQC(umi, subsets=list(Mito=is_mito))
umi <- addPerFeatureQC(umi)
colData(umi)

# ^^ Notice that the 'total' column here is signifying total transcript counts per cell (as rows are cells)
# Viewing counts(umi) or 'molecules' gives the total expression matrix, and this is effectively saying that for one column
# of the expression matrix (which would be one individual cell) the sum of all the counts of that column is 
# the total column. We can confirm that this is true. 

sum(molecules$NA19098.r1.A01) # sum = 63322
colData(umi)[1,14] # = 62669 

# Why is the sum here less? Remember we removed any NA genes, so this likely caused the decrease
# We can also test this by re-running the code and skipping the lines where we did this (line 24-26)
# Sure enough, doing it this way ensures that we got the same number! Hooray. So essentially what we 
# are seeing is 63322 - 62669 = 653 removed genes that were NA for this cell.

# ------------------------------------------------------------------------------

# PLOTS

# Since we ultimately want to determine what a significant number of UMI's per 
# cell is, we can do this with some visualization

# This is a histogram of total counts for all cells. 
hist(
  umi$total,
  breaks = 100
)
abline(v = 25000, col = "red")


# This is a histogram of number of unique umi's per cell
hist(
  umi_cell$detected,
  breaks = 100
)
abline(v = 6244, col = "red")
# This is saying that for most cells there are around 8500 umi's. This is  good, the higher the better


# It is difficult to tell by looking where the threshold should be. We will use adaptive threshold to determine points that 
# are more than 3 median absolute deviations (MADs) away from the median in any of the variables we use for QC.

# Detecting sums of features
qc.lib2 <- isOutlier(umi_cell$sum, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")

# Detecting threshold for cells detected per feature
qc.nexprs2 <- isOutlier(umi_cell$detected, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")

# Spike detection
qc.spike2 <- isOutlier(umi_cell$altexps_ERCC_percent, type="higher")
attr(qc.spike2, "thresholds")

# % Mitochondrial threshold
qc.mito2 <- isOutlier(umi_cell$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")

# Each metric is a true false matrix

# Compiling
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2), SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))


# Alternatively this could all be done in one scater command
reasons <- quickPerCellQC(umi_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))

# Let’s add another metadata column that would keep the information about whether a cell is discarded or not
umi$discard <- reasons$discard


# These plots let you compare different things based on if you would keep or discard a cell
plotColData(umi, x="sum", y="subsets_Mito_percent", colour_by="discard")
# ^ For example cells with high MT are typically dead or dying

plotColData(umi, x="sum", y="detected", colour_by="discard")
plotColData(umi, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")

# ------------------------------------------------------------------------------

# Highly expressed genes

# This might be too instense for this computer to run
plotHighestExprs(umi, exprs_values = "counts", 
                 feature_names_to_plot = "SYMBOL", colour_cells_by="detected")




# Below shows how to filter such that we only keep genes which were detected
# (expression value > 1) in 2 or more cells
keep_feature <- nexprs(umi,byrow = TRUE,detection_limit = 1) >= 2
rowData(umi)$discard <- ! keep_feature
table(rowData(umi)$discard)

# adding a logcounts assay
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)


saveRDS(umi, file = "umi.rds")

