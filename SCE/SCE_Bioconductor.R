library(SingleCellExperiment)
library(ggplot2)
library(scater)

tung_counts <- read.table("molecules_v2.txt", sep = "\t")

tung_annotation <- read.table("annotation_v2.txt", sep = "\t", header = TRUE)

tung <- SingleCellExperiment(
  assays = list(counts = as.matrix(tung_counts)),
  colData = tung_annotation
)

rm(tung_counts, tung_annotation)

assay(tung, "logcounts") <- log2(counts(tung) + 1)
# View logcounts with: logcounts(tung)[1:10, 1:4]

colData(tung)$mean_counts <- colMeans(counts(tung))

colData(tung)$total_counts <- colSums(counts(tung))

assay(tung, "cpm") <- sweep(counts(tung),2,tung$total_counts/1e6,'/')

gene_means <- rowMeans(counts(tung))

tung[gene_means > 0.01, ]

total_detected_per_cell <- colSums(counts(tung) > 0)
# tung[, total_detected_per_cell > 5000] Filter only for cells that have at least 5000 genes

cell_filter <- colSums(counts(tung)) >= 25000

# check how many TRUE/FALSE have
table(cell_filter)

gene_filter <- rowSums(counts(tung) > 5) > ncol(tung)/2

# check how many TRUE/FALSE have
table(gene_filter)

tung_filtered <- tung[gene_filter, cell_filter]

tung_filtered


cell_info <- as.data.frame(colData(tung))

head(cell_info)


# SPECIFY HERE WHAT YOU WANT TO COMPARE, EG. TREATMENT, or REPLICATE
# Manually change x to what you want. Below you can see x = treatment by default

ggplot(data = cell_info, aes(x = treatment, y = total_counts)) +
  geom_violin(fill = 'brown') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# producing the same plot
ggcells(tung, aes(x = treatment, y = total_counts)) + 
  geom_violin(fill = 'orange') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# plotting expression of one of our genes
# Choose a particular gene to look at
ggcells(tung, aes(x = treatment, y = ENSG00000000419.14), exprs_values = "logcounts") + 
  geom_violin(fill = 'coral2') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Adding variance column
colData(tung)$var_counts <- colVars(counts(tung))

ggcells(tung, aes(x = mean_counts, y = var_counts)) + geom_point(aes(colour = treatment)) + theme_bw()

