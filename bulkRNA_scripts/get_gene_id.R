library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(org.Hs.eg.db)
suppressPackageStartupMessages(library(DESeq2))
#df <- as.data.frame(readRDS("normalized_counts_from_deseq.rds"))
#df <- rownames_to_column(df, var = "name")
#df2 <- as.data.frame(readRDS("normalized_counts_from_deseq.rds"))
#df2 <- rownames_to_column(df2, var = "name")


df <- read.csv('results_abs_change_greater_than_2.csv')
df2 <- read.csv('results_abs_change_greater_than_2.csv')
df <- df %>%
  mutate(name = sub("\\..*", "", name))
gene_names <- as.character(mapIds(org.Hs.eg.db, keys=df$name, keytype="ENSEMBL", columns="SYMBOL",column="SYMBOL"))
df$SYMBOL <- gene_names
non_na_count <- sum(!is.na(df$SYMBOL))


write.csv(x=df, file='counts_with_id.csv', row.names = FALSE )




# Adding a description of each gene
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
IDs <- df$SYMBOL
genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = IDs, mart =ensembl)

genedesc <- genedesc%>%
  dplyr::rename(SYMBOL = external_gene_name)


merged_df <- merge(df, genedesc, by = 'SYMBOL', all = TRUE)

df_subset <- merged_df[grepl("membrane|cell surface", merged_df$description, ignore.case = TRUE), ]

write.csv(x=merged_df, file='counts_with_id_and_description.csv', row.names = FALSE )
write.csv(x=df_subset, file='counts_with_id_and_description_membrane_proteins.csv', row.names = FALSE )



