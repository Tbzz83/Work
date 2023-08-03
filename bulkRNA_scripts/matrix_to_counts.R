library(tidyverse)
counts <- read.table("matrix_gene_counts_salmon_no_subsample.tsv", sep = "\t", header = TRUE)%>%
  dplyr::select(-gene_symbol)%>%
  dplyr::rename(name = gene_id)

write.csv(counts, "counts.csv", row.names = FALSE)
