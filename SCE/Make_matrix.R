library(dplyr)
counts <- read.table("df_gene_counts_salmon_v2.tsv", sep = "\t", header = TRUE)%>%
  as.tibble()%>%
  dplyr::rename(sample_id = File)%>%
  mutate(sample_id = gsub("salmon_outdir_RNA-OME-Manuscript-", "", sample_id))%>%
  filter(gene_biotype == 'protein_coding')%>%
  select(c('sample_id', 'gene_id', 'countsFromAbundanceNo'))%>%
  pivot_wider(names_from = sample_id, values_from = 'countsFromAbundanceNo')%>%
  as.data.frame()
  #remove_rownames()%>%
  #column_to_rownames(var = 'gene_id')%>%
  
  
counts[is.na(counts)] = 0

write_tsv(counts, "molecules_v2.txt")

