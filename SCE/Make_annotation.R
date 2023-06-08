library(tidyverse)

df <- read.table("df_gene_counts_salmon_v2.tsv", sep = "\t", header = TRUE)%>%
  group_by(File)%>%
  select(File)

df <- unique(df$File)%>%
  as.tibble()%>%
  rename(batch = value)%>%
  mutate(batch = gsub("salmon_outdir_RNA-OME-Manuscript-", "", batch))%>%
  separate(batch, into = c("name","treatment","Cell"), remove = FALSE)
  
write_tsv(df, "annotation_v2.txt")
  
  
  
