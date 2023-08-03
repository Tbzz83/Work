library(tidyverse)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

normed <- read.csv('norm_counts.csv')

normed <- normed %>%
  mutate(X = sub("\\..*", "", X))

gene_names <- as.character(mapIds(org.Hs.eg.db, keys=normed$X, keytype="ENSEMBL", columns="SYMBOL",column="SYMBOL"))
normed$SYMBOL <- gene_names

names_to_include <- c("ENSG00000150938", "ENSG00000120306", "ENSG00000164935", "ENSG00000185070",
                      "ENSG00000110697", "ENSG00000149489", "ENSG00000218336",
                      "ENSG00000088726", "ENSG00000126106", "ENSG00000131634", "ENSG00000134825",
                      "ENSG00000137103",
                      "ENSG00000179363", "ENSG00000106609", "ENSG00000072954",
                      "ENSG00000145416", "ENSG00000183654", "ENSG00000173926", "ENSG00000144583",
                      "ENSG00000198060", "ENSG00000136536", "ENSG00000235169", "ENSG00000184785", "ENSG00000178947",
                      "ENSG00000205277", "ENSG00000173702", "ENSG00000169876", "ENSG00000176945",
                      "ENSG00000204544", "ENSG00000169894", 
                      "ENSG00000187054", "ENSG00000185873", "ENSG00000153802", "ENSG00000087128",
                      "ENSG00000079819", "ENSG00000082397", "ENSG00000115109" ,"ENSG00000166947"
)



df3 <- normed[normed$X %in% names_to_include, ]

library(biomaRt)
library(ggplot2)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
IDs <- df3$SYMBOL
genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = IDs, mart =ensembl)

genedesc <- genedesc%>%
  dplyr::rename(SYMBOL = external_gene_name)

design <- as_tibble(read.csv("design.csv"))


merged_df3 <- merge(df3, genedesc, by = 'SYMBOL', all = TRUE)

write.csv(x = merged_df3, file = "normalized_counts_for_important_genes.csv")

merged_df3 <- pivot_longer(merged_df3, cols = c("RNA_CTRL02", "RNA_CTRL01", "ADP1", "ADP2", 
                                                "ADP3", "ADP4", "UCP1", "UCP2", "UCP3", "UCP4", 
                                                "BMP1", "BMP2", "BMP3", "BMP4"),
                           names_to = 'sample', values_to = "normalized_count")

merged_df3 <- left_join(merged_df3, design, by = "sample")


p1 <- ggplot(merged_df3, aes(x = group, y = normalized_count, color = SYMBOL)) +
  geom_point() + facet_wrap( ~ SYMBOL, scales = "free_y")

ggsave(filename = "important_genes.jpg", plot = p1, width = 20)
