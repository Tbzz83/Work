library(ggplot2)
library(tidyverse)

design <- read.csv("design.csv")

df <- read.csv("norm_counts_id.csv")%>%
  mutate(RowMean = rowMeans(dplyr::select(., -c(X, SYMBOL)), na.rm = TRUE))%>%
  dplyr::select(-X.1)


df <- df %>%
  arrange(desc(RowMean))

# Choosing top 5 highest expressed genes
symbols_to_filter <- head(df$SYMBOL, 10)


df <- df%>%
  dplyr::select(design$sample, SYMBOL)%>%
  pivot_longer(cols = design$sample,
               names_to = "sample", 
               values_to = "normalized_count")

merged_df <- merge(df, design, by = "sample", all.x = TRUE)



# Define the values to filter
filtered_df <- subset(merged_df, SYMBOL %in% symbols_to_filter)

p <- ggplot(filtered_df, aes(x=SYMBOL, y=normalized_count, color = condition)) +
  geom_boxplot()+ theme_classic()
p
