library(ggplot2)
library(tidyverse)

design <- read.csv("design.csv")

df <- read.csv("counts_with_id_and_description.csv")

df <- df[order(desc(df$baseMeanB)),]

# Choosing top 5 highest expressed genes
symbols_to_filter <- head(df$SYMBOL, 5)


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
