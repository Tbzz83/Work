library(tidyverse)
df <- read.csv('results.csv')
df <- subset(df, abs(foldChange) > 2)
write.csv(x = df, file = "results_abs_change_greater_than_2.csv",  row.names = FALSE)

df <- df[!is.na(df$PAdj), ]
df$PAdj <- as.numeric(df$PAdj)
df2 <- df[df$PAdj < 0.05, ]

write.csv(x = df2, file = "results_abs_change_greater_than_2_PAdj_less_0.05.csv",  row.names = FALSE)


# ______________________________________________________________________________

