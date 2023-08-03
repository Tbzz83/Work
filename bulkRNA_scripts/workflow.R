#source("matrix_to_counts.r")

## IMPORTANT: Specify below into the following variables which groups in your 
## condition that you are comparing
comparison_group <- "A" 
control_group <- "B"
source("deseq2.r")
source("important_genes.r")
source("get_gene_id.r")
#source("Box_plots.r")
source("normal_counts_all_groups.r")
# Some plots don't automatically plot as they go, so here are those
plotPCA(vsd, intgroup=c("condition"))
plotPCA(rld, intgroup=c("condition"))
plotPCA(ntd, intgroup=c("condition"))

# box plot is called 'p'
#p

# Adding some plots for some genes we will want to look at
# Positive markers #
plotCounts(dse, gene=df2$name[grep("ENSG00000135318", df2$name)], intgroup="condition", main = 'CD73 (NT5E)')
plotCounts(dse, gene=df2$name[grep("ENSG00000154096", df2$name)], intgroup="condition", main = 'CD90 (THY1)')
plotCounts(dse, gene=df2$name[grep("ENSG00000106991", df2$name)], intgroup="condition", main = 'CD105 (ENG)')

# Negative markers #
plotCounts(dse, gene=df2$name[grep("ENSG00000170458", df2$name)], intgroup="condition", main = 'CD14')
plotCounts(dse, gene=df2$name[grep("ENSG00000174059", df2$name)], intgroup="condition", main = 'CD34')
plotCounts(dse, gene=df2$name[grep("ENSG00000081237", df2$name)], intgroup="condition", main = 'CD45 (PTPRC)')

# Genes of interest
plotCounts(dse, gene=df2$name[grep("ENSG00000150938", df2$name)], intgroup="condition", main = 'CRIM1')
plotCounts(dse, gene=df2$name[grep("ENSG00000120306", df2$name)], intgroup="condition", main = 'CYSTM1')
plotCounts(dse, gene=df2$name[grep("ENSG00000164935", df2$name)], intgroup="condition", main = 'DCSTAMP')
plotCounts(dse, gene=df2$name[grep("ENSG00000185070", df2$name)], intgroup="condition", main = 'FLRT2')
plotCounts(dse, gene=df2$name[grep("ENSG00000110697", df2$name)], intgroup="condition", main = 'PITPNM1')
plotCounts(dse, gene=df2$name[grep("ENSG00000149489", df2$name)], intgroup="condition", main = 'ROM1')



