## SET UP COUNT SIMULATION

set.seed(0)
G <- 5
N <- 30 
M <- 1000
initmean <- 5
initvar <- 10

# Generating random count matrix
mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)

# Assigning row and column names
rownames(mat) <- paste0('gene', 1:M)
colnames(mat) <- paste0('cell', 1:(N*G))
group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
names(group) <- colnames(mat)
heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

# Simulate upregulated or downregulated genes unique to each group
set.seed(0)
upreg <- 5
upregvar <- 10
ng <- 100

diff <- lapply(1:G, function(x) {
  diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
  mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})
names(diff) <- paste0('group', 1:G)

heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


# Let’s also simulate some differentially upregulated genes affecting groups of groups.
# This is typical of cell differentiation processes.

diff2 <- lapply(2:(G-1), function(x) {
  y <- x+G
  diff <- rownames(mat)[(((y-1)*ng)+1):(((y-1)*ng)+ng)]
  mat[diff, group %in% paste0("group", 1:x)] <<- mat[diff, group %in% paste0("group", 1:x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})

heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


## ACTUAL DIFFERENTIAL EXPRESSION

# Consider each group vs. all others. Can use T-test to test whether each gene is significantly 
# Upregulated in each group vs. all others

pv.sig <- lapply(levels(group), function(g){
  ingroup <- names(group)[group %in% g]
  outgroup <- names(group)[!(group %in% g)]
  pv <- sapply(1:M, function(i) {
    t.test(mat[i,ingroup], mat[i,outgroup], alternative='greater')$p.value
    #t.test(mat[i,ingroup], mat[i,outgroup])$p.value
  })
  names(pv) <- rownames(mat)
  pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni
  pv.sig
})
heatmap(mat[unique(unlist(pv.sig)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

perf <- function(pv.sig) {
  predtrue <- unique(unlist(pv.sig)) ## predicted differentially expressed genes
  predfalse <- setdiff(rownames(mat), predtrue) ## predicted not differentially expressed genes
  true <- c(unlist(diff), unlist(diff2)) ## true differentially expressed genes
  false <- setdiff(rownames(mat), true) ## true not differentially expressed genes
  TP <- sum(predtrue %in% true)
  TN <- sum(predfalse %in% false)
  FP <- sum(predtrue %in% false)
  FN <- sum(predfalse %in% true)
  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  prec <- TP/(TP+FP)
  fdr <- FP/(TP+FP)
  acc <- (TP+TN)/(TP+FP+FN+TN)
  return(data.frame(sens, spec, prec, fdr, acc))
}
print(perf(pv.sig))

# USE ANOVA

pv <- sapply(1:M, function(i) {
  mydataframe <- data.frame(y=mat[i,], ig=group)
  fit <- aov(y ~ ig, data=mydataframe)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(pv) <- rownames(mat)
pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni

heatmap(mat[pv.sig,], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

print(perf(pv.sig))


pv.sig.pair <- lapply(levels(group), function(g1) {
  g2 <- setdiff(levels(group), g1)
  unlist(lapply(g2, function(g) {
    ## test two groups vs. all others
    ingroup <- names(group)[group %in% c(g1, g)]
    outgroup <- names(group)[!(group %in% c(g1, g))]
    pv <- sapply(1:M, function(i) {
      t.test(mat[i,ingroup], mat[i,outgroup], alternative='greater')$p.value
    })
    names(pv) <- rownames(mat)
    pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni
    pv.sig
  }))
})

heatmap(mat[unique(unlist(pv.sig.pair)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

print(perf(pv.sig.pair))

pv.sig.trip <- lapply(levels(group), function(g1) {
  g2 <- setdiff(levels(group), g1)
  unlist(lapply(g2, function(gi) {
    g3 <- setdiff(levels(group), gi)
    unlist(lapply(g3, function(gj) {
      ## test three groups vs. all others
      ingroup <- names(group)[group %in% c(g1, gi, gj)]
      outgroup <- names(group)[!(group %in% c(g1, gi, gj))]
      pv <- sapply(1:M, function(i) {
        t.test(mat[i,ingroup], mat[i,outgroup], alternative='greater')$p.value
      })
      names(pv) <- rownames(mat)
      pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni
      pv.sig
    }))
  }))
})

heatmap(mat[unique(unlist(pv.sig.trip)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


print(perf(pv.sig.trip))


## average gene expression for each group
mat.group <- do.call(cbind, lapply(levels(group), function(g) {
  rowMeans(mat[,group==g])
}))
colnames(mat.group) <- levels(group)
## construct tree related groups
hc <- hclust(dist(t(mat.group)), method='complete')
plot(hc)


dend <- as.dendrogram(hc)
pv.sig.all <- c()

pv.recur <- function(dend) {
  g1 <- labels(dend[[1]])
  g2 <- labels(dend[[2]])
  #print(g1)
  #print(g2)
  ingroup <- names(group)[group %in% g1]
  outgroup <- names(group)[group %in% g2]
  pv <- sapply(1:M, function(i) {
    t.test(mat[i,ingroup], mat[i,outgroup])$p.value
  })
  names(pv) <- rownames(mat)
  pv.sig <- names(pv)[pv < 0.05/M/length(hc$height)] ## bonferonni
  #print(pv.sig)
  pv.sig.all <<- c(pv.sig.all, pv.sig) ## save
  
  ## recursion to go down tree if not leaf
  if(!is.leaf(dend[[1]])) {
    pv.recur(dend[[1]])
  }
  if(!is.leaf(dend[[2]])) {
    pv.recur(dend[[2]])
  }
}
pv.recur(dend)

heatmap(mat[unique(unlist(pv.sig.all)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)

print(perf(pv.sig.all))



