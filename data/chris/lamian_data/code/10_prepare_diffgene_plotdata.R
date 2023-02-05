rm(list=ls())
library(here)
library(ggplot2)
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/'
rdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/result/')
pdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot/pm/')
dir.create(pdir, recursive = T, showWarnings = F)
Res <- readRDS(paste0(rdir, 'lamian_pm.48cores.rds'))

statistics = Res[[1]]
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
str(diffgene)
write.csv(statistics[diffgene, ], paste0(pdir, 'differential_genes_statistics.csv'), row.names = T)
## --------------
## population fit
## --------------
## --------------
## population fit
## --------------
Res$expr <- readRDS(paste0(ddir, 'CD8_saverseuratall_CMsaver3_filtered.rds'))
Res$testvar = 2
# Res$design = Res$design[, c(1, 39, 2:38)]

## important !!!!
pt <- sort(Res$pseudotime)
pt2 <- seq(1, length(pt))
names(pt2) <- names(pt)
Res$pseudotime <- pt2
Res$cellanno <- Res$cellanno[names(Res$pseudotime), ]
## --------------


testobj <- Res
gene <- diffgene
type <- 'variable'
num.timepoint = 1e3
library(splines)
design = testobj$design
knotnum = testobj$knotnum
pseudotime = testobj$pseudotime
pseudotime = pseudotime[order(pseudotime)]
pt <- round(seq(1, max(pseudotime), length.out = min(num.timepoint, max(pseudotime)))) ## downsample


design <- unique(design)
rownames(design) <- paste0('respond',design[,'respond'],'_','treatment',design[,'treatment'])

fitlist <- lapply(gene, function(g){
  # beta <- lapply(testobj$parameter[g], function(i) {
  #   i$beta
  # })
  # names(beta) <- g
  tmp = matrix(testobj$parameter[[g]]$beta, ncol = knotnum[g]+4)
  beta = as.vector(tmp)
  
  x <- sapply(row.names(design), function(i) {
    kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE]) ###
  }, simplify = FALSE)
  
  
  if (knotnum[g] == 0) {
    # phi <- cbind(1, bs(pt))
    phi <- bs(pt, intercept = TRUE)
  } else {
    knots = seq(min(pt), max(pt), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
    # phi <- cbind(1, bs(pt, knots = knots))
    phi <- bs(pt,knots = knots, intercept = TRUE)
  }
  
  
    fit <- sapply(x, function(i) {
      if (ncol(phi) == nrow(i)){
        phi %*% i %*% beta
      } else {
        phi %*% t(i) %*% beta
      }
    })
  return(fit)
})

names(fitlist) <- gene
saveRDS(fitlist,file=paste0(rdir, 'plotdata.rds'))

