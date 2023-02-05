rm(list=ls())
library(here)
library(ggplot2)
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/'
rdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/result/')
pdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot_selgene/')
dir.create(pdir, recursive = T, showWarnings = F)
Res <- readRDS(paste0(rdir, 'lamian_pm.48cores.rds'))

selgene <- c('CCR7','IL7R','CD69','B3GAT1',
             'IFNG',
             'GZMB',
             'GNLY',
             'ITGAE',
             'GZMK',
             'CXCL13'
             ,'PDCD1'
             ,'CTLA4'
             ,'HAVCR2'
             ,'TIGIT'
             ,'LAG3'
             ,'C10orf54'
             ,'BTLA'
             ,'ADORA2A'
             ,'CD28'
             ,'ICOS'
             ,'CD40LG'
             ,'TNFRSF9'
             ,'CD27'
             ,'TNFRSF4'
             ,'TNFRSF18'
             , 'SIRPA'
             ,'CD274'
             ,'CD86'
             ,'LGALS9'
             ,'PVR'
             ,'CD276'
             ,'TNFRSF14'
             ,'CD80'
             ,'ICOSLG'
             ,'CD40'
             ,'TNFSF9'
             ,'CD70'
             ,'TNFSF4'
             ,'TNFSF18'
             ,'CD47'
             ,'PDCD1LG2'
             ,'TCF7',
             'KLRG1')
selgene <- selgene[selgene %in% rownames(Res[[1]])]
write.csv(Res[[1]][selgene, ], paste0(pdir, 'christ_selgene_statistics.csv'))
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
gene <- selgene
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
saveRDS(fitlist,file='/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot_selgene/plotdata.rds')

