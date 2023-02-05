rm(list=ls())
library(here)
library(ggplot2)
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/'
rdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/result/')
pdir <- paste0('/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot/')
dir.create(pdir, recursive = T, showWarnings = F)
Res <- readRDS(paste0(rdir, 'lamian_chisq.rds'))

statistics = Res[[1]]
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
str(diffgene)

## --------------
## population fit
## --------------
Res$expr <- readRDS(paste0(ddir, 'CD8_saverseuratall_CMsaver3_filtered.rds'))
Res$testvar = ncol(Res$design)
# Res$design = Res$design[, c(1, 39, 2:38)]

## important !!!!
pt <- sort(Res$pseudotime)
pt2 <- seq(1, length(pt))
names(pt2) <- names(pt)
Res$pseudotime <- pt2
Res$cellanno <- Res$cellanno[names(Res$pseudotime), ]
## --------------
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable', num.timepoint = max(Res$pseudotime))
colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- names(Res$pseudotime)
saveRDS(Res$populationFit, 'lamian_chisq_populationFit.rds')

## -----------
## clustering
## -----------
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = T, num.timepoint = max(Res$pseudotime))  
XDEType <- getXDEType(Res)
Res$XDEType <- XDEType

## autoclu
# clu <- clusterGene(Res, gene = names(XDEType)[!XDEType %in% c('nonXDE')], type = 'variable', scale.difference = F, method = 'kmeans', k.auto = TRUE)
clu <- clusterGene(Res, gene = names(Res$XDEType)[!Res$XDEType %in% c('nonXDE')], type = 'variable', scale.difference = T, method = 'kmeans', k.auto = F, k = 2) 
sink(paste0(pdir, '/lamian_chisq_XDEType_table.txt'))
print(table(XDEType))
table(clu)
sink()

# 47  54 143 161 207 454 249 
Res$cluster = clu

## --------------
## save diff gene
## --------------
gd <- apply(Res$covariateGroupDiff,1, max) - apply(Res$covariateGroupDiff, 1, min)
allg <- diffgene
res <- data.frame(gene = allg, statistics[allg, ],cluster = Res$cluster[allg], 
                  effect_size = gd[allg], stringsAsFactors = F)
res <- res[order(res[,2], -res[,4]), ]
res <- cbind(res, XDEType = XDEType[rownames(res)])
write.csv(res, paste0(pdir, 'lamian_chisq_idifferential_genes.csv'))

## -----------------------
## plotClusterMeanAndDiff
## -----------------------
pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'variable')
saveRDS(goRes, paste0(pdir, '/lamian_chisq_goRes.rds'))

nn <- sapply(names(goRes), function(i){
  tmp <- goRes[[i]]
  tmp <- tmp[tmp[, 'FDR'] < 0.25, ]
  write.csv(tmp, paste0(pdir, 'cluster', i, '_GO_FDR0.25.csv'))
  print(dim(tmp))
  tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(pdir, 'cluster', i, '_GO_FDR0.05.csv'))
  print(dim(tmp))
  return(0)
})

pdf(paste0(pdir, '/lamian_chisq_hm_GO_term5.pdf'), width = 6.8, height = 3.5)
print(plotGOEnrich(goRes))
dev.off()

pdf(paste0(pdir, '/lamian_chisq_hm_GO_term10.pdf'), width = 6.8, height = 7)
print(plotGOEnrich(goRes, n = 10))
dev.off()

pdf(paste0(pdir, '/lamian_chisq_hm_GO_term10_fdr0.25.pdf'), width = 5, height = 2)
print(plotGOEnrich(goRes, n = 10, fdr.cutoff = 0.25))
dev.off()

pdf(paste0(pdir, '/lamian_chisq_hm_GO_term10_sortbyFC_fdr0.25.pdf'), width = 9, height = 7)
print(plotGOEnrich(goRes, n = 10, sortByFDR = F,fdr.cutoff = 0.25))
dev.off()

# --------------------------------------
# compare original and fitted expression
# --------------------------------------
png(paste0(pdir, 'lamian_chisq_DiffFitHm5.png'),width = 2500,height = 2200,res = 300)
plotDiffFitHm5(Res, subsampleCell = TRUE)
dev.off()

Res$expr <- NULL
saveRDS(Res, paste0(rdir, paste0('lamian_chisq_numeric_res_with_clu.rds')))

id = seq(from = 1, to = length(Res$pseudotime), length.out = ncol(Res$populationFit[[1]]))
colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- colnames(Res$expr)[id]

pdf(paste0(pdir, 'lamian_chisq_cluster_mean.pdf'), width = 4.2, height = 3.5)
plotClusterMean(testobj = Res, type = 'Variable', facet = TRUE, cluster = Res$cluster)
dev.off()

######################################################################
## plot cluster mean
## split up cluster to cluster a and b (e.g. cluster 1 to 1a, and 1b)
######################################################################
abtype = unlist(sapply(Res$XDEType[names(Res$cluster)], function(i){
  if (i == 'trendSig') {
    'a'
  } else if (i == 'bothSig'){
    'b'
  } else {
    ''
  }
}))
clutype = paste0(Res$cluster, abtype)
names(clutype) <- names(Res$cluster)
clutype = clutype[!clutype %in% seq(1,6)]

pdf(paste0(pdir, 'lamian_chisq_cluster_mean_abtype.pdf'), width = 3.8, height = 3.8)
plotClusterMean(testobj = Res, cluster = clutype, type = 'Variable', facet = TRUE, facet_scales = 'free', facet_nrow = 7)
dev.off()

pdf(paste0(pdir, 'lamian_chisq_cluster_group_difference.pdf'), width = 2.2, height = 7.5)
plotClusterDiff(testobj = Res,
                            gene = names(Res$cluster),
                            cluster = Res$cluster,
                            each = T,
                            facet_scales = 'fixed',
                            facet_variable = 'cluster',
                            facet_nrow = 8,
                            sep = '',
                            reverse = F)
dev.off()


### -----------------------------
### plot individual gene example
### -----------------------------
g <- names(clu)[clu == 1][1]
g
for (g in names(clu)){
  pdf(paste0(pdir,g,'.pdf'), width = 2.6, height = 2)
#plotGeneSampleAndPopulation(Res,i,variable='respond')
  plotGene(Res, g, variable = 'respond', continuous = F, axis.text.blank = T, plot.point = FALSE, line.size = 0.3)
  dev.off()
}

