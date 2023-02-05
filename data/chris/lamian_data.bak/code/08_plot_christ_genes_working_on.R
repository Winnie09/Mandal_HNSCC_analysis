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
Res <- readRDS(paste0(rdir, 'lamian_chisq.rds'))

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
Res$populationFit <- getPopulationFit(Res, gene = selgene, type = 'variable', num.timepoint = max(Res$pseudotime))
colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- names(Res$pseudotime)
saveRDS(Res$populationFit, 'lamian_chisq_populationFit_selgene.rds')

### -----------------------------
### plot individual gene example
### -----------------------------
g <- names(clu)[clu == 1][1]
g
for (g in selgene){
  pdf(paste0(pdir,g,'.pdf'), width = 2.6, height = 2)
#plotGeneSampleAndPopulation(Res,i,variable='respond')
  plotGene(Res, g, variable = 'respond', continuous = F, axis.text.blank = T, plot.point = FALSE, line.size = 0.3)
  dev.off()
}

### =============================================================
### try to add line types for pre/post: pre: dashed, post: solid
### =============================================================

variable = 'respond'; continuous = F; axis.text.blank = T; plot.point = FALSE; line.size = 0.3;
testobj = Res; gene = g; variable.text = NULL; free.scale = TRUE; facet.sample = FALSE; plot.point = FALSE; line.alpha = 1; point.alpha=1; point.size=0.5; sep = NA; palette = 'Dark2'; ncol = NULL;   cellProp = FALSE; x.lab = 'Pseudotime'; y.lab = 'Expression'; use.palette = FALSE; legend.position = 'right'; axis.text.range.only = FALSE; ylim = NA; xlim = NA; axis.size = 8; title.size = 8; axis.text.size = 8

library(splines)
library(ggplot2)
library(gridExtra)
library(viridis)   
library(RColorBrewer) ##
if ('expr.ori' %in% names(testobj)) expression <- testobj[['expr.ori']] else expression <- testobj[['expr']]
pseudotime <- testobj[['pseudotime']]
cellanno <- testobj[['cellanno']]
colnames(cellanno) <- c('Cell', 'Sample')
predict.values <- predict_fitting(testobj, gene= gene, test.type = testobj$test.type)
pseudotime = pseudotime[colnames(expression)]
cellanno <- cellanno[match(colnames(expression), cellanno[,1]), ]
# predict.values <- predict.values[, colnames(expression),drop=F]
knotnum <- testobj$knotnum
knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
design <- testobj[['design']]
cellanno <- data.frame(Cell = as.character(cellanno[,1]), 
                       Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
variable.d <- if(is.null(variable)) 1 else variable
if (!is.null(variable.text) & variable.d != 1) {
  design[,variable.d] <- ifelse(design[, variable.d] == 0, variable.text[1], variable.text[2])
}
a <- if (free.scale) 'free' else 'fixed'
if (length(gene) == 1){
  print('plotting one gene ...')
  # pd <- data.frame(expr = expression[gene, ], 
  #                  Sample = cellanno[,2], 
  #                  Variable = design[match(cellanno[,2], rownames(design)), variable.d],  ##
  #                  pseudotime = pseudotime[colnames(expression)])
  # pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
  linedlist <- lapply(unique(cellanno[,2]), function(p){
    tmpcellid <- which(cellanno[,2]==p)
    tmpdf <- data.frame(
      expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][gene, tmpcellid],
      Sample = p,
      Variable = design[rownames(design) == p, variable.d], ##
      pseudotime=pseudotime[tmpcellid])
  })
  ld <- do.call(rbind, linedlist)
  ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
  ld <- ld[order(ld$pseudotime), ] ## add 20200812
  
  p <- ggplot() + 
    geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size)   
  p <- p + 
    theme_classic() +
    xlab(x.lab) + ylab(y.lab) + 
    labs(color = variable) +
    theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), legend.position = legend.position) +
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
  
  p <- p + ggtitle(gene)
  
  if (continuous & !use.palette){
    p <- p + scale_color_viridis(discrete = TRUE, direction = -1)   
  } else {
    if (length(unique(ld[, 'Variable'])) > 8){
      p <- p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))  
    } else {
      p <- p + scale_color_manual(values = rev(brewer.pal(8, palette)[1:length(unique(ld[, 'Variable']))]))
    } 
  }
  
  p <- p + theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), axis.title = element_text(size = axis.size), axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.size), axis.text.y = element_text(size = axis.text.size), plot.title = element_text(size = title.size), legend.position = legend.position) +
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) 
  
}
  
  
