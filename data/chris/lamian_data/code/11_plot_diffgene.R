library(reshape2)
library(ggplot2)
rdir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot/pm/'
rdir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/plot/pm/pdf/'
d <- readRDS(paste0(rdir, 'plotdata.rds'))

for (g in names(d)) {
  pd <- melt(d[[g]])
  pd$respond <- ifelse(sub('_.*','',pd[,2])=='respond1','Responder','Non-responder')
  pd$treatment <- ifelse(sub('.*_','',pd[,2])=='treatment1','Post-treatment','Pre-treatment')
  pdf(paste0(pdir,g,'.pdf'),width=4,height=2.3)
  print(ggplot(pd,aes(x=Var1,y=value,col=respond,linetype=treatment)) + geom_line() + theme_classic() + scale_color_manual(values=c('Non-responder'='royalblue','Responder'='brown2'),name='') + scale_linetype_manual(values=c('Post-treatment'='solid','Pre-treatment'='dashed'),name='') + xlab('Pseudotime') + ylab('Gene expression'))
  dev.off()
}



