rm(list=ls())
ddir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/'
rdir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/result/'

### filtered lowly variable genes in gene expression matrix
expr <- readRDS(paste0(ddir, 'CD8_saverseuratall_CMsaver3.rds'))
expr <- expr[rowMeans(expr>0.1)>0.01, ]
saveRDS(expr, paste0(ddir, 'CD8_saverseuratall_CMsaver3_filtered.rds'))

### read in data
expr <- readRDS(paste0(ddir, 'CD8_saverseuratall_CMsaver3_filtered.rds'))
d <- readRDS(paste0(ddir, 'CMsaver3_metadata_incpseudotime.rds'))

### prepare data in the right format
d <- d[colnames(expr), ]
pt = d[,1]
names(pt) <- rownames(d)

meta <- d[, c(4,6,2)]
meta <- meta[!duplicated(meta), ]
rownames(meta) <- paste0(meta[,1], '_', meta[,3])
meta[,1] = 1
meta[,2] = ifelse(meta[,2]=='Responder', 1, 0)
meta[,3] = ifelse(meta[,3]=='Post', 1, 0)
design <- as.matrix(meta)
colnames(design) <- c('intercept', 'respond', 'treatment')

cellanno = data.frame(cell = rownames(d), sample = paste0(d[, 4], '_', d[,2]))
length(intersect(names(pt), cellanno[,1]))

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[names(pt), ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]

saveRDS(pt, paste0(ddir, 'pt.rds'))
saveRDS(cellanno, paste0(ddir, 'cellanno.rds'))
saveRDS(design, paste0(ddir, 'design.rds'))


### prepare hdf5 data
source('/home/whou10/scratch16/whou10/trajectory_variability/h5func/saveh5.R')
path <- paste0(ddir, 'hnscc.h5')
saveh5(expr,pt,cellanno,path)




