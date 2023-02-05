library(rhdf5)
ddir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/'
rdir <- '/home/whou10/scratch4/whou10/RajMandal/data/christ/lamian_data/result/'

pt <- readRDS(paste0(ddir,'pt.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))

## hdf5
library(parallel)
source('/home/whou10/scratch16/whou10/trajectory_variability/h5func/01_function.R')
cdir <- '/home/whou10/scratch16/whou10/trajectory_variability/h5func/'
af <- list.files(cdir,pattern = 'multi')
for (f in af) source(paste0(cdir,f))

res <-
  testpt(
    expr = paste0(ddir, 'hnscc.h5'),
    cellanno = cellanno,
    pseudotime = pt,
    design = design,
    testvar = ncol(design),
    test.type = 'Variable',
    demean = FALSE,
    overall.only = F,
    test.method = 'permutation',
    ncores = 48
  )
saveRDS(res,file=paste0(rdir, 'lamian_pm.48cores.rds'))
