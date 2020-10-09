#!/home/rpeng/bin/Rscript --no-save --no-site-file --no-init-file --default-packages=NULL

start <- proc.time()

adata <- .readRDS("data/all-sites-collapsed.rds")

source("poisson-full-gibbs.R")
source("loadResults.R")

delta <- 1
bdelta <- 0.2

loadResults("COPDAE")

rm(loadResults); gc()

## Run collapsed models 
pfgibbs(adata,
        outcome = "COPDAE",
        pollutant = "Lag(pm25tmean, 0:13)",
        bMLE = b,
        vMLE = v,
        seed = 100,
        maxit = 100000,
        delta = delta,
        bdelta = bdelta)

print(proc.time() - start)

