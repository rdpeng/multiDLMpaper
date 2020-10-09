#!/home/rpeng/bin/Rscript --no-save

start <- proc.time()

source("poisson-gibbs.R")

gibbsState <- .readRDS("gibbsState.rds")

set.seed(500)

## Run collapsed models 
pgibbs(gibbsState,
       maxit = 80000,
       deleteCache = TRUE
       )

print(proc.time() - start)

