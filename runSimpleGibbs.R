######################################################################
## Load the results

start <- proc.time()

source("simple-gibbs.R")
source("loadResults.R")

outcomeName <- "COPDAE"
covModel <- "exp"
DLag <- "DLag 0:13"
sigmaG <- 0.005

dbfile <- paste("stategibbs", gsub(" +", ".", outcomeName), covModel,
                gsub(" +", ".", DLag), as.character(sigmaG),
                sep = "-")

loadResults(outcomeName, model = DLag)

g <- gibbs(b, v, maxit = 25000, verbose = TRUE, covModel = covModel,
           sigmaG = sigmaG, dbfile = dbfile, do.plot = FALSE)

print(proc.time() - start)
