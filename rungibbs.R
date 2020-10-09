######################################################################
## Load the results

start <- proc.time()

source("approx-gibbs.R")
source("loadResults.R")

outcomeName <- "COPDAE"
covModel <- "exp"
DLag <- "DLag 0:13"
sigmaE <- 0.005
sigmaG <- 0.005

dbfile <- paste("stategibbs", gsub(" +", ".", outcomeName), covModel,
                gsub(" +", ".", DLag), as.character(sigmaG),
                sep = "-")

loadResults(outcomeName, model = DLag)

g <- gibbs(b, v, maxit = 30000, verbose = TRUE, covModel = covModel,
           sigmaE = sigmaE, sigmaG = sigmaG, dbfile = dbfile, do.plot = FALSE)

print(proc.time() - start)
