cargs <- commandArgs()
infile <- cargs[match("--args", cargs) + 1]

source("gibbs-plots.R")
g <- readMu(infile)

pdf()
plotGibbs(g)
plotEta(g)
plotGamma(g)
dev.off()
