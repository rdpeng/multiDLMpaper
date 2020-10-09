source("data/fipsNames.R")
source("BDLMapprox.R")
source("loadResults.R")
source("approx-gibbs.R")

loadResults("COPDAE", nkeep = 30)


y <- unlist(b) * 1000
std <- sqrt(unlist(lapply(v, diag))) * 1000

bdlm <- lapply(seq(along=b), function(i) BDLMeta(b[[i]], v[[i]], c(-0.4, -0.4)))
y <- as.vector(sapply(bdlm, "[[", "mu")) * 1000
std <- as.vector(sapply(bdlm, "[[", "std")) * 1000
    
yy <- as.vector(sapply(bdlm, "[[", "av.theta")) * 1000
bf <- as.vector(sapply(bdlm, "[[", "bfac"))
x <- rep(0:13, length(b))
ff <- gl(30, 14, labels = fipsNames(names(b)))
f <- gl(length(bdlm), 100, labels = fipsNames(names(b)))


etaMat <- as.matrix(expand.grid(seq(-0.35, -0.05, length = 10),
                                seq(-0.37, 0, length = 10)))
x1 <- unname(rep(etaMat[,1], length(bdlm)))
x2 <- unname(rep(etaMat[,2], length(bdlm)))

pal <- colorRampPalette(brewer.pal(4, "Reds"))


## Posterior for gamma
## dev.set(2)
levelplot(bf ~ x1 * x2 | f, as.table = T, col.regions = pal(100),
          ncut=30, par.strip.text=list(cex=.7),
          xlab = expression(gamma[1]), ylab = expression(gamma[2]))


## The data
## dev.set(3)
xyplot(y ~ x | ff, as.table=T, type = "l",
       panel=function(x,y,subscripts,...){
           panel.xyplot(x,y,...)
           panel.abline(h=0,lty=3)
           llines(x,y+2*std[subscripts],lty=2)
           llines(x,y-2*std[subscripts],lty=2)
       }, subscripts=T,par.strip.text=list(cex=.7),
       ylab = expression("% increase in admissions for a 10 " * mu * g/m^3 * " increase in " * PM[2.5]), xlab= "Lag",
       ylim=range(y-2*std,y+2*std))

## Smoothed DL functions
## dev.set(4)
xyplot(yy ~ x | ff, as.table=T, type = "l",
       panel = function(x,y,...) {
           panel.abline(h=0,lty=3)
           panel.xyplot(x,y,...)
       }, par.strip.text = list(cex=0.7), xlab = "Lag",
       ylab = expression("% increase in admissions for a 10 " * mu * g/m^3 * " increase in " * PM[2.5]))
       

ebdlm <- lapply(seq(along = b), function(i) BDLMeta(b[[i]], v[[i]], c(-.4,-.5)))
yyy <- unlist(lapply(ebdlm, "[[", "mu")) * 1000

xyplot(yyy ~ x | ff, as.table = TRUE, type = "l",
       panel = function(x, y, ...) {
           panel.abline(h = 0, lty = 3)
           panel.xyplot(x, y, ...)
       }, par.strip.text = list(cex = 0.7),
       xlab = "Lag",
       ylab = expression("% increase in admissions for a 10 " * mu * g/m^3 * " increase in " * PM[2.5]))
