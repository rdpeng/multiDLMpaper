library(MCAPSdb)
library(tsModel)

sites <- dget("siteList.R")
sites <- sites[-match("24005", sites)]

data <- lapply(sites, function(site) {cat(site, "\n"); readSite(site) })

source("modeling.R")

postProcess <- function(glmObject) {
    ## Extract coefficients
    cc <- summary(glmObject)$coefficients

    ## Extract covariance matrix
    V <- vcov(glmObject)

    ## Only keep coefficients corresponding to pollutant variable
    i <- grep("pm25tmean", rownames(cc), fixed = TRUE)
    

    if(length(i) == 0) ## No match
        stop("pollutant coefficients not found in model fit")
    rval <- list(beta = cc[i, 1], var = V[i, i])

    if(is.null(rval) || length(rval) == 0)
        stop("problem with return value from 'fitSingleSite'")
    else
        rval
}

source("EMpool.R")

poll <- c("Lag(pm25tmean, 0:3, agecat)", "xLag(pm25tmean, 0:3, agecat)")

fits <- lapply(data, function(d) { cat(as.character(d$fips[1]), "\n"); g <- fitSingleSite0(d, "admitp5", poll, "denomp5"); postProcess(g) })

b <- lapply(fits, "[[", "beta"); v <- lapply(fits, "[[", "var")

e <- MVEMpool(b, v, tol = 1e-6)

plot(e)

######################################################################

d <- .readRDS("results/pm25.dist-lag.rds")
regions <- read.csv("data/fips-region.csv")

library(pooling)

dd <- subset(d, dfTime == 8 & model == "DLag 0:6"
             & outcome == "COPD496"
             & !is.na(beta) & !is.na(var))

## dd <- merge(dd, regions[, c("fips", "region7")], by = "fips")

b <- lapply(dd$beta, function(x) unserialize(x))
v <- lapply(dd$var, function(x) unserialize(x))

e <- dMVEMpool(b, v, tol = 1e-6)


library(lattice)

y <- as.vector(t(e$mhat)) * 1000
x <- rep(0:6, length(b))
std <- sqrt(unlist(lapply(seq(along = v), function(i) diag(e$vhat[,,i])))) * 1000
f <- gl(length(b), 7)
g <- factor(rep(dd$region7, each = 7))
rng <- range(y - 2*std, y + 2*std)
ff <- factor(interaction(sort(g), f[order(g)]))
ff <- factor(ff, levels = sort(levels(ff)))

p <- xyplot(y ~ x | ff, subscripts = TRUE, 
            as.table = TRUE, type = "l", ylim = rng,
            par.strip.text = list(cex = 0.5),
            panel = function(x, y, subscripts, ...) {
                panel.xyplot(x, y, ...)
                llines(x, y - 2*std[subscripts], lty = 2)
                llines(x, y + 2*std[subscripts], lty = 2)
            })
print(p)



bm <- do.call("rbind", b)
xqr <- qr(bm)
Ri <- solve(qr.R(qr))

vo <- lapply(v, function(x) crossprod(Ri, x) %*% Ri)
bo <- lapply(b, function(x) drop(crossprod(Ri, x)))

e <- MVEMpool(bo, vo, tol = 1e-7)


######################################################################
## BDLMs

library(pooling)
source("BDLMapprox.R")

d <- .readRDS("results/pm25.dist-lag.rds")

dd <- subset(d, dfTime == 8 & model == "DLag 0:27"
             & outcome == "COPDAE"
             & !is.na(beta) & !is.na(var))

b <- lapply(dd$beta, function(x) unserialize(x))
v <- lapply(dd$var, function(x) unserialize(x))

bdlm <- lapply(seq(along = b), function(i) {
    g <- BDLM(b[[i]], v[[i]])
    list(theta = g$av.theta, var = g$av.cov)
})

bd <- lapply(bdlm, "[[", "theta")
vd <- lapply(bdlm, "[[", "var")


wmean <- function(x, v) {
    xm <- do.call("rbind", x)
    vinv <- lapply(v, solve)
    x1 <- lapply(seq(along = x), function(i) {
        crossprod(x[[i]], vinv[[i]])
    })
    va <- array(unlist(vinv), c(dim(vinv[[1]]), length(vinv)))
    svinv <- apply(va, c(1, 2), sum)
    x2 <- colSums(do.call("rbind", x1))
    V <- solve(svinv)

    list(b = drop(crossprod(x2, V)),
         v = V)
}

w <- wmean(bd, vd)

with(w, {
    y <- b*1000
    std <- sqrt(diag(v) + diag(h)) * 1000
    rng <- range(y - 2*std, y + 2*std)
    xpts <- 0:27
    plot(xpts, y, ylim = rng, xlab = "Lag", ylab=expression(hat(beta)%*%1000),
         type = "b", pch = 20)
    lines(xpts, y - 2*std, lty = 2)
    lines(xpts, y + 2*std, lty = 2)
    abline(h = 0, lty = 3)
})













dmvnorm <- function(x, mu, Sigma) {
    dS <- determinant(Sigma, log = TRUE)$modulus
    qq <- crossprod(x - mu, solve(Sigma)) %*% (x - mu)
    drop(-0.5 * (dS + qq))
}

source("BDLMapprox.R")

makeObj <- function(mu, b) {
    function(p) { 
        V <- dlCov(p, 21, 0.004)
        LL <- sapply(b, function(x) {
            dmvnorm(x, mu, V)
        })
        -sum(LL)
    }
}

######################################################################
## DL models

d <- .readRDS("results/pm25.dist-lag.rds")

library(pooling)
outcomes <- c("COPD", "COPDAE",
              "ischemic heart disease", "heart failure")

r <- lapply(outcomes, function(out) {
    cat(out, "\n")
    dd <- subset(d, dfTime == 8 & model == "DLag 0:6"
                 & outcome == out
                 & !is.na(beta) & !is.na(var))
    
    b <- lapply(dd$beta, function(x) unserialize(x))
    v <- lapply(dd$var, function(x) unserialize(x))
    MVEMpool(b, v, tol = 1e-7)
})

rng <- range(unlist(lapply(r, function(x) with(x, range(mu - 2*sqrt(diag(Vmu)), mu + 2*sqrt(diag(Vmu))))))) * 1000

total <- round(sapply(r, function(x) sum(x$mu)) * 1000, 1)
totalSD <- round(sapply(r, function(x) sqrt(sum(x$Vmu))) * 1000, 1)

## png(file = "meetings/medicare-retreat-2006-02-14/DL0-6.png", width = 600, height = 600)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 2) + .1, las = 1)
for(i in seq(along = r)) {
    x <- r[[i]]
    std <- sqrt(diag(x$Vmu))
    plot(0:6, x$mu*1000, type = "b", pch = 20,
         xlab = "Lag", ylim = rng,
         ylab = expression("% incr. in hosp. with 10 " * mu * g/m^3 * " incr. in " *PM[2.5]),
         main = outcomes[i])
    abline(h = 0, lty = 2)
    lines(0:6, 1000*(x$mu - 2*std), lty = 2)
    lines(0:6, 1000*(x$mu + 2*std), lty = 2)
    text(5, 2, label = paste(total[i], " (", totalSD[i], ")", sep = ""))
}
         
## dev.off()



######################################################################
## Profile eta

ngridE <- 3
etaMat <- as.matrix(expand.grid(eta1 = seq(-0.7, -0.05, length = ngridE),
                                eta2 = seq(-0.7, -0.05, length = ngridE)))
etaMat <- etaMat[nrow(etaMat):1, ]

em.results <- lapply(seq(nrow(etaMat)), function(i) {
    cat("Eta:", etaMat[i, ], "\n")
    e <- iterMu(b, v, maxit = 500, tol = 1e-6, verbose = FALSE, eta = etaMat[i, ])
    plotMu(e)
    e
})


######################################################################

library(filehash)
source("gibbs-plots.R")

## COPDAE
dirs <- file.path("cluster", formatC(seq(1,18,2), flag="0", width = 2))

## HF
dirs <- file.path("cluster", formatC(seq(2,18,2), flag="0", width = 2))

system.time({
    tot <- lapply(dirs, function(d) {
        cat(d, "\n")
        db <- dbInit(file.path(d, "stategibbs"))
        r <- lapply(db, function(x) list(mu = x$mu))
        unname(r)
    })
})
