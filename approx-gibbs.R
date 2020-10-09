######################################################################
## Construct constraint covariance matrix

## Use "exp" or "exp2"

dlCov <- function(param, n, sigmasq, model) {
    varFun <- function(lambda, u) {  ## variance of distributed lag u
        switch(model,
               exp = exp(-lambda * u),
               exp2 = exp(-lambda * u^2),
               delay = exp(-lambda * ifelse(u <= 2, 0, u - 2))
               )
    }
    corFun <- function(gamma, u) {  ## weight function for correlation matrix
        switch(model,
               exp = exp(-gamma * u),
               exp2 = exp(-gamma * u^2),
               delay = exp(-gamma * ifelse(u <= 2, 0, u - 2))
               )
    }
    genCorMat <- function(gamma, n) {
        ## correlation matrix for distributed lag vector
        if(!isTRUE(n > 1))
            stop("argument 'n' should be > 1")
        D <- diag(corFun(gamma, 0:(n - 1)))
        M <- matrix(1, nrow = n, ncol = n)
        I <- diag(1, n)
        cov <- tcrossprod(D) + tcrossprod((I-D) %*% M, I-D)
        cov2cor(cov)
    }
    if(length(param) == 1)
        param <- rep(param, 2)
    
    ## prior covariance matrix for n x 1 distributed lag vector
    ## param is 2 parameters determining covariance of distributed lag vector
    eta1 <- param[1]
    eta2 <- param[2]
    R <- genCorMat(eta2, n)
    V <- diag(sqrt(varFun(eta1, 0:(n - 1))))
    sigmasq * (V %*% R %*% V)
}

######################################################################
## Gibbs sampler


gibbs <- function(b, v, maxit = 5000, verbose = TRUE,
                  sigmaE = 0.005, sigmaG = 0.005, seed = 100,
                  ngridE = 100, ngridG = 60, fixEta = FALSE,
                  covModel = "exp", dbfile = "stategibbs",
                  do.plot = TRUE, deleteCache = FALSE) {
    makeKey <- function(i) {
        paste("state", formatC(i, flag = "0", width = ceiling(log10(maxit))),
              sep = "")
    }
    assign("do.plot", do.plot, globalenv())
    
    library(MASS)
    library(filehash)
    set.seed(seed)
    
    ## Initialize database of results
    if(file.exists(dbfile)) {
        if(deleteCache) {
            file.remove(dbfile)
            cat("removing existing cache file\n")
        }
        else
            stop(sprintf("cache file '%s' already exists", dbfile))
    }
    dbCreate(dbfile, "DB1")
    db <- dbInit(dbfile, "DB1")
            
    d <- length(b[[1]])  ## dimension
    etaMat <- switch(covModel,
                     exp = ,
                     delay = seq(0.2, 0.8, length = ngridE),
                     exp2 = seq(0.01, 0.06, length = ngridE)
                     )
    gammaMat <- switch(covModel,
                       exp = ,
                       delay = seq(0.05, 0.8, length = ngridG),
                       exp2 = seq(0.01, 0.06, length = ngridG)
                       )
    n <- length(b)
    eta0 <- replicate(n, rep(mean(etaMat), 2), simplify = FALSE)
    gamma0 <- rep(mean(gammaMat), 2)
    
    gibbsState <- list(theta = b, mu = rep(0, d), eta = eta0,
                       eta0 = eta0, covModel = covModel,
                       sigmaE = sigmaE, gamma = gamma0, gamma0 = gamma0,
                       sigmaG = sigmaG, n = n, etaMat = etaMat,
                       d = d, b = b, v = v, 
                       gammaMat = gammaMat)
    params <- if(!fixEta)
        c("eta", "theta", "mu", "gamma")
    else
        c("theta", "mu", "gamma")

    if(do.plot) {
        graphics.off()
        x11(width = 5, height = 5)
        x11(width = 5, height = 5, xpos = -20) 
        ## if(!fixEta)
        ##     x11(width = 5, height = 5, xpos = -450)
        x11(width = 5, height = 5, ypos = -100)        
    }
    for(i in seq(maxit)) {
        if(verbose) cat("Iteration", i, "\n")

        for(j in seq(along = params)) {
            if(verbose) cat("\tUpdating", params[j], "\t")
            gibbsState[params[j]] <- list(NULL)
            gibbsState <- updateGibbs(gibbsState)
            cat("\n")
        }
        dbInsert(db, makeKey(i), gibbsState[params])
    }
    db
}

updateGibbs <- function(state) {
    use <- which(sapply(state, is.null))
    name <- names(state)[use]

    ## If there's no error, update the state object.  If there's an
    ## error, don't update the state, just move to the next parameter.
    val <- switch(name,
               theta = updateTheta(state),
               mu = updateMu(state),
               eta = updateEta2(state),
               gamma = updateGamma2(state)
               )
    if(name == "gamma")
        state[["gamma0"]] <- state[["gamma"]]
    if(name == "eta")
        state[["eta0"]] <- state[["eta"]]
    state
}

library(RColorBrewer)
library(lattice)

plotEta <- function(weights, etaMat) {
    dev.set(dev.next())
    pal <- colorRampPalette(brewer.pal(5, "Blues"))
    y <- weights
    x1 <- rep(etaMat[,1], length(weights))
    x2 <- rep(etaMat[,2], length(weights))
    p <- levelplot(y ~ x1 * x2, as.table = TRUE, col.regions = pal(100),
                   xlab = expression(eta[1]), ylab = expression(eta[2]))
    print(p)
}

## Update eta multivariately
updateEta <- function(state) {
    with(state, {
        nr <- nrow(etaMat)
        thetam <- do.call("rbind", theta)
        lweights <- sapply(seq(nr), function(i) {
            L <- ldmvnorm(thetam, mu, dlCov(etaMat[i, ], d, sigmaE^2, covModel))
            sum(L)
        })
        weights <- local({
            e <- exp(lweights - max(lweights))
            e / sum(e)
        })
        if(do.plot) 
            plotEta(weights, etaMat)
        select <- etaMat[sample(nr, 1, prob = weights), ]
        cat("eta:", select)
        select
    })
}

## Update eta one at a time; use common eta across counties
updateEta2.orig <- function(state) {
    with(state, {
        nr <- length(etaMat)
        eta.1 <- eta0
        thetam <- do.call("rbind", theta)

        if(do.plot) {
            dev.set(dev.next())
            par(mfrow = c(2, 1), mar = c(4,4,1,1))
        }
        for(k in 1:2) {  ## cycle over eta
            lweights <- sapply(seq(nr), function(i) {  ## cycle over etaMat
                eta.1[k] <- etaMat[i]
                OE <- dlCov(eta.1, d, sigmaE^2, covModel)
                L <- ldmvnorm(thetam, mu, OE)
                sum(L)
            })
            weights <- local({
                e <- exp(lweights - max(lweights))
                e / sum(e)
            })
            if(do.plot) {
                plot(etaMat, weights,
                     xlab = substitute(eta[s], list(s=k)),
                     type = "h", ylab = "Conditional probability",
                     ylim = c(0, 1))
            }
            eta.1[k] <- sample(etaMat, 1, prob = weights)
        }
        cat("eta:", eta.1)
        eta.1
    })
}

## Update eta one at a time; use county-specific eta
updateEta2 <- function(state) {
    with(state, {
        nr <- length(etaMat)
        eta.1 <- eta0
        
        for(j in seq(along = eta.1)) {
            cat(j)
            for(k in 1:2) {  ## cycle over eta
                lweights <- sapply(seq(nr), function(i) {  ## cycle over etaMat
                    eta.1[[j]][k] <- etaMat[i]
                    OE <- dlCov(eta.1[[j]], d, sigmaE^2, covModel)
                    L <- ldmvnorm(theta[[j]], mu, OE)
                    sum(L)
                })
                weights <- local({
                    e <- exp(lweights - max(lweights))
                    e / sum(e)
                })
                eta.1[[j]][k] <- sample(etaMat, 1, prob = weights)
            }
        }
        cat("\n")
        eta.1
    })
}




updateGamma1 <- function(state) {
    with(state, {
        nr <- length(gammaMat)
        gamma.1 <- gamma0
        
        lweights <- sapply(seq(nr), function(i) {
            gamma.1[varyGamma] <- gammaMat[i]
            ldmvnorm(mu, 0, dlCov(gamma.1, d, sigmaG^2, covModel))
        })
        weights <- exp(lweights - max(lweights))
        weights <- weights / sum(weights)

        if(do.plot) {
            dev.set(dev.next())
            plot(gammaMat, weights,
                 xlab = substitute(gamma[s], list(s=varyGamma)),
                 type = "h", ylab = "Conditional probability")
        }
        select <- gammaMat[sample(nr, 1, prob = weights)]
        cat("gamma:", select)
        select
    })
}

updateGamma2 <- function(state) {
    with(state, {
        nr <- length(gammaMat)
        gamma.1 <- gamma0

        if(do.plot) {
            dev.set(dev.next())
            par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
        }
        for(i in 1:2) {
            lweights <- sapply(seq(nr), function(j) {
                gamma.1[i] <- gammaMat[j]
                ldmvnorm(mu, 0, dlCov(gamma.1, d, sigmaG^2, covModel))
            })
            weights <- exp(lweights - max(lweights))
            weights <- weights / sum(weights)
            
            if(do.plot) {
                plot(gammaMat, weights,
                     xlab = substitute(gamma[s], list(s=i)),
                     type = "h", ylab = "Conditional probability",
                     ylim = c(0, 1))
            }
            select <- gammaMat[sample(nr, 1, prob = weights)]
            gamma.1[i] <- select
        }
        cat("gamma:", gamma.1)
        gamma.1
    })
}

updateMu <- function(state) {
    with(state, {
        OE <- dlCov(eta, d, sigmaE^2, covModel) / n
        OG <- dlCov(gamma, d, sigmaG^2, covModel)
        thetabar <- colMeans(do.call("rbind", theta))
        B <- OG %*% solve(OE + OG)

        V <- (diag(1, d) - B) %*% OG
        m <- drop(B %*% thetabar)

        if(do.plot) {
            dev.set(dev.next())
            plot(seq(0, length(m) - 1), 1000 * m, type = "b",
                 ylim = c(-2, 5), pch = 20,
                 xlab = "Lag", ylab = expression(mu %*% 1000))
            abline(h = 0, lty = 3)
        }        
        drop(mvrnorm(1, m, V))
    })
}

plotTheta <- function(theta) {
    dev.set(dev.next())
    d <- length(theta[[1]])
    n <- length(theta)
    y <- unlist(theta) * 1000
    x <- rep(seq(0, d - 1), n)
    f <- gl(n, d)
    p <- xyplot(y ~ x | f, as.table = TRUE, type = "l",
                panel = function(x, y, ...) {
                    panel.abline(h = 0, lty = 2)
                    panel.xyplot(x, y, ...)
                }, xlab = "Lag",
                ylab = expression(hat(beta) %*% 1000))
    print(p)
}

updateTheta.orig <- function(state) {
    with(state, {
        I <- diag(1, d)
        OE <- dlCov(eta, d, sigmaE^2, covModel)
        r <- lapply(seq(n), function(i) {
            ## B <- OE %*% solve(v[[i]] + OE)
            B <- OE %*% solve(qr(v[[i]] + OE, LAPACK = TRUE))
            V <- (I - B) %*% OE
            mvrnorm(1, mu + drop(B %*% (b[[i]] - mu)), V)
        })
        if(do.plot) {
            dev.set(dev.next())
            matplot(0:13, do.call("cbind", r), type = "l",
                    ylim = c(-.02, .02))
        }
        r
    })    
}

## Use county-specific eta
updateTheta <- function(state) {
    with(state, {
        I <- diag(1, d)
        r <- lapply(seq(n), function(i) {
            browser()
            OE <- dlCov(eta[[i]], d, sigmaE^2, covModel)
            ## B <- OE %*% solve(v[[i]] + OE)
            B <- OE %*% solve(qr(v[[i]] + OE, LAPACK = TRUE))
            V <- (I - B) %*% OE
            mvrnorm(1, mu + drop(B %*% (b[[i]] - mu)), V)
        })
        if(do.plot) {
            dev.set(dev.next())
            matplot(0:13, do.call("cbind", r), type = "l",
                    ylim = c(-.02, .02))
        }
        r
    })    
}

ldmvnorm <- function(x, mean, sigma) {
    if(is.vector(x))
       x <- matrix(x, ncol = length(x))

    sigmaI <- solve(qr(sigma, LAPACK = TRUE))
    distval <- mahalanobis(x, center = mean, cov = sigmaI, inverted = TRUE)
    logdet <- as.numeric(determinant(sigma, log = TRUE)$modulus)
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
    logretval
}



######################################################################

plotLik <- function(mu, sigmaG = 0.05) {
    library(lattice)
    library(RColorBrewer)
    pal <- colorRampPalette(brewer.pal(4, "Greens"))

    mmu <- if(is.list(mu))
        do.call("rbind", mu)
    else
        mu
    
    f <- function(eta) {
        sum(ldmvnorm(mmu, rep(0, length(mu)), dlCov(eta, 14, sigmaG^2)))
    }
    etaMat <- as.matrix(expand.grid(seq(-0.5, -0.1, length = 30),
                                    seq(-0.5, -0.1, length = 30)))
    L <- sapply(seq(nrow(etaMat)), function(i) f(etaMat[i, ]))
    x1 <- etaMat[,1]; x2 <- etaMat[,2]
    p <- levelplot(exp(L-max(L)) ~ x1 * x2, col.regions = pal(100))
    print(p)
}
        
