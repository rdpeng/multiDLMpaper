## Simple two level model, no heterogeneity

######################################################################
## Construct constraint covariance matrix

## Use "exp" or "exp2"

dlCov <- function(param, n, sigmasq, model = "exp") {
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

gibbs <- function(b, v, maxit = 5000, verbose = TRUE, sigmaG = 0.005,
                  seed = 200, ngridG = 100, covModel = "exp",
                  dbfile = "stategibbs", do.plot = TRUE) {
    makeKey <- function(i) {
        paste("state", formatC(i, flag = "0", width = ceiling(log10(maxit))),
              sep = "")
    }
    assign("do.plot", do.plot, globalenv())
    assign("verbose", verbose, globalenv())
    
    library(MASS)
    set.seed(seed)
    
    ## Initialize database of results
    library(filehash)
    dbCreate(dbfile, "DB1")
    db <- dbInit(dbfile, "DB1")

    d <- length(b[[1]])  ## dimension
    gammaMat <- switch(covModel,
                       exp = ,
                       delay = seq(0.05, 0.8, length = ngridG),
                       exp2 = seq(0.01, 0.06, length = ngridG)
                       )
    gamma0 <- rep(mean(gammaMat), 2)

    gibbsState <- list(theta = b,
                       mu = rep(0, d), 
                       covModel = covModel,
                       gamma = gamma0,
                       gamma0 = gamma0,
                       sigmaG = sigmaG,
                       n = length(b),
                       d = d,
                       b = b,
                       v = v,
                       gammaMat = gammaMat)
    params <- c("mu", "gamma")

    if(do.plot) {
        graphics.off()
        x11(width = 5, height = 5)
        x11(width = 5, height = 5, xpos = -20) 
    }
    for(i in seq(maxit)) {
        if(verbose) cat("Iteration", i, "\n")

        for(j in seq(along = params)) {
            if(verbose) cat("\tUpdating", params[j], "\t")
            gibbsState[params[j]] <- list(NULL)
            gibbsState <- updateGibbs(gibbsState)
            if(verbose) cat("\n")
        }
        dbInsert(db, makeKey(i), gibbsState[params])
    }
    db
}

updateGibbs <- function(state) {
    use <- which(sapply(state, is.null))
    name <- names(state)[use]
    
    state[[name]] <- switch(name,
                            mu = updateMu(state),
                            gamma = updateGamma2(state)
                            )
    if(name == "gamma")
        state[["gamma0"]] <- state[["gamma"]]
    state
}

updateMu <- function(state) {
    with(state, {
        OG <- dlCov(gamma, d, sigmaG^2, covModel)
        W <- lapply(v, solve)
        thetabar <- rowSums(sapply(seq(along = theta), function(i) {
            drop(W[[i]] %*% theta[[i]])
        }))
        D <- rowsum(do.call("rbind", W), rep(seq(d), n))
        Vw <- solve(D)
        thetaw <- drop(Vw %*% thetabar)
        B <- OG %*% solve(Vw + OG)
        m <- drop(B %*% thetaw)
        V <- (diag(1, d) - B) %*% OG
        
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
        if(verbose) cat("gamma:", gamma.1)
        gamma.1
    })
}

ldmvnorm <- function(x, mean, sigma) {
    if(is.vector(x))
       x <- matrix(x, ncol = length(x))
    
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- as.numeric(determinant(sigma, log = TRUE)$modulus)
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
    logretval
}
