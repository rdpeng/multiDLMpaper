## Full Poisson likelihood Gibbs sampler

######################################################################
## Construct constraint covariance matrix

## Use "exp" or "exp2"

dlCov <- function(param, n, sigmasq, model, rng, shift) {
    varFun <- function(lambda, u) {  ## variance of distributed lag u
        ## lambda is in [0, 1]
        switch(model,
               exp = exp(-(rng * lambda + shift) * u),
               exp2 = exp(-lambda * u^2),
               delay = exp(-lambda * ifelse(u <= 2, 0, u - 2))
               )
    }
    corFun <- function(gamma, u) {  ## weight function for correlation matrix
        ## gamma is in [0, 1]
        switch(model,
               exp = exp(-(rng * gamma + shift) * u),
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
## Reduce dataset size

expression({
    adata <- .readRDS("data/all-sites.rds")
    bdata <- lapply(adata, function(x) {
        data <- as(x, "data.frame")
        cat(fips <- as.character(data$fips[1]), "\n")
        if(fips == "24005") {
            data$tmpd <- with(data, (tmax35 + tmin35) / 2)
            data$dptp <- data$dptp35
            data$rmtmpd <- with(data, runMean(tmpd, 1:3))
            data$rmdptp <- with(data, runMean(dptp, 1:3))
        }
        data[, c("date", "agecat", "dow", "tmpd", "rmtmpd", "dptp",
                 "rmdptp", "COPDAE", "COPDAEd", "admitp5", "denomp5",
                 "pm25tmean")]
    })
    .saveRDS(bdata, file = "data/all-sites-reduced.rds", compress = TRUE)
    
})


######################################################################
## County-specific time series model

fitSingleSite <- function(data,  ## A data frame for a site
                          outcome,  ## character, name of outcome
                          pollutant,  ## character, name of pollutant variable
                          theta = NULL, ## DL coefficients
                          denom,  ## character, denominator variable name
                          df.Time = 6 * 4,  ## df for smooth function of time
                          df.time = 1 * 4,  ## df for smooth function x age category
                          df.Temp = 6,  ## df for temp smooth function
                          df.Dew = 3,
                          subset = TRUE, singleAgeCat = FALSE,
                          ...) {  
    library(splines)
    library(tsModel)
    if(!is.character(outcome))
        stop("'outcome' should be character")
    if(missing(denom)) {
        denom <- if(length(grep("admit", outcome, fixed = TRUE) > 0))
            sub("admit", "denom", outcome, fixed = TRUE)
        else
            paste(outcome, "d", sep = "")
    }        
    xterms <- c("dow",
                if(singleAgeCat) {
                    NULL
                } else {
                    "agecat"
                },
                paste("ns(tmpd,", df.Temp, ")"),
                paste("ns(rmtmpd,", df.Temp, ")"),
                paste("ns(dptp,", df.Dew, ")"),
                paste("ns(rmdptp,", df.Dew, ")"),
                paste("ns(date,", df.Time, ")"),
                if(singleAgeCat) {
                    NULL
                } else {
                    paste("I(ns(date,",df.time,")*(agecat==\"75p\"))")
                },
                sprintf("offset(log(%s))", denom),
                if(is.null(theta)) {
                    pollutant
                } else {
                    paste("offset(", pollutant, " %*% theta)")
                }
                )
    form <- reformulate(xterms, response = outcome)
    glm(form, data = data, family = poisson, na.action = na.omit,
        control = glm.control(epsilon = 1e-10, maxit = 1000), subset = subset)
}


## Construct a function that evaluates the profile likelihood for
## theta.

makeProfLLik <- function(data, outcome, pollutant) {
    ## data <- as(data, "data.frame")
    fit0 <- fitSingleSite(data, outcome, pollutant)
    
    ## Get non aliased coefficients
    raw.coefs <- coef(fit0)
    aliased <- is.na(raw.coefs)
    coefs <- raw.coefs[!aliased]

    ## Model matrix
    mm <- model.matrix(fit0)
    mm <- mm[, !aliased]
    
    ## Get offset
    mf <- model.frame(fit0)
    off <- model.offset(mf)
    y <- model.response(mf)

    theta.idx <- grep("pm25tmean", names(coefs), fixed = TRUE)
    theta.mm <- mm[, theta.idx]
    beta.mm <- mm[, -theta.idx]

    rm(fit0, data, mm, mf)
    
    function(theta) {
        ## fit <- fitSingleSite(data, outcome, pollutant, theta)
        ## LL <- logLik(fit)
        ## as.numeric(LL)
        g <- glm.fit(beta.mm, y, offset = off + theta.mm %*% theta,
                     family = poisson())
        with(g, rank - aic / 2)
    }
}


## Make function that evaluates the estimated (plug-in) log-likelihood
## for theta.

makeEstLLik <- function(data, outcome, pollutant) {
    data <- as(data, "data.frame")
    ## data <- collapse(data)
    fit0 <- fitSingleSite(data, outcome, pollutant, singleAgeCat = FALSE)

    ## Get terms predictions and remove PM2.5 term
    p <- predict(fit0, type = "terms")
    not.used <- grep("pm25tmean", colnames(p))
    p <- p[, -not.used]

    ## Get constant term (this is what predict.lm does)
    mm <- model.matrix(fit0)
    const <- sum(colMeans(mm) * coef(fit0), na.rm = TRUE)

    ## Get offset
    mf <- model.frame(fit0)
    off <- model.offset(mf)

    ## Sum nuisance terms
    others <- off + rowSums(p) + const
    
    idx.lag <- grep("pm25tmean", colnames(mm), fixed = TRUE)
    lagMat <- sweep(mm[, idx.lag], 2, colMeans(mm[, idx.lag]))
    y <- model.response(mf)

    ## Remove big items from the environment
    rm(fit0, data, mm, mf, p)
    
    function(theta) {
        fit <- tryCatch({
            link <- drop(lagMat %*% theta) + others
            L <- dpois(y, exp(link), log = TRUE)
            sum(L)
        }, interrupt = function(int) {
            int
        }, error = function(e) {
            cat(as.character(e))
            e
        }, warning = function(w) {
            cat(as.character(w))
            w
        })
        if(inherits(fit, "interrupt"))
            stop("computation interrupted")
        if(inherits(fit, "condition"))
            return(-Inf)
        fit
    }
}


## Construct a function that evaluates the log-likelihood for theta
## and beta (nuisance parameters)

makeLL <- function(data, outcome, pollutant) {
    data <- as(data, "data.frame")
    ## data <- collapse(data)  ## Collapse age categories for now
    fit0 <- fitSingleSite(data, outcome, pollutant, singleAgeCat = FALSE)

    ## Get non aliased coefficients
    raw.coefs <- coef(fit0)
    aliased <- is.na(raw.coefs)
    coefs <- raw.coefs[!aliased]

    ## Model matrix
    mm <- model.matrix(fit0)
    mm <- mm[, !aliased]
    
    ## Get offset
    mf <- model.frame(fit0)
    off <- model.offset(mf)
    y <- model.response(mf)

    theta.idx <- grep("pm25tmean", names(coefs), fixed = TRUE)
    theta.mm <- mm[, theta.idx]
    beta.mm <- mm[, -theta.idx]

    beta0 <- coefs[-theta.idx]
    var0 <- vcov(fit0)[names(beta0), names(beta0)]

    rm(fit0, data, mm, mf)

    LL <- function(theta, beta) {
        link <- drop( theta.mm %*% theta + beta.mm %*% beta + off)
        L <- dpois(y, exp(link), log = TRUE)
        sum(L)
    }
    list(LL = LL, theta0 = coefs[theta.idx], beta0 = beta0, var0 = var0)
}

######################################################################
## Main 'pgibbs()' function


pgibbs <- function(sitedata,
                   outcome = "COPDAE",
                   pollutant = "Lag(pm25tmean, 0:13, agecat)",
                   bMLE = NULL,
                   vMLE = NULL,
                   mle.proposal = FALSE,
                   maxit = 5000,
                   sigmaG = 0.005,
                   ngridG = 100,
                   covModel = "exp",
                   d = 14,  ## dimension of theta
                   seed = 100,
                   delta = 1,
                   bdelta = 1,
                   gamma0 = NULL,
                   verbose = TRUE,
                   dbfile = "statepgibbs", 
                   do.plot = FALSE,
                   deleteCache = FALSE) {
    library(MASS)
    library(filehash)

    if(file.exists("random-seed.rds"))
        assign(".Random.seed", .readRDS("random-seed.rds"), globalenv())
    else
        set.seed(seed)

    initGibbsState <- function() {
        s <- new.env(hash = TRUE, parent = globalenv())
        n <- length(sitedata)

        if(mle.proposal) {
            s$bMLE <- bMLE
            s$vMLE <- vMLE
        }
        s$mle.proposal <- mle.proposal
        s$mu <- rep(0, d)
        s$gamma <- gamma0
        s$gammaVals <- gammaVals
        s$sigmaG <- sigmaG

        s$sitedata <- sitedata
        s$outcome <- outcome
        s$pollutant <- pollutant
        
        s$covModel <- covModel
        s$n <- n  ## number of sites
        s$d <- d  ## dimension

        out <- if(file.exists("cacheProfLik.rds")) {
             if(verbose) cat("using cached likelihoods\n")
            .readRDS("cacheProfLik.rds")
        }
        else {
            if(verbose) cat("creating likelihood functions\n")
            pLL <- lapply(seq(along = sitedata), function(i) {
                if(verbose) cat(i, " ")
                data <- sitedata[[i]]
                ## makeEstLLik(data, outcome, pollutant)
                ## makeLL(data, outcome, pollutant)
                makeProfLLik(data, outcome, pollutant)
            })
            if(verbose) cat("\n")
            if(verbose) cat("caching likelihoods for later use\n")
            .saveRDS(pLL, file = "cacheProfLik.rds", compress = TRUE)
           pLL
        }
        s$profileLL <- out
        s$maxit <- maxit
        s$delta <- if(length(delta) == 1) 
            rep(delta, n)
        else {
            if(length(delta) != n)
                stop("'delta' should be length 1 or length n")
            delta
        }
        s$bdelta <- bdelta
        s
    }
    makeKey <- function(i) {
        sprintf("state%s", formatC(i, flag = "0", width = ceiling(log10(maxit+1))))
    }
    ## Plot intermediate results?
    assign("DO.PLOT", do.plot, globalenv())
  
    ## Setup database of results
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

    ## File for keeping track of 'theta' acceptances
    con <- file("theta.accept", "w")
    close(con)

    ## Initialize eta, gamma uniform grids
    if(covModel != "exp")
        stop("can only use 'exp' model for now")
    gammaVals <- seq(0, 1, length = ngridG)
    
    ## Pick a random initial point
    if(is.null(gamma0))
        gamma0 <- gammaVals[sample(ngridG, 2)]
       
    ## Initialize Gibbs state object
    gibbsState <- initGibbsState()
    params <- c("mu", "gamma")

    if(DO.PLOT) {
        if(length(dev.list()) > 0)
            graphics.off()
        x11(width = 5, height = 5)
        x11(width = 4, height = 4, xpos = -20) 
        x11(width = 5, height = 5, xpos = -450)
        x11(width = 5, height = 5, ypos = -100)
    }
    for(i in seq(maxit)) {
        gibbsState$iteration <- i
        
        if(verbose) cat("Iteration", i, "\n")

        for(j in seq(along = params)) {
            if(verbose) cat("\tUpdating", params[j], "\t")
            gibbsState <- updateGibbs(gibbsState, params[j])
            if(verbose) cat("\n")
        }
        dbInsert(db, makeKey(i), mget(params, gibbsState))
    }
    state <- if(exists(".Random.seed", globalenv())) 
        get(".Random.seed", globalenv())
    else 
        NULL
    .saveRDS(state, file = "random-seed.rds")
}

updateGibbs <- function(state, name) {
    state[[name]] <- switch(name,
                            mu = updateMu(state),
                            gamma = updateGamma2(state)
                            )
    state
}

updateTheta <- function(state) {
    library(mvtnorm)
    newtheta <- with(state, {
        OE <- dlCov(eta, d, sigmaE^2, covModel, rng = 0.6, shift = 0.2)        

        sampleTheta <- function(i) {
            ## LL <- logLik[[i]]
            LL <- profileLL[[i]]
            step <- delta[[i]]^2  ## Scaling multiplier for proposal
            theta0 <- theta[[i]]
            
            ## Proposal
            if(mle.proposal) {
                B <- OE %*% solve(vMLE[[i]] + OE)
                I <- diag(1, d)
                V <- (I - B) %*% OE
                theta.m <- mu + drop(B %*% (bMLE[[i]] - mu))
                thetastar <- mvrnorm(1, theta.m, step * V)
            }
            else {
                ## Simple random walk Metropolis
                thetastar <- mvrnorm(1, theta0, step * OE)
            }
            ## Metropolis (log) ratio
            num <- (LL(thetastar) + ldmvnorm(thetastar, mu, OE))
            denom <- (LL(theta0) + ldmvnorm(theta0, mu, OE))
            
            if(mle.proposal) {
                num <- num + ldmvnorm(theta0, theta.m, step * V)
                denom <- denom + ldmvnorm(thetastar, theta.m, step * V)
            }
            lr <- num - denom
            
            ## Accept/reject
            lu <- log(runif(1))
            accept.thetastar <- is.finite(lr) && !is.na(lr) && !is.nan(lr) && lu < lr

            cat(as.integer(accept.thetastar))
            
            if(accept.thetastar) 
                list(theta = thetastar, accept = accept.thetastar)
            else
                list(theta = theta0, accept = accept.thetastar)
        }
        out <- lapply(seq(n), sampleTheta)
        cat("\n")
        out
    })
    accept <- sapply(newtheta, "[[", "accept")

    ## Print to file
    cat(as.integer(accept), "\n", file = "theta.accept", append = TRUE)
    
    newtheta <- lapply(newtheta, "[[", "theta")

    cat("\n\t\tacceptance %:", round(100 * mean(accept), 1))

    if(DO.PLOT) {
        dev.set(dev.next())
        par(mar = c(4, 4, 1, 1))
        matplot(seq(0, length(newtheta[[1]]) - 1),
                    do.call("cbind", newtheta), type = "l",
                xlab = "Lag", ylab = expression(theta), ylim = c(-.02, .02))
    }
    newtheta
}

updateBeta <- function(state) {
    newbeta <- with(state, {
        zero <- rep(0, length(beta))
        
        sampleBeta <- function(i) {
            ## browser()
            D <- varbeta[[i]]
            LL <- logLik[[i]]

            ## Proposal (random walk Metropolis)
            beta0 <- beta[[i]]
            betastar <- mvrnorm(1, beta0, bdelta^2 * D)

            accept.betastar <- tryCatch({
                ## Metropolis (log) ratio
                num <- LL(theta[[i]], betastar) + ldmvnorm(betastar, zero, D)
                denom <- LL(theta[[i]], beta0) + ldmvnorm(beta0, zero, D)
                lr <- num - denom
                
                ## Accept/reject
                lu <- log(runif(1))
                is.finite(lr) && !is.na(lr) && !is.nan(lr) && lu < lr
            }, error = function(e) {
                FALSE
            })
            if(accept.betastar)
                list(beta = betastar, accept = accept.betastar)
            else
                list(beta = beta0, accept = accept.betastar)
        }
        out <- lapply(seq(n), sampleBeta)
        out
    })
    accept <- sapply(newbeta, "[[", "accept")
    newbeta <- lapply(newbeta, "[[", "beta")

    cat("\n\t\tacceptance %:", round(100 * mean(accept), 1))

    newbeta
}

updateMu <- function(state) {
    with(state, {
        I <- diag(1, d)
        OE <- dlCov(eta, d, sigmaE^2, covModel, rng = 0.6, shift = 0.2) / n
        OG <- dlCov(gamma, d, sigmaG^2, covModel, rng = 0.75, shift = 0.05)
        thetabar <- colMeans(do.call("rbind", theta))
        B <- OG %*% solve(OE + OG)

        V <- (I - B) %*% OG
        m <- drop(B %*% thetabar)  ## 0 + B %*% (thetabar - 0)

        if(DO.PLOT) {
            dev.set(dev.next())
            par(mar = c(4, 4, 1, 1))
            plot(seq(0, length(m) - 1), 1000 * m, type = "b",
                 ylim = c(-3, 5), pch = 20,
                 xlab = "Lag", ylab = expression(mu %*% 1000))
            abline(h = 0, lty = 3)
        }        
        drop(mvrnorm(1, m, V))
    })
}

updateGamma2 <- function(state) {
    with(state, {
        nr <- length(gammaVals)
        gamma.1 <- gamma

        if(DO.PLOT) {
            dev.set(dev.next())
            par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
        }
        for(i in 1:2) {
            lweights <- sapply(seq(nr), function(j) {
                ## browser()
                gamma.1[i] <- gammaVals[j]
                ldmvnorm(mu, 0, dlCov(gamma.1, d, sigmaG^2, covModel, 0.7, 0.05))
            })
            weights <- exp(lweights - max(lweights))
            weights <- weights / sum(weights)
            
            if(DO.PLOT) {
                plot(gammaVals, weights,
                     xlab = substitute(gamma[s], list(s=i)),
                     type = "h", ylab = "Conditional probability",
                     ylim = c(0, 1))
            }
            select <- gammaVals[sample(nr, 1, prob = weights)]
            gamma.1[i] <- select
        }
        cat("gamma:", gamma.1)
        gamma.1
    })

}


## Update eta one at a time; use common eta across counties
updateEta2 <- function(state) {
    with(state, {
        nr <- length(etaVals)
        eta1 <- eta
        thetam <- do.call("rbind", theta)

        ## browser()
        
        if(DO.PLOT) {
            dev.set(dev.next())
            par(mfrow = c(2, 1), mar = c(4,4,1,1))
        }
        for(k in 1:2) {  ## cycle over eta
            lweights <- sapply(seq(nr), function(i) {  ## cycle over etaVals
                eta1[k] <- etaVals[i]
                OE <- dlCov(eta1, d, sigmaE^2, covModel, 0.6, 0.2)
                L <- ldmvnorm(thetam, mu, OE)
                sum(L)
            })
            weights <- local({
                e <- exp(lweights - max(lweights))
                e / sum(e)
            })
            if(DO.PLOT) {
                plot(etaVals, weights,
                     xlab = substitute(eta[s], list(s=k)),
                     type = "h",
                     ylab = "Conditional probability",
                     ylim = c(0, 1)
                     )
            }
            eta1[k] <- sample(etaVals, 1, prob = weights)
        }
        cat("eta:", round(eta1, 4))
        eta1
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
## Old code

local({
    ## Update eta one at a time; use common eta across counties
    updateEta2 <- function(state, MetropolisEta = FALSE) {
        with(state, {
            nr <- length(etaVals)
            eta1 <- eta
            thetam <- do.call("rbind", theta)

            if(DO.PLOT) {
                dev.set(dev.next())
                par(mfrow = c(2, 1), mar = c(4,4,1,1))
            }
            a <- 4; b <- 1  ## For beta prior
            for(k in 1:2) {  ## cycle over eta
                if(MetropolisEta) {
                    
                    eta0 <- eta1
                    eta1[k] <- runif(1, min(etaVals), max(etaVals))
                    L1 <- (sum(ldmvnorm(thetam, mu,
                                        dlCov(eta1, d, sigmaE^2, covModel, 0.6, 0.2)))
                           + dbeta(eta0[k], a, b))
                    L2 <- (sum(ldmvnorm(thetam, mu,
                                        dlCov(eta0, d, sigmaE^2, covModel, 0.6, 0.2)))
                           + dbeta(eta1[k], a, b))
                    lr <- L1 - L2
                    cat(sprintf("eta*: %.2f (%d)\t", eta1[k], round(lr)))
                    lu <- log(runif(1))
                    accept <- if(!is.finite(lr) || is.na(lr) || is.nan(lr))
                        FALSE
                    else
                        lu < lr
                    if(!accept)
                        eta1[k] <- eta0[k]  ## Reset to the old value
                }
                else {
                    lweights <- sapply(seq(nr), function(i) {  ## cycle over etaVals
                        eta1[k] <- etaVals[i]
                        OE <- dlCov(eta1, d, sigmaE^2, covModel, 0.6, 0.2)
                        L <- ldmvnorm(thetam, mu, OE)
                        sum(L)
                    })
                    weights <- local({
                        e <- exp(lweights - max(lweights))
                        e / sum(e)
                    })
                    if(DO.PLOT) {
                        plot(etaVals, weights,
                             xlab = substitute(eta[s], list(s=k)),
                             type = "h",
                             ylab = "Conditional probability",
                             ylim = c(0, 1)
                             )
                    }
                    eta1[k] <- sample(etaVals, 1, prob = weights)
                }
            }
            cat("eta:", round(eta1, 4))
            eta1
        })
    }
})
