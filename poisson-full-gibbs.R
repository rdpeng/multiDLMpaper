## Full Poisson likelihood Gibbs sampler

######################################################################
## Construct constraint covariance matrix

## Use "exp" or "exp2"
## For gamma parameters we use rng = 0.7 and shift = 0.05
## For eta parameters we use rng = 0.6 and shift = 0.2

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


## Construct a function that evaluates the log-likelihood for theta
## and beta (nuisance parameters)

makeLL <- function(data, outcome, pollutant) {
    data <- as(data, "data.frame")
    ## data <- collapse(data)  ## Collapse age categories for now
    fit0 <- fitSingleSite(data, outcome, pollutant, singleAgeCat = TRUE)

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

    ## browser()

    ## Vprior <- diag(1000, length(beta0))
    ## M <- Vprior %*% solve(var0 + Vprior)
    ## I <- diag(1, length(beta0))
    ## beta0 <- M %*% beta0
    ## var0 <- (I - M) %*% Vprior

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


pfgibbs <- function(sitedata,
                    outcome = "COPDAE",
                    pollutant = "Lag(pm25tmean, 0:13)",
                    bMLE = NULL,
                    vMLE = NULL,
                    mle.proposal = TRUE,
                    maxit = 80000,
                    sigmaE = 0.005,
                    sigmaG = 0.005,
                    covModel = "exp",
                    d = 14,  ## dimension of theta
                    seed = 100,
                    delta = 1,
                    bdelta = 1,
                    eta0 = NULL,
                    gamma0 = NULL,
                    verbose = TRUE,
                    dbfile = "statepfgibbs", 
                    deleteCache = TRUE,
                    singleAgeCat = TRUE) {
    library(MASS)
    library(filehash)

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
        s$eta <- eta0
        s$gamma <- gamma0
        s$sigmaE <- sigmaE
        s$sigmaG <- sigmaG

        s$sitedata <- sitedata
        s$outcome <- outcome
        s$pollutant <- pollutant
        
        s$covModel <- covModel
        s$n <- n  ## number of sites
        s$d <- d  ## dimension

        out <- if(file.exists("cacheLogLik.rds")) {
             if(verbose) cat("using cached log likelihoods\n")
            .readRDS("cacheLogLik.rds")
        }
        else {
            if(verbose) cat("creating log likelihood functions\n")
            pLL <- lapply(seq(along = sitedata), function(i) {
                if(verbose) cat(i, " ")
                data <- sitedata[[i]]
                makeLL(data, outcome, pollutant)
            })
            if(verbose) cat("\n")
            if(verbose) cat("caching log likelihoods for later use\n")
            .saveRDS(pLL, file = "cacheLogLik.rds", compress = TRUE)
           pLL
        }
        s$logLik <- lapply(out, "[[", "LL")
        s$beta <- lapply(out, "[[", "beta0")
        s$betaInit <- s$beta
        s$varbetaInit <- lapply(out, "[[", "var0")
        ## s$theta <- lapply(out, "[[", "theta0")
        s$theta <- replicate(n, rep(0, d), simplify = FALSE)
        
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

    ## Pick a random initial point
    if(is.null(eta0))
        eta0 <- runif(2)
    if(is.null(gamma0))
        gamma0 <- runif(2)
       
    ## Initialize Gibbs state object
    gibbsState <- initGibbsState()
    params <- c("theta", "beta", "mu", "gamma", "eta")
    ## params <- c("theta", "mu", "gamma", "eta")

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
}

updateGibbs <- function(state, name) {
    state[[name]] <- switch(name,
                            theta = updateTheta(state),
                            beta = updateBeta(state),
                            mu = updateMu(state),
                            eta = updateEta3(state),
                            gamma = updateGamma3(state)
                            )
    state
}

## This version uses a profile likelihood for theta

updateTheta <- function(state) {
    library(mvtnorm)
    newtheta <- with(state, {
        OE <- dlCov(eta, d, sigmaE^2, covModel, rng = 0.6, shift = 0.2)        

        sampleTheta <- function(i) {
            LL <- logLik[[i]]
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
            num <- (LL(thetastar, beta[[i]]) + ldmvnorm(thetastar, mu, OE))
            denom <- (LL(theta0, beta[[i]]) + ldmvnorm(theta0, mu, OE))
            
            if(mle.proposal) {
                num <- num + ldmvnorm(theta0, theta.m, step * V)
                denom <- denom + ldmvnorm(thetastar, theta.m, step * V)
            }
            lr <- num - denom
            
            ## Accept/reject
            lu <- log(runif(1))
            accept.thetastar <- is.finite(lr) && !is.na(lr) && lu < lr

            ## cat(as.integer(accept.thetastar))
            
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

    newtheta
}

updateBeta <- function(state) {
    newbeta <- with(state, {
        sampleBeta <- function(i) {
            ## browser()

            ## Prior variance for beta
            Vbeta <- diag(100, length(beta[[i]]))
            zero <- rep(0, length(beta[[i]]))
            D <- varbetaInit[[i]]
            ## M <- D %*% solve(D + Vbeta)
            ## I <- diag(1, length(beta[[i]]))
           
            ## Proposal 
            beta0 <- beta[[i]]
            betastar <- mvrnorm(1, beta0, bdelta^2 * D)
            ## pMu <- drop(M %*% betaInit[[i]])
            ## pSigma <- bdelta^2 * (I-M) %*% D
            ## betastar <- mvrnorm(1, pMu, pSigma)

            LL <- logLik[[i]]

            accept.betastar <- tryCatch({
                ## Metropolis (log) ratio
                num <- LL(theta[[i]], betastar) + ldmvnorm(betastar, zero, Vbeta)
                ## num <- num + ldmvnorm(beta0, pMu, pSigma)
                denom <- LL(theta[[i]], beta0) + ldmvnorm(beta0, zero, Vbeta)
                ## denom <- denom + ldmvnorm(betastar, pMu, pSigma)
                lr <- num - denom
                
                ## Accept/reject
                lu <- log(runif(1))
                is.finite(lr) && !is.na(lr) && !is.nan(lr) && lu < lr
            }, error = function(e) {
                cat(as.character(e), "\n")
                FALSE
            })

            ## cat(as.integer(accept.betastar))
            
            if(accept.betastar)
                list(beta = betastar, accept = accept.betastar)
            else
                list(beta = beta0, accept = accept.betastar)
        }
        out <- lapply(seq(n), sampleBeta)
        cat("\n")
        out
    })
    accept <- sapply(newbeta, "[[", "accept")
    newbeta <- lapply(newbeta, "[[", "beta")

    cat(as.integer(accept), "\n", file = "beta.accept", append = TRUE)

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

        drop(mvrnorm(1, m, V))
    })
}

updateGamma3 <- function(state) {
    with(state, {
        gamma.0 <- gamma
        proposal <- numeric(2)
        
        for(k in 1:2) {
            gamma.1 <- gamma.0

            ## Proposal
            prop <- runif(1, max(gamma.0[k] - 0.2, 0), min(gamma.0[k] + 0.2, 1))
            gamma.1[k] <- proposal[k] <- prop
            ## gamma.1[k] <- runif(1)

            ## Metropolis (log) ratio
            num <- ldmvnorm(mu, 0, dlCov(gamma.1, d, sigmaG^2, covModel, 0.7, 0.05))
            denom <- ldmvnorm(mu, 0, dlCov(gamma.0,d,sigmaG^2,covModel,0.7, 0.05))

            lr <- num - denom
            lu <- log(runif(1))
            accept <- is.finite(lr) && !is.na(lr) && lu < lr

            if(accept) 
                gamma.0[k] <- gamma.1[k]
        }
        cat(sprintf("gamma: %.4f %.4f (%.4f %.4f)\n",
                    gamma.0[1], gamma.0[2], proposal[1], proposal[2]))
        ## cat("gamma:", round(gamma.0, 4))
        gamma.0
    })
}

## Update eta one at a time; use common eta across counties
updateEta3 <- function(state) {
    with(state, {
        eta.0 <- eta
        thetam <- do.call("rbind", theta)
        proposal <- numeric(2)
        
        for(k in 1:2) {  ## cycle over eta
            eta.1 <- eta.0

            ## Proposal
            ## shift <- 0.05
            shift <- 0.025  ## For CVD only
            prop <- runif(1, max(eta.0[k] - shift, 0), min(eta.0[k] + shift, 1))
            eta.1[k] <- proposal[k] <- prop

            ## Metropolis (log) ratio
            OE1 <- dlCov(eta.1, d, sigmaE^2, covModel, 0.6, 0.2)
            num <- sum( ldmvnorm(thetam, mu, OE1) )
            OE0 <- dlCov(eta.0, d, sigmaE^2, covModel, 0.6, 0.2)
            denom <- sum( ldmvnorm(thetam, mu, OE0) )

            lr <- num - denom
            lu <- log(runif(1))
            accept <- is.finite(lr) && !is.na(lr) && lu < lr

            if(accept)
                eta.0[k] <- eta.1[k]
        }
        cat(sprintf("eta: %.4f %.4f (%.4f %.4f)\n",
                    eta.0[1], eta.0[2], proposal[1], proposal[2]))
        
        ## cat("eta:", round(eta.0, 4))
        eta.0
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

