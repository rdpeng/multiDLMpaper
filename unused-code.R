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


