## Full Poisson likelihood Gibbs sampler

################################################################################
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
## Main 'pgibbs()' function


pgibbs <- function(gibbsState,
                   maxit = 80000,                   
                   verbose = TRUE,
                   dbfile = "statepgibbs", 
                   deleteCache = FALSE,
                   singleAgeCat = TRUE,
                   sigmaE = NULL,
                   delta = NULL) {
        library(MASS)

        ## Setup database of results
        if(file.exists(dbfile)) {
                if(deleteCache) {
                        message("removing existing cache file")
                        file.remove(dbfile)
                }
                else
                        stop(sprintf("cache file '%s' already exists", dbfile))
        }
        con <- gzfile(dbfile, "wb")
        on.exit(close(con))

        if(!is.null(sigmaE))
                gibbsState$sigmaE <- sigmaE
        if(!is.null(delta))
                gibbsState$delta <- delta

        ## File for keeping track of 'theta' acceptances
        file.create("theta.accept")

        params <- c("theta", "mu", "gamma", "eta")

        for(i in seq_len(maxit)) {
                if(file.exists("KILL_PGIBBS")) {
                        message("saving full 'gibbsState' object for restarting")
                        serialize(gibbsState, con)
                        message("saving '.Random.seed'")
                        dput(.Random.seed, file = "RANDOM_SEED.R")
                        break
                }
                if(verbose) cat("Iteration", i, "\n")

                for(j in seq_along(params)) {
                        if(verbose) cat("\tUpdating", params[j], "\t")
                        gibbsState <- updateGibbs(gibbsState, params[j])
                        if(verbose) cat("\n")
                }
                out <- gibbsState[params]
                serialize(out, con)
                flush(con)
        }
        message(sprintf("exiting 'pgibbs' at %d iterations", i - 1))
}

updateGibbs <- function(state, name) {
        state[[name]] <- switch(name,
                                theta = updateTheta(state),
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
                OE <- dlCov(eta, d, sigmaE^2, covModel, rng = 2, shift = 0.2)

                out <- lapply(seq_len(n), function(i) {
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
                                denom <- denom + ldmvnorm(thetastar,theta.m,step*V)
                        }
                        lr <- num - denom
                        
                        ## Accept/reject
                        lu <- log(runif(1))
                        accept.thetastar <- is.finite(lr) && !is.na(lr) && lu < lr

                        cat(as.integer(accept.thetastar))
                        
                        if(accept.thetastar) 
                                list(theta = thetastar, accept = accept.thetastar)
                        else
                                list(theta = theta0, accept = accept.thetastar)
                })
                cat("\n")
                out
        })
        accept <- sapply(newtheta, "[[", "accept")

        ## Print to file
        cat(as.integer(accept), "\n", file = "theta.accept", append = TRUE)
        cat("\n\t\tacceptance %:", round(100 * mean(accept), 1))
        
        newtheta <- lapply(newtheta, "[[", "theta")
        newtheta
}

updateMu <- function(state) {
        with(state, {
                I <- diag(1, d)
                OE <- dlCov(eta, d, sigmaE^2, covModel, rng = 2, shift = 0.2) / n
                OG <- dlCov(gamma, d, sigmaG^2, covModel, rng = 1, shift = 0.05)
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
                        ## shift <- 0.2
                        shift <- 0.05
                        prop <- runif(1, max(gamma.0[k] - shift, 0),
                                      min(gamma.0[k] + shift, 1))
                        gamma.1[k] <- proposal[k] <- prop

                        ## Metropolis (log) ratio
                        num <- ldmvnorm(mu, 0, dlCov(gamma.1, d, sigmaG^2,
                                                     covModel, 0.7, 0.05))
                        denom <- ldmvnorm(mu, 0, dlCov(gamma.0, d, sigmaG^2,
                                                       covModel, 0.7, 0.05))

                        lr <- num - denom
                        lu <- log(runif(1))
                        accept <- is.finite(lr) && !is.na(lr) && lu < lr

                        if(accept) 
                                gamma.0[k] <- gamma.1[k]
                }
                cat(sprintf("gamma: %.4f %.4f (%.4f %.4f)\n",
                            gamma.0[1], gamma.0[2], proposal[1], proposal[2]))
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
                        ## shift <- 0.01
                        shift <- 0.025
                        prop <- runif(1, max(eta.0[k] - shift, 0),
                                      min(eta.0[k] + shift, 1))
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



