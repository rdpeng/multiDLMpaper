## Full Poisson likelihood Gibbs sampler

## [2008-01-09] This file runs a modified model which doesn't have a
## second level to it (i.e. no 'eta' parameter)

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
                   singleAgeCat = TRUE) {
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

        ## File for keeping track of 'theta' acceptances
        file.create("mu.accept")

        params <- c("mu", "gamma")

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
                                mu = updateMu(state),
                                gamma = updateGamma3(state)
                                )
        state
}


updateMu <- function(state) {
        with(state, {
                zero <- rep(0, length(mu))
                OG <- dlCov(gamma, d, sigmaG^2, covModel, rng=0.75, shift=0.05)
                mustar <- mvrnorm(1, mu, 0.05^2 * OG)

                LL0 <- sapply(seq_len(n), function(i) {
                        LL <- profileLL[[i]]
                        LL(mu)
                })
                LLstar <- sapply(seq_len(n), function(i) {
                        LL <- profileLL[[i]]
                        LL(mustar)
                })
                sumLL0 <- sum(LL0) + ldmvnorm(mu, zero, OG)
                sumLLstar <- sum(LLstar) + ldmvnorm(mustar, zero, OG)
                lr <- sumLLstar - sumLL0
                lu <- log(runif(1))
                accept <- is.finite(lr) && !is.na(lr) && lu < lr

                cat(sprintf("mu: %d %.2f", as.integer(accept), lr))
                cat(as.integer(accept), "\n", file = "mu.accept", append = TRUE)
                
                if(accept)
                        mustar
                else
                        mu
        })
}

updateGamma3 <- function(state) {
        with(state, {
                gamma.0 <- gamma
                proposal <- numeric(2)
                
                for(k in 1:2) {
                        gamma.1 <- gamma.0

                        ## Proposal
                        prop <- runif(1, max(gamma.0[k] - 0.2, 0),
                                      min(gamma.0[k] + 0.2, 1))
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

ldmvnorm <- function(x, mean, sigma) {
        if(is.vector(x))
                x <- matrix(x, ncol = length(x))

        sigmaI <- solve(qr(sigma, LAPACK = TRUE))
        distval <- mahalanobis(x, center = mean, cov = sigmaI, inverted = TRUE)
        logdet <- as.numeric(determinant(sigma, log = TRUE)$modulus)
        logretval <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
        logretval
}



