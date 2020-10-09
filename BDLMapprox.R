###############################################################################
## Copyright (C) 2005, 2006, Leah J. Welty <lwelty@northwestern.edu>
## with modifications by Roger D. Peng <rpeng@jhsph.edu
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
###############################################################################

BDLMeta <- function(b, v, eta, sigmaE = 0.005) {
    d <- length(b)
    OE <- dlCov(eta, d, sigmaE^2)
    B <- OE %*% solve(v + OE)
    V <- B %*% v
    list(mu = drop(B %*% b), V = V, std = sqrt(diag(V)))
}

BDLMglm <- function(x, pollutant, sigma = 0.004, eta.mat = NULL, ...) {
    cfs <- summary(x)$coefficients
    i <- grep(pollutant, rownames(cfs), fixed = TRUE)
    theta.hat <- cfs[i, "Estimate"]
    V <- vcov(x)[i, i]
    BDLM(theta.hat, V, sigma, eta.mat)
}

BDLM <- function(x,  ## MLEs from Poisson model
                 cov.hat,  ## Estimated covariance matrix of theta
                 sigma = 0.005,  ## Prior variance
                 eta.mat = NULL) {
    theta.hat <- x
    hyper.mat <- if (is.null(eta.mat)) {
        as.matrix(expand.grid(seq(-0.35, -0.05, length = 10),
                              seq(-0.37, 0, length = 10)))
    }
    else if(!is.matrix(eta.mat) && length(eta.mat) == 2)
        matrix(eta.mat, byrow = TRUE, ncol = 2)
    else {
        stopifnot(ncol(eta.mat) == 2, is.matrix(eta.mat))
        eta.mat
    }
    n.hyperpar <- nrow(hyper.mat)
    nLags <- length(theta.hat)
    inv.cov.hat <- solve(cov.hat)
    
    ## initialize posterior mean and covariance vectors
    val <- vector(length = n.hyperpar)
    postmeans <- matrix(nrow = n.hyperpar, ncol = nLags)
    postvar <- array(dim = c(nLags, nLags, n.hyperpar))
    
    for (k in 1:n.hyperpar) {
        param <- as.numeric(hyper.mat[k,])

        ## Construct prior covariance matrix
        pr.cov <- dlCov(param = param, n = nLags, sigmasq = sigma^2)

        ## compute inv of sum of inverses s.t. don't invert pr.cov
        inv.mat <- solve(pr.cov %*% inv.cov.hat + diag(1, nLags))
        inv.sum.cov.inv <- inv.mat %*% pr.cov
        val[k] <- 0.5 * determinant(inv.mat, log = TRUE)$modulus -  
            0.5 * crossprod(theta.hat, inv.cov.hat - inv.cov.hat %*% inv.sum.cov.inv %*% inv.cov.hat) %*% theta.hat
        postmeans[k, ] <- inv.sum.cov.inv %*% inv.cov.hat %*% theta.hat
        postvar[, , k] <- inv.sum.cov.inv
    }
    bfac <- vector(length = n.hyperpar)  # weights

    for (k in 1:n.hyperpar) 
        bfac[k] <- sum(exp(val - val[k]))^(-1)

    av.theta <- colSums(postmeans * bfac)
    av.cov <- matrix(0, ncol = nLags, nrow = nLags)

    for (k in 1:n.hyperpar) 
        av.cov <- av.cov + postvar[, , k] * bfac[k]^2

    structure(list(post.means = postmeans, post.var = postvar, bfac = bfac,
                   av.theta = av.theta, av.cov = av.cov),
              class = "BDLMapprox")
}

plot.BDLMapprox <- function(x, y, m1000 = TRUE, ...) {
    n <- length(x$av.theta)
    mult <- if(m1000)
        1000
    else
        1
    y <- x$av.theta * mult
    std <- sqrt(diag(x$av.cov)) * mult
    xpts <- seq(0, n - 1)
    rng <- range(y - 2*std, y + 2*std)
    plot(xpts, y, type = "b", ylim = rng, xlab = "Lag",
         ylab = expression(hat(beta)))
    lines(xpts, y + 2*std, lty = 2)
    lines(xpts, y - 2*std, lty = 2)

    invisible()
}


######################################################################
## Support functions

dlCov <- function(param, n, sigmasq) {
    varFun <- function(lambda, u) {
        ## variance of distributed lag u
        exp(lambda * u)
        ## (u + 1)^lambda
    }

    corFun <- function(gamma, u) {
        ## weight function for correlation matrix
        exp(gamma * u)
    }

    genCorMat <- function(gamma, n) {
        ## correlation matrix for distributed lag vector
        if(!isTRUE(n > 1))
            stop("argument 'n' should be > 1")
        
        D <- diag(corFun(gamma, 0:(n - 1)))
        M <- matrix(1, nrow = n, ncol = n)
        I <- diag(1, n)
        cov <- D %*% t(D) + (I-D) %*% M %*% t(I-D)
        cov2cor(cov)
    }

    ## prior covariance matrix for n x 1 distributed lag vector
    ## param is 2 parameters determining covariance of distributed lag vector
    eta1 <- param[1]
    eta2 <- param[2]
    R <- genCorMat(eta2, n)
    V <- diag(sqrt(varFun(eta1, 0:(n - 1))))
    sigmasq * (V %*% R %*% V)
}
