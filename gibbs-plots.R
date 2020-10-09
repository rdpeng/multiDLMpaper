######################################################################
## Results methods

readDB <- function(file) {
        maxit <- 100000
        results <- vector("list", length = maxit)
        con <- gzfile(file, "rb")
        on.exit(close(con))
        status <- TRUE

        for(i in seq_len(maxit)) {
                status <- try({
                        results[[i]] <- unserialize(con)
                }, silent = TRUE)
                if(inherits(status, "try-error"))
                        break
        }
        results[seq_len(i - 1)]
}

convertDB <- function(file) {
        db <- dbInit(file, "DB1")
        message("reading keys")
        keys <- sort(dbList(db))
        con <- gzfile(paste(file, "gz", sep = "."), "wb")
        on.exit(close(con))

        message("rewriting database")
        for(i in seq_along(keys)) {
                obj <- dbFetch(db, keys[i])
                serialize(obj, con)
        }
}

readEW <- function(file, n = 100000) {
        con <- gzfile(file, "rb")
        on.exit(close(con))
        
        EW <- array(dim = c(2, 14, n))
        siteInfo <- dget("siteInfo.R")
        east <- siteInfo$region == "East"
        
        for(i in seq_len(n)) {
                status <- try({
                        out <- unserialize(con)
                        theta <- do.call("rbind", out$theta)
                        EW[1, , i] <- colMeans(theta[east, ])
                        EW[2, , i] <- colMeans(theta[!east, ])
                }, silent = TRUE)
                if(inherits(status, "try-error"))
                        break
        }
        EW[, , seq_len(i - 1)]
}

readNS <- function(file, n = 100000) {
        con <- gzfile(file, "rb")
        on.exit(close(con))
        
        NS <- array(dim = c(2, 14, n))
        siteInfo <- dget("siteInfo.R")
        north <- siteInfo$lat > 36.5
        
        for(i in seq_len(n)) {
                status <- try({
                        out <- unserialize(con)
                        theta <- do.call("rbind", out$theta)
                        NS[1, , i] <- colMeans(theta[north, ])
                        NS[2, , i] <- colMeans(theta[!north, ])
                }, silent = TRUE)
                if(inherits(status, "try-error"))
                        break
        }
        NS[, , seq_len(i - 1)]
}


compareRegion <- function(ew, ...) {
        
}

plotRdiff <- function(ew, ...) {
        m <- apply(ew, c(1, 2), mean)
        qq <- apply(ew[1, , ] - ew[2, , ], 1, quantile, prob = c(0.025, 0.975))
        m <- m * 1000
        qq <- qq * 1000
        plot(0:13, m[1, ] - m[2, ], type = "o", ylim = range(qq), ...)
        lines(0:13, qq[1,], lty = 3)
        lines(0:13, qq[2,], lty = 3)
        abline(h = 0, lty = 2)
}

readMu <- function(db){
        library(filehash)
        if(is.character(db)) 
                db <- dbInit(db)
        keys <- sort(dbList(db))
        unname(lapply(keys, function(k) {
                with(dbFetch(db, k), {
                        list(mu = mu, eta = eta, gamma = gamma)
                })
        }))
}

plotGibbs <- function(x, m1000 = TRUE, burnin = 1000, ylim = range(lo, hi),
                      ...) {
        f <- ifelse(m1000, 1000, 1)
        m <- if(!is.matrix(x))
                do.call("rbind", lapply(x, "[[", "mu"))
        else
                x
        m <- if(nrow(m) < (burnin + 100))
                m
        else
                m[-seq(burnin), ]
        mu <- colMeans(m) * f
        xpts <- seq(0, ncol(m) - 1)
        hi <- apply(m, 2, quantile, prob = 0.975) * f
        lo <- apply(m, 2, quantile, prob = 0.025) * f

        plot(xpts, mu, type = "b", pch = 20, ylim = ylim, ,
             main = paste("Chain:", nrow(m), "iterations"), xlab = "Lag (days)",
             ylab = expression(hat(theta) %*% 1000), ...)
        abline(h = 0, lty = 3)
        lines(xpts, lo, lty = 2)
        lines(xpts, hi, lty = 2)
}

total <- function(g) {
        sapply(g, function(x) sum(x$mu))
}

marginals <- function(g) {
        list(mu = do.call("rbind", lapply(g, "[[", "mu")),
             gamma = sapply(g, "[[", "gamma"))
}

histGibbsGamma <- function(g) {
        ga <- sapply(g, "[[", "gamma")
        hist(ga, xlab = expression(gamma),
             ylab = "Posterior distribution")
}

plotGibbsGamma <- function(g) {
        library(KernSmooth)
        library(RColorBrewer)
        pal <- colorRampPalette(brewer.pal(4, "Blues"))
        gamma <- do.call("rbind", lapply(g, "[[", "gamma"))
        k <- bkde2D(gamma, 0.03)
        with(k, image(x1, x2, fhat, col = pal(100)))
        box()
}

plotGibbsEta <- function(g) {
        library(KernSmooth)
        library(RColorBrewer)
        pal <- colorRampPalette(brewer.pal(4, "Blues"))
        eta <- do.call("rbind", lapply(g, "[[", "eta"))
        k <- bkde2D(eta, 0.03)
        with(k, image(x1, x2, fhat, col = pal(100)))
        box()
}

plotEta <- function(g) {
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        par(mfrow = c(2, 1), mar = c(2, 4, 1, 1))

        eta <- t(sapply(g, "[[", "eta"))
        plot(eta[,1], type = "l", xlab = "", ylab = expression(eta[1]))
        plot(eta[,2], type = "l", xlab = "", ylab = expression(eta[2]))
}

plotGamma <- function(g) {
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        par(mfrow = c(2, 1), mar = c(2, 4, 1, 1))

        gamma <- t(sapply(g, "[[", "gamma"))
        plot(gamma[,1], type = "l", xlab = "", ylab = expression(gamma[1]))
        plot(gamma[,2], type = "l", xlab = "", ylab = expression(gamma[2]))
}

trackTheta <- function(file = "theta.accept", nc = 94) {
        m <- matrix(scan(file, integer(0)), byrow = TRUE, ncol = nc)
        print(round(100 * colMeans(m)))
        invisible(m)
}
