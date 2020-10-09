db2rds <- function(filename, n = 80000,
                   outfile = paste(filename, "rds", sep = ".")) {
        con <- gzfile(filename, "rb")
        on.exit(close(con))

        message("allocating data structures")
        mu <- matrix(nrow = n, ncol = 14)
        eta <- matrix(nrow = n, ncol = 2)
        gamma <- matrix(nrow = n, ncol = 2)
        theta.a <- array(dim = c(94, 14, n))

        message("reading database...")
        for(i in 1:n) {
                obj <- unserialize(con)
                mu[i, ] <- obj$mu
                eta[i, ] <- obj$eta
                gamma[i, ] <- obj$gamma
                theta.a[, , i] <- do.call("rbind", obj$theta)
        }
        r <- list(mu = mu, eta = eta, gamma = gamma)

        message("averaging thetas...")
        r$theta <- list(mean = apply(theta.a, c(1, 2), mean),
                        lo = apply(theta.a, c(1, 2), quantile, prob = 0.025),
                        hi = apply(theta.a, c(1, 2), quantile, prob = 0.975))
        message("saving output...")
        .saveRDS(r, file = outfile, compress = TRUE)
}
