if(!file.exists("etaFix.rda")) {
    etaFix <- local({
        source("loadResults.R", local = TRUE)
        source("BDLMapprox.R", local = TRUE)

        cat("Loading results...\n")
        loadResults("COPDAE", 1000)
        cat("Applying BDLM...\n")
        bdlm <- lapply(seq(along = b), function(i) BDLM(b[[i]], v[[i]]))

        cat("Averaging eta's...\n")
        etaMat <- as.matrix(expand.grid(seq(-0.35, -0.05, length = 10),
                                        seq(-0.37, 0, length = 10)))
        
        eta <- lapply(bdlm, function(x) colSums(etaMat * x$bfac))
        names(eta) <- names(b)
        eta
    })
    save(etaFix, file = "etaFix.rda")
} else {
    cat("Loading 'etaFix.rda'\n")
    load("etaFix.rda", globalenv())
}
