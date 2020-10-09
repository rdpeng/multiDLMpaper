## 'x' and 'v' are both lists

wmean <- function(x, v) {
    xm <- do.call("rbind", x)
    vinv <- lapply(v, solve)
    x1 <- lapply(seq(along = x), function(i) {
        crossprod(x[[i]], vinv[[i]])
    })
    va <- array(unlist(vinv), c(dim(vinv[[1]]), length(vinv)))
    svinv <- apply(va, c(1, 2), sum)
    x2 <- colSums(do.call("rbind", x1))
    V <- solve(svinv)

    vm <- do.call("rbind", lapply(v, diag))
    hetV <- pmax(apply(xm, 2, var) - colMeans(vm), 0)
    bh <- sapply(1:ncol(xm), function(i) {
        weighted.mean(xm[,i], 1 / (vm[,i] + hetV[i]))
    })
    vh <- sapply(1:ncol(vm), function(i) {
        1 / sum(1 / (vm[,i] + hetV[i]))
    })
    
    r <- list(b = drop(crossprod(x2, V)),
              v = V,
              bh = bh,
              vh = vh,
              sdh = sqrt(vh),
              std = sqrt(diag(V)))
    structure(r, class = "wmeanDistLag")
}

plot.wmeanDistLag <- function(x, ...) {
    xpts <- seq(along = x$b) - 1
    y <- x$b * 1000
    lo <- with(x, b - 2 * std) * 1000
    hi <- with(x, b + 2 * std) * 1000

    plot(xpts, y, ylim = range(lo, hi), type = "l")
    lines(xpts, lo, lty = 2)
    lines(xpts, hi, lty = 2)
}

itlnise <- function(b, v, seed = 12231) {
    library(tlnise)
    xm <- do.call("rbind", b)
    vm <- do.call("rbind", lapply(v, diag))

    initTLNise()
    g <- sapply(1:ncol(xm), function(i) {
        out <- tlnise(xm[,i], vm[,i], prnt = FALSE, seed = sample(1:1000, 1))
        out$gamma[, 1:2]
    })
    r <- list(b = g[1,],
              std = g[2,])
    structure(r, class = "wmeanDistLag")
}
