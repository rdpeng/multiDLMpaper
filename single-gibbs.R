######################################################################
## County-specific time series model

fitSingleSite <- function(data,  ## A data frame for a site
                          outcome,  ## character, name of outcome
                          pollutant,  ## character, name of pollutant variable
                          theta = NULL, ## pollutant coefficients
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
    fips <- as.character(data$fips[1])

    if(fips == "24005") {
        data <- transform(data, tmpd = (tmax35 + tmin35) / 2, dptp = dptp35)
        data <- transform(data, rmtmpd = runMean(tmpd, 1:3),
                          rmdptp = runMean(dptp, 1:3))
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
                    sprintf("offset(%s * %s)", pollutant, theta)
                }
                )
    form <- reformulate(xterms, response = outcome)
    glm(form, data = data, family = poisson, na.action = na.omit,
        control = glm.control(epsilon = 1e-10, maxit = 1000), subset = subset)
}

## Make function that evaluates the profile log-likelihood for theta.
## Actually, this function computes an estimated likelihood.

makeProfLLik <- function(data, outcome, pollutant) {
    data <- as(data, "data.frame")
    fit0 <- fitSingleSite(data, outcome, pollutant)
    LL0 <- as.numeric(logLik(fit0))

    function(theta) {
        fit <- fitSingleSite(data, outcome, pollutant, theta)
        LL <- logLik(fit)
        as.numeric(LL) - LL0
    }
}

makeProfNLik <- function(mu, sigma) {
    function(theta) {
        dnorm(theta, mu, sigma, log = TRUE) - dnorm(mu, mu, sigma, log = TRUE)
    }
}

makeEstLLik <- function(data, outcome, pollutant) {
    data <- as(data, "data.frame")
    fit0 <- fitSingleSite(data, outcome, pollutant)
    LL0 <- logLik(fit0)

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
    lagMat <- sweep(mm[, idx.lag, drop = FALSE], 2,
                    colMeans(mm[, idx.lag, drop = FALSE]))
    y <- model.response(mf)

    ## Remove big items from the environment
    rm(fit0, data, mm, mf, p)
    
    function(theta) {
        link <- drop(lagMat %*% theta) + others
        L <- dpois(y, exp(link), log = TRUE)
        sum(L) - as.numeric(LL0)
    }
}

psinglegibbs <- function() {

}
