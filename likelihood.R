################################################################################
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
                          subset = TRUE, singleAgeCat = TRUE,
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
            control = glm.control(epsilon = 1e-10, maxit = 1000),
            subset = subset)
}


## Construct a function that evaluates the profile likelihood for
## theta.

makeProfLLik <- function(data, outcome, pollutant, singleAgeCat = FALSE) {
        ## data <- as(data, "data.frame")
        fit0 <- fitSingleSite(data, outcome, pollutant,
                              singleAgeCat = singleAgeCat)
        
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

        ## part of model matrix corresponding to PM25 
        theta.mm <- mm[, theta.idx]
        ## part of model matrix corresponding to other variables
        beta.mm <- mm[, -theta.idx]

        rm(fit0, data, mm, mf)
        
        function(theta) {
                ## Hold theta fixed and maximize over beta (nuisance parameters)
                g <- glm.fit(beta.mm, y, offset = off + theta.mm %*% theta,
                             family = poisson())
                with(g, rank - aic / 2)
        }
}


pLL <- makeProfLLik(d, "COPDAE", "Lag(pm25tmean, 0:13)", TRUE)

f <- fitSingleSite(d,"COPDAE", "Lag(pm25tmean,0:13)")
b <- coef(f)[50:63]
v <- vcov(f)[names(b), names(b)]
b <- unname(b)

ldmvnorm <- function(x, mean, sigma) {
        if(is.vector(x))
                x <- matrix(x, ncol = length(x))

        sigmaI <- solve(qr(sigma, LAPACK = TRUE))
        distval <- mahalanobis(x, center = mean, cov = sigmaI, inverted = TRUE)
        logdet <- as.numeric(determinant(sigma, log = TRUE)$modulus)
        logretval <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
        logretval
}

makeApproxPLL <- function(b, v) {
        function(theta) {
                ldmvnorm(theta, b, v)
        }
}


apLL <- makeApproxPLL(b, v)
