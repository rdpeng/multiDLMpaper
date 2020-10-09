######################################################################
## 2005-04-27
## Make tables with output from new error checking models

setwd("~/projects/MCAPS2")

## Configure this!
## We need to specify the location of the results files
## (i.e. results.p3.rda)

indir <- "raw-results"

######################################################################

outcomes <- dget("outcomeNames.R")

local({
    infile <- file.path(indir, paste("results", outcomes[1], "rda", sep = "."))
    load(infile)

    if(!exists("polls"))
        stop(sQuote("polls"), " object does not exist in ",
             sQuote(basename(infile)), " from which to get model names")
    assign("modelNames", names(polls), globalenv())
    assign("dfVals", dfVec, globalenv())
})

readOutcome <- function(outcome) {
    ## load 'results' and other objects
    infile <- file.path(indir, paste("results", outcome, "rda", sep = "."))
    load(infile)
    get("results")
}

procDFResult <- function(dfResult) {
    if(inherits(dfResult, "condition") || is.null(dfResult)) 
        data.frame(beta = NA, var = NA, condition = as.character(dfResult))
    else {
        r <- as.data.frame(dfResult)
        data.frame(r, condition = as.character(NA))
    }
}

procLagModel2 <- function(lagModel) {
    r <- lapply(lagModel, function(county) {
        data.frame(do.call("rbind", county), dfTime = dfVals)
    })
    countyNames <- names(lagModel)
    data.frame(do.call("rbind", r),
               fips = rep(countyNames, each = length(dfVals)))
}

procLagModel1 <- function(lagModel) {
    lapply(lagModel, function(county) {
        lapply(county, procDFResult)
    })
}

results <- lapply(outcomes, function(outcome) {
    cat(outcome, "\n")
    results <- readOutcome(outcome)
    
    ## Set to NA models that had errors/warnings or no estimates
    results <- lapply(results, procLagModel1)
    results <- lapply(results, procLagModel2)
    lagVar <- factor(rep(modelNames, sapply(results, nrow)))
    results <- data.frame(do.call("rbind", results), lag = lagVar)
})

outcomeVar <- factor(unlist(rep(names(outcomes), sapply(results, nrow))))

results <- data.frame(do.call("rbind", results), outcome = outcomeVar)
row.names(results) <- seq(nrow(results))

## Reorder columns
i <- match("condition", names(results))
results <- results[, c(seq(ncol(results))[-i], i)]


## Reorder factor levels
## results <- transform(results, lag = factor(as.character(lag), levels = c("Lag 0", "Lag 1", "Lag 2", "Lag -1", "Lag -2", "Mean 0:2")))
                              

######################################################################
## Write out table files


## Default single lag model
outfile <- file.path("results", "pm25.single-lag.rds")



## .saveRDS(results, file = outfile, compress = TRUE)

