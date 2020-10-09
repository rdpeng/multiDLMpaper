######################################################################
## 2005-04-05
## Make tables with output from new error checking models

setwd("~/projects/multiDLM")

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

## Start binding things together
procLagModel2 <- function(lagModel) {
    r <- lapply(lagModel, function(county) {
        data.frame(do.call("rbind", county),
                   dfTime = rep(dfVals, sapply(county, NROW)))
    })
    countyNames <- names(lagModel)
    data.frame(do.call("rbind", r), fips = rep(countyNames, sapply(r, NROW)))
}

## Process individual degrees of freedom within county, within lag
## model
procLagModel1 <- function(model) {
    ## Within a model, there are many counties
    r <- lapply(model, function(county) {
        ## Within a county, there are many degrees of freedom
        lapply(county, procDFResult)
    })
    use <- !sapply(r, inherits, what = "condition")

    if(any(!use))
        message("removing some counties with no results")
    r[use]
}

procDFResult <- function(dfResult) {
    ## Within a degree of freedom, there may be many estimates, if the
    ## model is a distributed lag model.  If we have DL model, we want
    ## the individual lag estimates, and the total effect (the sum)
    
    if(inherits(dfResult, "condition") || is.null(dfResult)) {
        d <- data.frame(beta = NA, var = NA,
                        condition = conditionMessage(dfResult))
        return(d)
    }  
    r <- with(dfResult, {
        data.frame(beta = serialize(beta, NULL, ascii = TRUE),
                   var = serialize(var, NULL, ascii = TRUE))
    })        
    data.frame(r, condition = as.character(NA))
}

results <- lapply(outcomes, function(outcome) {
    cat(outcome, "\n")
    results <- readOutcome(outcome)
    
    ## Set to NA models that had errors/warnings or no estimates
    results <- lapply(results, procLagModel1)
    results <- lapply(results, procLagModel2)
    modelVar <- factor(rep(modelNames, sapply(results, NROW)))
    results <- data.frame(do.call("rbind", results), model = modelVar)
    results
})

outcomeVar <- factor(rep(names(outcomes), sapply(results, NROW)))
results <- data.frame(do.call("rbind", results), outcome = outcomeVar)
row.names(results) <- seq(nrow(results))

## Reorder columns
i <- match("condition", names(results))
results <- results[, c(seq(ncol(results))[-i], i)]

## Convert beta/var factor to character; dfTime to integer
results$beta <- as.character(results$beta)
results$var <- as.character(results$var)
results$dfTime <- as.integer(results$dfTime)

attr(results, "runModels.Rout") <- readLines("runModels.Rout")

## .saveRDS(results, file = "results/pm25.dist-lag.rds", compress = TRUE)
## .saveRDS(results, file = "results/pm25.tsdecomp.rds", compress = TRUE)
## .saveRDS(results, file = "results/pm25.step.rds", compress = TRUE)
## .saveRDS(results, file = "results/pm25.orthoLag.rds", compress = TRUE)

######################################################################

expr <- expression({
    ## Only use counties in Medicare study
    medsites <- read.csv("../medicare/tables/fips-region.csv", colClasses = "factor")
    r <- subset(results, fips %in% medsites$fips)
    results <- transform(r, fips = factor(fips))

    sites <- levels(results$fips)

    ## Get some descriptives
    sinfo <- getDBSiteData()$siteInfo
    sitedata <- transform(subset(sinfo, fips %in% sites,
                                 c(fips, lat, long, pop100)),
                          lat = lat / 1e6, long = long / 1e6)

    ## Get avg denominators
    denom <- numeric(length(sites))
    names(denom) <- sites

    for(i in seq(along = sites)) {
        d <- collapse(readSite(sites[i], FALSE))
        denom[i] <- mean(d$COPDd, na.rm = TRUE)
        cat(sites[i], denom[i], "\n")
    }

    sitedata <- merge(sitedata, data.frame(fips = sites, denom = round(denom, 1)),
                      by = "fips")
})


