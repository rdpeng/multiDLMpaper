################################################################################
## Copyright 2005, Roger D. Peng <rpeng@jhsph.edu>
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
################################################################################

################################################################################
## Load packages and data

## The 'tsModel' package is needed to create lags of variables or
## perhaps running means.

library(tsModel)
library(filehash)

sitedata <- .readRDS("data/all-sites-collapsed.rds")

## Post process glm object; extract various entities
postProcess <- function(glmObject) {
        ## Extract coefficients
        cc <- summary(glmObject)$coefficients
        
        ## Extract covariance matrix
        V <- vcov(glmObject)

        ## Only keep coefficients corresponding to pollutant variable
        i <- grep("pm25tmean", rownames(cc), fixed = TRUE)

        if(length(i) == 0) ## No match
                stop("pollutant coefficients not found in model fit")
        rval <- list(beta = cc[i, 1], var = V[i, i])

        if(is.null(rval) || length(rval) == 0)
                stop("problem with return value from 'fitSingleSite'")
        else
                rval
}

## The file 'modeling.R' contains the function 'fitSingleSite0()'
## which fits a county-specific model to each county.  By default,
## 'fitSingleSite()' is assigned to be 'fitSingleSite0()'.  

## Choose the fitting function

source("modeling.R")

fitSingleSite <- get("fitSingleSite0", mode = "function")

## The file 'runModels.R' calls 'fitSingleSite()' many times to fit
## models for each outcome, lag, and dfTime.

################################################################################

## Make results directory/DB

basedir <- "stage1"
dir.create(basedir, showWarnings = FALSE)

dbCreate(file.path(basedir, "stage1.db"), "DB1")
db <- dbInit(file.path(basedir, "stage1.db"), "DB1")

## Check for file 'siteList.R'.  If it doesn't exist, complain.
siteList <- if(file.exists("siteList.R")) {
        message("Using 'siteList.R'")
        dget("siteList.R")
} else {
        stop("'siteList' not found")
}
stopifnot(length(siteList) > 0)

## The vector of outcomes is set.  
outcomes <- dget("outcomeNames.R")

cat("Using outcomes:", paste(outcomes, collapse = " "), "\n")


## The vector of degrees of freedom for the smooth function of time is
## set (for sensitivity analysis).  Note these are the numbers of
## degrees of freedom *per year* of data (of which there are four
## years).

dfVec <- c(2, 4, 6, 8, 10, 12, 14)

## Expressions for the pollutant variables and lags.

polls <- c(
           "Lag(pm25tmean, 0)", 
           "Lag(pm25tmean, 1)", 
           "Lag(pm25tmean, 2)",
           "Lag(pm25tmean, 0:6)",
           "Lag(pm25tmean, 0:13)",
           "Lag(pm25tmean, 0:13) + pm25tmean:runMean(pm25tmean, 1:3)"
           )
names(polls) <- c("Lag 0", "Lag 1", "Lag 2", "DLag 0:6", "DLag 0:13",
                  "iDLag 0:13")

if(isTRUE(any(is.na(names(polls)))) || is.null(names(polls))) {
        stop("problem with 'polls'")
}

################################################################################
## Running the models:  The big loop

## As the loop goes through the counties, there may be errors in the
## model fitting.  Right now, these errors are caught in the
## 'tryCatch' block and are handled so that the execution of the
## entire loop is not interrupted.

models <- expand.grid(outcome = outcomes, poll = polls,
                      siteName = siteList, d0 = dfVec)

.saveRDS(models, file = file.path(basedir, "models.rds"))

overwrite <- FALSE

results <- lapply(seq_len(nrow(models)), function(i) {
        outcome <- as.character(models[i, "outcome"])
        poll <- as.character(models[i, "poll"])
        siteName <- as.character(models[i, "siteName"])
        d0 <- as.numeric(models[i, "d0"])
        d <- d0 * 4 ## There are 4 years of data

        key <- file.path(outcome, poll, siteName, d0)
        
        cat(outcome, poll, siteName, d0, "\t")
        
        if(dbExists(db, key) && !overwrite) {
                cat("[cached]\n")
                return(TRUE)
        }
        rval <- tryCatch({
                f <- fitSingleSite(sitedata[[siteName]], outcome = outcome,
                                   pollutant = poll, df.Time = d)
                postProcess(f)
        }, error = function(cond) {
                cat("\t", trunc(d / 4), as.character(cond))
                cond$call <- NULL
                cond
        }, warning = function(cond) {
                cat("\t", trunc(d / 4), as.character(cond))
                cond$call <- NULL
                cond
        })
        dbInsert(db, key, rval)
        
        cat("\n")

        !inherits(rval, "condition")
})
                  
        
## .saveRDS(results, file = "modelResults.rds")
