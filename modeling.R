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

fitSingleSite0 <- function(data,  ## A data frame for a site
                           outcome,  ## character, name of outcome
                           pollutant,  ## character, name of pollutant variable
                           denom,  ## character, denominator variable name
                           df.Time = 6 * 4,  ## df for smooth function of time
                           df.time = 1 * 4,  ## df for smooth function x age category
                           df.Temp = 6,  ## df for temp smooth function
                           df.Dew = 3,
                           subset = TRUE, singleAgeCat = TRUE,
                           fitModel = TRUE, ...) {  
    library(splines)

    stopifnot(is.character(outcome), is.character(pollutant))
    if(missing(denom)) {
        denom <- if(length(grep("admit", outcome, fixed = TRUE) > 0))
            sub("admit", "denom", outcome, fixed = TRUE)
        else
            paste(outcome, "d", sep = "")
    }        
    form <- reformulate(c("dow",
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
                          paste("offset(log(", denom, "))"),
                          pollutant),
                        response = outcome)
    if(fitModel)
        glm(form, data = data, family = quasipoisson, na.action = na.exclude,
            control = glm.control(epsilon = 1e-10, maxit = 1000), subset = subset)
    else
        form
}


createWeekend <- function(dataframe) {
    stopifnot("date" %in% names(dataframe))
    weekend <- with(dataframe, {
        factor(as.integer(weekdays(date) %in% c("Saturday", "Sunday")),
               labels = c("No", "Yes"))
    })
}

selectMonitorDailyValue <- function(x, y, z) {
    ## x:  closest monitor
    ## y:  25 km away
    ## z:  35 km away
    m <- cbind(x, y, z)
    stopifnot(is.matrix(m))
    v <- apply(m, 1, function(r) {
        if(all(is.na(r)))
            NA
        else if(all(!is.na(r)))
            r[1]
        else {
            rn <- r[-which(is.na(r))]
            rn[1]
        }
    })
    stopifnot(length(v) == nrow(m))
    v
}
