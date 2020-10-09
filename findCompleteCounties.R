######################################################################
## Find counties with everyday data
######################################################################


library(MCAPSdb2)

allsites <- listSites()
keep <- logical(length(allsites))
names(keep) <- allsites

for(i in seq(along = allsites)) {
        data <- readSite(allsites[i], asDataFrame = FALSE)
        metpoll <- data$exposure
        pm25 <- metpoll$pm25tmean
        years <- as.POSIXlt(metpoll$date)$year + 1900
        
        nmiss <- tapply(pm25, years, function(x) sum(is.na(x)))
        cat(allsites[i], ":", nmiss, " ")
        
        keep[i] <- sort(nmiss)[2] < 90

        cat(keep[i], "\n")
}

## These are the sites we'll use
sites <- allsites[keep]

dput(sites, file = "siteList.R")

## Get some descriptives
sinfo <- getMetaData("siteInfo")
sitedata <- transform(subset(sinfo, fips %in% sites,
                             c(fips, lat, long, pop100)),
                      lat = lat / 1e6, long = long / 1e6)

## Get avg covariates
vars <- lapply(sites, function(site) {
        cat(site, "\n")
        d <- collapse(readSite(site, FALSE))
        with(d, {
                c(denom = mean(d$COPDd, na.rm = TRUE),
                  temp = mean(tmpd, na.rm = TRUE),
                  pm25 = mean(pm25tmean, na.rm = TRUE),
                  pm10 = mean(pm10tmean, na.rm = TRUE),
                  o3 = mean(o3tmean, na.rm = TRUE))
        })
})

varsdf <- data.frame(fips = sites, do.call("rbind", vars))
varsdf <- transform(varsdf, denom = round(denom, 1))

sitedata <- merge(sitedata, varsdf, by = "fips")


## Write out to file
siteInfo <- sitedata
dput(siteInfo, file = "siteInfo.R")
