################################################################################
## Collapse datasets

library(APHealth)
library(tsModel)
adata <- .readRDS("data/all-sites.rds")
bdata <- lapply(adata, function(x) {
        data <- collapse(x)
        cat(fips <- as.character(data$fips[1]), "\n")
        if(fips == "24005") {
                data$tmpd <- with(data, (tmax35 + tmin35) / 2)
                data$dptp <- data$dptp35
                data$rmtmpd <- with(data, runMean(tmpd, 1:3))
                data$rmdptp <- with(data, runMean(dptp, 1:3))
        }
        data$IHD <- data$admitp3
        data$IHDd <- data$denomp3
        data$RESPALL <- with(data, admitp9 + admitp8)
        data$RESPALLd <- data$denom
        data[, c("date", "dow", "tmpd", "rmtmpd", "dptp", "rmdptp",
                 "COPDAE", "COPDAEd",
                 "IHD", "IHDd",
                 "RESPALL", "RESPALLd",
                 "admitp5", "denomp5",
                 "CVD", "CVDd",
                 "RESP", "RESPd",
                 "pm25tmean", "fips")]
})

.saveRDS(bdata, file = "data/all-sites-collapsed.rds", compress = TRUE)

