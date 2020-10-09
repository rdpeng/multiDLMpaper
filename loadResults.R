######################################################################
## Load the distributed lag stage 1 results

loadResults <- local({
        results0 <- .readRDS("results/pm25.dist-lag.rds")
        sinfo <- dget("siteInfo.R")
        slist <- dget("siteList.R")
        
        function(outcome, nkeep = Inf, model = "DLag 0:13", dfTime = 8,
                 env = parent.frame()) {
                if(!is.character(outcome))
                        stop("'outcome' should be character")
                if(!(outcome %in% levels(results0$outcome)))
                        stop(sprintf("outcome '%s' invalid", outcome))

                use <- (results0$model == model & results0$dfTime == dfTime
                        & results0$outcome == outcome & !is.na(results0$beta)
                        & !is.na(results0$var) & results0$fips %in% slist)
                results <- results0[use, ]
                
                b <- lapply(results$beta, function(x) unname(unserialize(x)))
                v <- lapply(results$var, function(x) unname(unserialize(x)))
                names(b) <- names(v) <- as.character(results$fips)
                
                nkeep <- min(nkeep, length(b))
                ord <- order(sinfo$pop100, decreasing = TRUE)
                fipsList <- sinfo[ord, "fips"][seq(nkeep)]
                b <- b[fipsList]
                v <- v[fipsList]
                
                isNULL <- sapply(b, is.null)
                
                assign("b", b[!isNULL], env)
                assign("v", v[!isNULL], env)
        }
})
