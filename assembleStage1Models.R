################################################################################
## Get results from 'stage1' directory and insert into data frame

basedir <- "stage1"
models <- .readRDS(file.path(basedir, "models.rds"))

beta <- character(nrow(models))
var <- character(nrow(models))

makePath <- function(models, i) {
        file.path(as.character(models[i, "outcome"]),
                  as.character(models[i, "poll"]),
                  as.character(models[i, "siteName"]),
                  as.character(models[i, "d0"]))
}

db <- dbInit(file.path(basedir, "stage1.db"), "DB1")

for(i in seq_len(nrow(models))) {
        key <- makePath(models, i)
        results <- dbFetch(db, key)
        if(!inherits(results, "condition")) {
                beta[i] <- rawToChar(serialize(results$beta,NULL,ascii=TRUE))
                var[i] <- rawToChar(serialize(results$var,NULL,ascii=TRUE))
        }
        else {
                beta[i] <- NA
                var[i] <- NA
        }
}

models$beta <- beta
models$var <- var

attr(models, "out.attrs") <- NULL

.saveRDS(models, file = file.path("results", "pm25.stage1.dist-lag.rds"))
