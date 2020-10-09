## Use simple-gibbs.R here

params <- list(outcomeName = c("COPDAE", "heart failure",
               "cerebrovascular disease", "ischemic heart disease"),
               covModel = c("exp", "exp2", "delay"),
               sigmaG = c(0.005, 0.01),
               DLag = c("DLag 0:13", "DLag 0:6"))

makeClusterDirs <- function(n, basedir = "cluster2", offset = 0) {
    if(!file.exists(basedir)) {
        message("creating ", sQuote(basedir))
        dir.create(basedir)
    }
    start <- offset + 1
    end <- offset + n
    dirNames <- file.path(basedir, formatC(seq(start, end), flag = "0", width = 2))
    
    for(i in seq(along = dirNames)) {
        dir.create(dirNames[i])
        file.symlink("~/projects/multiDLM/cluster.sh", dirNames[i])
        file.symlink("~/projects/multiDLM/loadResults.R", dirNames[i])
        file.symlink("~/projects/multiDLM/siteList.R", dirNames[i])
        file.symlink("~/projects/multiDLM/simple-gibbs.R", dirNames[i])
        file.symlink("~/projects/multiDLM/results", dirNames[i])
        file.symlink("~/projects/multiDLM/siteInfo.R", dirNames[i])
    }
    invisible(dirNames)
}


writeParamsFile <- function(ptable, dirnames) {
    for(i in seq(along = dirnames)) {
        con <- file(file.path(dirnames[i], "rungibbs.R"), "w")
        writeLines("start <- proc.time()", con)
        writeLines("source(\"simple-gibbs.R\")", con)
        writeLines("source(\"loadResults.R\")", con)

        writeLines("\n", con)
        for(j in seq(ncol(ptable))) {
            vname <- names(ptable)[j]
            val <- ptable[i, vname]

            if(is.numeric(val))
                writeLines(paste(vname, " <- ", val, "", sep = ""), con)
            else
                writeLines(paste(vname, " <- \"", val, "\"", sep = ""), con)
        }
        writeLines("\n", con)
        writeLines("loadResults(outcomeName, model = DLag)", con)
        writeLines("g <- gibbs(b, v, maxit = 25000, verbose = TRUE, covModel = covModel, sigmaG = sigmaG, do.plot = FALSE)", con)
        writeLines("print(proc.time() - start)", con)
        close(con)
    }
}

start <- function(params, offset = 0) {
    if(is.list(params))
        ptable <- expand.grid(params)
    else if(is.data.frame(params))
        ptable <- params
    else
        stop("'params' should be list or data.frame")
    dirnames <- makeClusterDirs(nrow(ptable), offset = offset)
    writeParamsFile(ptable, dirnames)
    cwd <- getwd()

    for(i in seq(along = dirnames)) {
        cat("Starting", dirnames[i], "\n")
        setwd(dirnames[i])
        system("qsub -cwd cluster.sh")
        setwd(cwd)
        Sys.sleep(1)
    }
    TRUE
}
