infile <- commandArgs(trailingOnly = TRUE)

source("db2rds.R")

db2rds(infile)
