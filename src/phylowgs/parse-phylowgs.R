## Parse PhyloWGS
library(jsonlite)
library(ReporteRs)
library(plotrix)
library(data.table)
library(gtools)
library(RColorBrewer)

config.path <- paste(getwd(), "/config.R",sep="")
source(config.path, local = TRUE)

samples <- processFile(samples)
all.patients <- list()

## prep sample names
for ( i in 1:length(samples)){
  assign(paste(prefix,i,sep=""), samples[[i]])
  all.patients <-list.append(all.patients, get(paste(prefix,i,sep="")))
}

