## Parse PhyloWGS
library(jsonlite)
library(ReporteRs)
library(plotrix)
library(data.table)
library(gtools)
library(RColorBrewer)

src.dir <- paste(getwd(), "/src/phylowgs/ssSignature_Local_Functions.R",sep="" )
utils.dir <- paste(getwd(), "/src/utils.R",sep="" )
input.dir <- paste(getwd(), "/inputs",sep="" )
samples <- paste(input.dir,'/samples.txt',sep="")
prefix <- "samples.pt"

source(src.dir, local = TRUE)
source(utils.dir, local = TRUE)

samples <- processFile(samples)
all.patients <- list()

## prep sample names
for ( i in 1:length(samples)){
  assign(paste(prefix,i,sep=""), samples[[i]])
  
  all.patients <-list.append(all.patients, get(paste(prefix,i,sep="")))
}

