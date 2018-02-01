src.path <- paste(getwd(), "/src/phylowgs/ssSignature_Local_Functions.R",sep="" )
utils.path <- paste(getwd(), "/src/utils.R",sep="" )
input.dir <- paste(getwd(), "/inputs",sep="" )
samples <- paste(input.dir,'/samples.txt',sep="")
prefix <- "samples.pt"

source(src.path, local = TRUE)
source(utils.path, local = TRUE)
