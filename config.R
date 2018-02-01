src.path <- paste(getwd(), "/src/phylowgs/ssSignature_Local_Functions.R",sep="" )
utils.path <- paste(getwd(), "/src/utils.R",sep="" )
input.dir <- paste(getwd(), "/inputs",sep="" )
samples <- paste(input.dir,'/samples.txt',sep="")
sample.prefix <- "samples.pt"
patiet.prefix <- "Patient"

# phylowgs best indecies
index   <- c("589","2144","2097","2197","1971","2499","1250")
maf.dir <- paste(getwd(), "/inputs/mafs",sep="" )

source(src.path, local = TRUE)
source(utils.path, local = TRUE)
