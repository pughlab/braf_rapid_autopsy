src.path <- paste(getwd(), "/src/phylowgs/ssSignature_Local_Functions.R",sep="" )
utils.path <- paste(getwd(), "/src/utils.R",sep="" )
input.dir <- paste(getwd(), "/inputs",sep="" )
ref.dir <- paste(getwd(), "/refs",sep="" )

samples <- paste(input.dir,'/samples.txt',sep="")
cancer.genes.path <-  paste(ref.dir,'/cancer.bed',sep="")

sample.prefix <- "samples.pt"
patiet.prefix <- "Patient"

# phylowgs best indecies
index   <- c("589","2144","2097","2197","1971","2499","1250")

samplelist<-c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7")
names(index) <- samplelist


maf.dir <- paste(getwd(), "/inputs/mafs",sep="" )
witness.data <- "~/witness_dec6" #hard-coded for now


source(src.path, local = TRUE)
source(utils.path, local = TRUE)