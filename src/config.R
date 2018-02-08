library(rprojroot)
root_dir = rprojroot::find_rstudio_root_file()


src.path <- paste(root_dir, "/src/phylowgs/ssSignature_Local_Functions.R",sep="" )
utils.path <- paste(root_dir, "/src/phylowgs/utils.R",sep="" )
input.dir <- paste(root_dir, "/inputs",sep="" )
output.dir <- paste(root_dir, "/outputs",sep="" )
ref.dir <- paste(root_dir, "/refs",sep="" )

samples <- paste(input.dir,'/samples.txt',sep="")
cancer.genes.path <-  paste(ref.dir,'/cancer.bed',sep="")

sample.prefix <- "samples.pt"
patient.prefix <- "Patient"

# phylowgs best indecies
index   <- c("589","2144","2097","2197","1971","2499","1250")
recalculate_indecies <- TRUE

samplelist<-c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7")
names(index) <- samplelist


maf.dir <- paste(input.dir, "/mafs",sep="" )
witness.data <- "~/witness_dec6" #hard-coded for now


source(src.path, local = TRUE)
source(utils.path, local = TRUE)