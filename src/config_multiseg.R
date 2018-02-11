library(rprojroot)
root_dir = rprojroot::find_rstudio_root_file()

input.dir <- paste(root_dir, "/inputs/segmentation",sep="" )
maf.dir <- paste(input.dir, "/mafs",sep="" )
seg.dir <- paste(input.dir, "/sequenza",sep="" )
patient.list <- c("RAP1","RAP2","RAP3","RAP4","RAP5","RAP6","RAP7")

output.dir.all <- paste(root_dir, "/outputs/multisegs",sep="" )
