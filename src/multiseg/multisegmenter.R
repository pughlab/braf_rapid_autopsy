library(sequenza)
library(copynumber)
library(plyr)
library(dplyr)
library(data.table)
library(MASS)

library(rprojroot)
root_dir = rprojroot::find_rstudio_root_file()
config.path <- paste(root_dir, "/src/config_multiseg.R",sep="")
source(config.path, local = TRUE)

for (patient in patient.list){
  
  sig.dir.patient <- paste(seg.dir, '/', patient, sep='')
  seq.gz.files <- list.files(sig.dir.patient, pattern=glob2rx("RAP-00*.data.gz"), full.names = TRUE)
  
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  
  seg.list <- list()
  i <-1
  for (data.file in seq.gz.files){
    
    seqz.data <- read.seqz(data.file)
    
    gc.stats <- gc.sample.stats(data.file)
    gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
    seqz.data$adjusted.ratio <- seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)]
    
    d <- seqz.data[,c(1,2,15)]
    d$log2R <- log2(d$adjusted.ratio)
    d$chromosome <- gsub("chr","",d$chromosome)
    d$Genome_Position <- paste(d$chromosome , "_",d$position, sep="")
    
    seg.list[[i]] <- d
    i <- i+1
  }
  
  pt.all <- Reduce(function(x, y) merge(x, y, by = "Genome_Position"), seg.list)
  pt.all <- pt.all[with(pt.all, order(as.numeric(chromosome.x), position.x)), ]
  pt.all <- cbind(pt.all[,c(2,3)], pt.all[, grepl("log",names(pt.all))] )
  colnames(pt.all) <- c("chrom", "pos",gsub(".data.gz","",basename(seq.gz.files)))
  
  pt.wins <- winsorize(data=pt.all, verbose=FALSE)
  multi.seg <- multipcf(data=pt.wins, verbose=FALSE, gamma = 200)
  breaks   <- multi.seg[, c("chrom", "start.pos", "end.pos", "arm")]
  breaks$chrom <- paste("chr",breaks$chrom,sep="")
  
  
  for (data.file in seq.gz.files){
    
    output.dir<-paste(output.dir.all, strsplit(basename(seq.gz.files[[1]]),"_")[[1]][1], sep="/" )
    
    data  <- sequenza.extract(data.file , breaks = breaks)
    CP.example <- sequenza.fit(data, 
                               ploidy = seq(1.5 , 4.5, 0.1))
    sequenza.results(data , CP.example, 
                     sample.id = paste(basename(data.file),sep="") , 
                     out.dir= output.dir )
    
  }
  
  
}



