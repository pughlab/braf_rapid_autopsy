## Parse PhyloWGS
library(jsonlite)
library(ReporteRs)
library(plotrix)
library(data.table)
library(gtools)
library(RColorBrewer)

################################################################## Functions begin ##################################################################
parsePhyloWGS <- function(mutass.file,
                          tree.index,
                          ssms.file="",
                          maf.file="",
                          cnvs.file = "") { 
  
  CNV <- !(cnvs.file=="")
  ssm_loci <- read.table(ssms.file,
                         header=TRUE,
                         sep="\t",
                         stringsAsFactors = FALSE)
  row.names(ssm_loci) <- ssm_loci$id
  
  if(maf.file!=""){
    maf <- read.table(maf.file,
                      header=TRUE,
                      sep="\t",
                      quote="",
                      stringsAsFactors = FALSE)
    
    maf$Genome_Position <-paste(maf$Chromosome, maf$Start_Position, sep="_")
    maf <- maf[!duplicated(maf$Genome_Position),]
    row.names(maf) <- paste(maf$Chromosome, maf$Start_Position, sep="_")
  }
  
  if(CNV) {
    cnv_loci <- read.table(cnvs.file,
                           header=TRUE,
                           sep="\t",
                           stringsAsFactors = FALSE)
    
    row.names(cnv_loci) <- cnv_loci$cnv
    
    if(nrow(cnv_loci)==1) {
      warning("No subclonal copy-number alterations. All CNVs are in Node 1.")
      CNV=FALSE
    } else {
      ## parse CNV loci
      
      ## sometimes c0 has multiple CNVs, so this is a temporary work-around
      CNs <- unlist(strsplit(cnv_loci$physical_cnvs[1:nrow(cnv_loci)],","))
      CNs <- unlist(strsplit(CNs,"="))
      CNs <- unlist(strsplit(CNs,";"))
      
      cnv_loci$Chromosome[1:nrow(cnv_loci)] <- CNs[seq(2,length(CNs),by=12)]
      cnv_loci$Start.Position[1:nrow(cnv_loci)] <- CNs[seq(4,length(CNs),by=12)]
      cnv_loci$End.Position[1:nrow(cnv_loci)] <- CNs[seq(6,length(CNs),by=12)]
      cnv_loci$A[1:nrow(cnv_loci)] <- CNs[seq(8,length(CNs),by=12)]
      cnv_loci$B[1:nrow(cnv_loci)] <- CNs[seq(10,length(CNs),by=12)]
    }
  }
  
  mutass_dir <- gsub(".zip","",mutass.file)
  
  if (!file.exists(mutass_dir)) {
    system(paste("mkdir -p",mutass_dir))
    
    setwd(mutass_dir)
    system(paste("unzip",mutass.file))
  }
  
  json <- read.table(paste(mutass_dir,"/",tree.index,".json",sep=""),
                     header=FALSE,
                     sep="\t",
                     quote="",
                     stringsAsFactors = FALSE)
  
  json <- json$V1
  
  json.list <- fromJSON(json)
  
  
  for(i in 1:length(json.list$mut_assignments)) {
    ssms <- json.list$mut_assignments[[i]]$ssms
    if(CNV) {
      cnvs <- json.list$mut_assignments[[i]]$cnvs
      
      if(length(cnvs)>0){
        cnv_loci_list_tmp <- list(cnv_loci[cnvs,])
      } else{
        cnv_loci_list_tmp <- NA
      }
      names(cnv_loci_list_tmp) <- paste("Node_",i,sep="")
    }
    
    tab <- ssm_loci[ssms,]
    
    if(maf.file!="") {
      tab <- maf[tab$gene,]
      tab$gene <- row.names(tab)
    }
    
    ssm_loci_list_tmp <- list(tab)
    names(ssm_loci_list_tmp) <- paste("Node_",i,sep="")
    
    if (i == 1) {
      ssm_loci_list <- ssm_loci_list_tmp
      if(CNV) {
        cnv_loci_list <- cnv_loci_list_tmp
      } 
      
    } else {
      ssm_loci_list <- append(ssm_loci_list,ssm_loci_list_tmp)
      if(CNV) {
        cnv_loci_list <- append(cnv_loci_list,cnv_loci_list_tmp)
      }
    }
    
  }
  
  if(CNV) {
    return(list(ssm_loci_list,cnv_loci_list))
  } else {
    return(ssm_loci_list)
  }
}

################################################################## Functions end ##################################################################
config.path <- paste(getwd(), "/config.R",sep="")
source(config.path, local = TRUE)

samples <- processFile(samples)
all.patients <- list()

## prep sample names
for ( i in 1:length(samples)){
  assign(paste(sample.prefix,i,sep=""), samples[[i]])
  all.patients <-list.append(all.patients, get(paste(sample.prefix, i, sep="")))
}


for (j in 1:length(all.patients)){
  pt.maf <- c()
  maf_name <- paste(patiet.prefix,j,".maf",sep="")
  sample_name <- paste(patiet.prefix,j,sep="")
  
  df.patient.all <- c()
  tumors <- all.patients[[j]]
  
  cbio.maf <- read.table(paste(maf.dir,maf_name, sep="/"), sep="\t", quote="", comment.char="#", header=TRUE)
  cbio.maf$Genome_Position <- paste(cbio.maf$Chromosome, "_",cbio.maf$Start_Position, sep= "")
  
  BRAF_Nodes  <- parsePhyloWGS(mutass.file = paste(witness.data,"/data/",
                                                   sample_name,"/",sample_name,".mutass.zip", sep=""),
                               tree.index = index[sample_name],
                               ssms.file = paste(witness.data,"/data/ssms/",
                                                 sample_name,"_ssms.txt", sep=""),
                               cnvs.file = paste(witness.data,"/data/cnvs/",
                                                 sample_name,"_cnv.txt", sep=""),
                               maf.file=paste(maf.dir,
                                              maf_name, sep="/"))
  
  
  BRAF_Nodes_ssms <- BRAF_Nodes[[1]]
  BRAF_Nodes_cnvs <- BRAF_Nodes[[2]]
}

