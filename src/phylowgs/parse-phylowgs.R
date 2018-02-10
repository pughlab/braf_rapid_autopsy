## Parse PhyloWGS
library(jsonlite)
library(ReporteRs)
library(plotrix)
library(data.table)
library(gtools)
library(RColorBrewer)
library(rprojroot)

root_dir = rprojroot::find_rstudio_root_file()

config.path <- paste(root_dir, "/src/config_phylowgs.R",sep="")
source(config.path, local = TRUE)

samples <- processFile(samples)
all.patients <- list()



## prep sample names
for ( i in 1:length(samples)){
  assign(paste(sample.prefix,i,sep=""), samples[[i]])
  all.patients <-list.append(all.patients, get(paste(sample.prefix, i, sep="")))
}

mega.maf <- c()
for (j in 1:length(all.patients)){

  pt.maf <- c()
  maf_name <- paste(patient.prefix, j,".maf",sep="")
  patient_name <- paste(patient.prefix, j, sep="")
  
  df.patient.all <- c()
  tumors <- all.patients[[j]]
  mutass.file = paste(witness.data,"/data/",patient_name,"/",patient_name,".mutass.zip", sep="")
  ssms.file = paste(witness.data,"/data/ssms/",patient_name,"_ssms.txt", sep="")
  cnvs.file = paste(witness.data,"/data/cnvs/",patient_name,"_cnv.txt", sep="")
  maf.file = paste(maf.dir,maf_name, sep="/")
  ccf.file <- paste(ccf.dir,"/" ,patient_name, "_CCF.txt",sep="")
  
  if(recalculate_indecies){
    # find top 50 indecies
    top.indeces <- findTopIndecies(mutass.file ,ssms.file , cnvs.file, maf.file)
    private.snvs <- findPrivateSNVs(maf.file)
    
    best.index <- calcInconsistenyScore(mutass.file ,top.indeces, private.snvs, ssms.file )
    print(paste('******* Best index for ', patient_name, ':',best.index))
  }else{
    best.index = index[patient_name]
    
  }

  BRAF_Nodes  <- parsePhyloWGS(mutass.file , best.index, ssms.file, cnvs.file, maf.file)
  BRAF_Nodes_ssms <- BRAF_Nodes[[1]]
  BRAF_Nodes_cnvs <- BRAF_Nodes[[2]]
  
  annotateCNVs(BRAF_Nodes_cnvs, patient_name, output.dir)
  
  tmp.maf <- generateCBio(tumors, BRAF_Nodes_ssms, BRAF_Nodes_cnvs, maf.file, ccf.file)
  mega.maf <- rbind(mega.maf,tmp.maf)
}

