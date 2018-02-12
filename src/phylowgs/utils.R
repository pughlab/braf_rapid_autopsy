library(rlist)

processFile <- function(filepath) {
  
  con = file(filepath, "r")
  
  lines <- list()
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    line<-unlist(strsplit(line, split = ","))
    lines<-list.append(lines,line)
   
  }
  
  close(con)
  return(lines)
}

#############################################################################################################################

parsePhyloWGS <- function(mutass.file,
                          tree.index,
                          ssms.file="",
                          cnvs.file = "",
                          maf.file="") { 
  
 
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

#############################################################################################################################

annotateCNVs <- function(BRAF_Nodes_cnvs, sample_name, outdir){
  
  annotate.dir <- paste(outdir, '/annotations',sep='')
  dir.create(annotate.dir)
  
  for(k in 1:length(BRAF_Nodes_cnvs)){
    
    linex <- c()
    line2 <- c()
    
    node <- paste("Node_",k,sep="")
    outfile <- paste(annotate.dir,'/',sample_name,"_" ,node,".bed" ,sep="" )
    #absCN <- as.numeric(unlist(BRAF_Nodes_cnvs[[k]]["A"])) + as.numeric(unlist(BRAF_Nodes_cnvs[[k]]["B"]))
    
    line2$chr <- as.numeric(paste(unlist(BRAF_Nodes_cnvs[[k]]["Chromosome"])))
    line2$start <- unlist(BRAF_Nodes_cnvs[[k]]["Start.Position"])
    line2$end <- unlist(BRAF_Nodes_cnvs[[k]]["End.Position"])
    line2$absCN <- as.numeric(unlist(BRAF_Nodes_cnvs[[k]]["A"])) + as.numeric(unlist(BRAF_Nodes_cnvs[[k]]["B"]))
    
    linex <- c(linex,line2)
    linex <- as.data.frame(linex)
    linex <- linex[with(linex, order(linex[,1], linex[,2])), ]
    
    
    if(!is.na(BRAF_Nodes_cnvs[[k]]["Chromosome"])){
      
      write.table(linex , outfile,
                  quote=FALSE, col.names = FALSE, row.names=FALSE, sep="\t")
      
    }
  }
  
}


#############################################################################################################################

findTopIndecies <- function(mutass.file,
                      ssms.file="",
                      maf.file="",
                      cnvs.file = ""
                      , private.ssm){
  
  sum.fn <- gsub("mutass.zip","summ.json",mutass.file)
  
  sum.file <- fromJSON(sum.fn, flatten = TRUE)
  
  sum.file.list <-sum.file$trees[0:2449]
  
  
  df <- c()
  indecies <- length(sum.file$trees) - 1
  for (i in 1:indecies){
    
    if(length(sum.file$trees[[i]]) > 3){
      df <- as.data.frame(rbind(df, c(i,sum.file$trees[[i]]$llh, sum.file$trees[[i]]$branching_index,
                                      sum.file$trees[[i]]$linearity_index, sum.file$trees[[i]]$clustering_index )))
      colnames(df) <- c("index", "llh", "branch_index", "linearity_index", "clustering_index")
    }else{
      
      my.command <- paste("sum.file.list$", "`", i, "`", "$structure$`0`",sep="")
      
      node0.len <- length(eval(parse(text=my.command)))
      
      #no polyclonal accepted
      if(node0.len == 1){
        df <- as.data.frame(rbind(df, c(i,sum.file$trees[[i]]$llh)))
        colnames(df) <- c("index", "llh")
      }
     
    }
    
  }
  num.samples <- length(sum.file$trees[[1]]$populations[[1]]$cellular_prevalence)
  
  total.ssms <- 0
  for (k in 1: length(sum.file$trees[[1]]$populations)){
    total.ssms <- total.ssms + sum.file$trees[[1]]$populations[[k]]$num_ssms
  }
  
  df$nlgLH <- (-df$llh/(total.ssms*num.samples))/log(2)
  df <- df[,c("index", "nlgLH")]
  df1 <- df[order(-df$nlgLH),]
  index <- as.integer(df1$index[[1]])
  
  top50.index <- as.character(df1[1:50,]$index)
  return(top50.index)
  
}

#############################################################################################################################

findPrivateSNVs <- function(maf.file=""){
  if(maf.file!=""){
    maf <- read.table(maf.file,
                      header=TRUE,
                      sep="\t",
                      quote="",
                      stringsAsFactors = FALSE)
    maf$Genome_Position<-paste(maf$Chromosome, maf$Start_Position, sep="_")
    n_occur <- data.frame(table(maf$Genome_Position))
    private_mutations <- n_occur[n_occur$Freq == 1,]
    
  }
  return(as.character(private_mutations[,1]))
}

#############################################################################################################################
calcInconsistenyScore <- function(mutass.file,top.indeces, private.snvs, ssms.file){
  
  mutass_dir <- gsub(".zip","",mutass.file)
  
  ssm_loci <- read.table(ssms.file,
                         header=TRUE,
                         sep="\t",
                         stringsAsFactors = FALSE)
  row.names(ssm_loci) <- ssm_loci$id
  
  private.mut.ids <- ssm_loci[ssm_loci$gene %in% private.snvs,]$id
  df<-c()
  
  for (i in top.indeces){
    jsonfn <- paste(mutass_dir,"/",i,".json",sep="")
    json <- fromJSON(jsonfn, flatten = TRUE)
    
    common<-Reduce(intersect, list(json$mut_assignments$`1`$ssms, private.mut.ids))
    
    inconsistency.score <- length(common)
    #print(inconsistency.score)
    df <- as.data.frame(rbind(df,c(as.integer(i),as.integer(inconsistency.score))))
  }
  df <- df[order(df$V2),]
  
  index <- as.integer(df$V1[[1]])
  return (index)
  
}


#############################################################################################################################
generateCBio <- function( tumors, BRAF_Nodes_ssms, BRAF_Nodes_cnvs, maf.file, ccf.file){
  cbio.maf <- read.table(maf.file, sep="\t", quote="", comment.char="#", header=TRUE)
  cbio.maf$Genome_Position <- paste(cbio.maf$Chromosome, "_",cbio.maf$Start_Position, sep= "")
  cbio.maf[,c("CLONAL_NODE", "CANCER_CELL_FRACTION")] <- 0
  CCF.mat<- read.table(ccf.file)
  
  
  for(k in 1:length(tumors)){
    tumor_barcode <- trimws(tumors[k])
    sample_maf <- cbio.maf[cbio.maf$Tumor_Sample_Barcode == tumor_barcode,]
    
    colnames(CCF.mat) <- c("NodeID", unlist(tumors))
    
    rname <- c()
    for(ii in 1:length(BRAF_Nodes_cnvs)){
      rname <- c(rname, paste("Node",ii,sep=""))
    }
    row.names(CCF.mat) <- rname
    CCF.mat2 <- CCF.mat[,-c(1)]
    CCF.mat3 <- t(CCF.mat2)
    
    for(m in 1:length(BRAF_Nodes_cnvs)){
      node <- paste("Node",m,sep="")
      genes.this.node <- as.list(as.character(unlist(list(BRAF_Nodes_ssms[[m]]["Hugo_Symbol"]))))
      
      if( nrow(cbio.maf[cbio.maf$Tumor_Sample_Barcode == tumor_barcode & cbio.maf$Hugo_Symbol %in% genes.this.node,]) > 0){
        cbio.maf[cbio.maf$Tumor_Sample_Barcode == tumor_barcode & cbio.maf$Hugo_Symbol %in% genes.this.node,]$CANCER_CELL_FRACTION <- as.numeric(CCF.mat3[tumor_barcode, m])
        cbio.maf[cbio.maf$Tumor_Sample_Barcode == tumor_barcode & cbio.maf$Hugo_Symbol %in% genes.this.node,]$CLONAL_NODE <- node
      }
      
    }
    
  }
  return(cbio.maf)
}

removeBottomIndecies<-function(mutass.file, top.indeces){
  mutass_dir <- gsub(".zip","",mutass.file)
  
  files <- list.files(mutass_dir , pattern = ".json$")
  top.indeces <- paste(top.indeces,".json",sep="")
  
  for (fn in files){
    if (!(fn %in% top.indeces)){
      
      file.remove(paste(mutass_dir,"/",fn,sep=""))
    }
  }
  
}