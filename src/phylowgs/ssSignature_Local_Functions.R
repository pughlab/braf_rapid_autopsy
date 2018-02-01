### Single Sample mutatation spectrum analysis (ssMSA)
library(SomaticSignatures)
library(deconstructSigs)  ## Load thid package second. Order is important because function "signatures" exists in both.
library(grid)
library(Rmisc)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gridExtra)
library(gplots)

ssSignature <- function(mutect.file,
                        output_dir,
                        pat = ".call_stats.keep",
                        input_tool = c("mutect","oncotator"),
                        samplename="",
                        mut_threshold = 50, ## 50 arbitrary number suggested from deconstructSigs paper
                        ppt_pie = FALSE,
                        pie_rad = 0.8,
                        tri.counts.method = "V5UTR",
                        just_pie = FALSE,
                        hg38 = FALSE) { 
  
  
  ## Add data frame of trinucleotide 
  if(tri.counts.method == "V5UTR"){
    tri.counts.method <- read.table("~/Desktop/PIES/tri.counts.exomeUTRv5.txt",
                                    sep="\t",
                                    header=TRUE)
    
    #print("Normalizing to Sureselect Exome+UTR V5 calculated tri-nucleotide occurences.")
  }
  
  
  #if (is.character(mutect.file)) cat(paste("Processing", names(mutect.file),"\n"))
  #if (is.data.frame(mutect.file) & tolower(input_tool) == "oncotator") cat(paste("Processing",mutect.file$Tumor_Sample_Barcode[1],"\n"))
  
  ### Manually currated signature definitions from http://cancer.sanger.ac.uk/cosmic/signatures
  
  etiologies_loc <- "~/Desktop/PIES/signature_etiologies.txt";
  
  etiologies <- read.table(etiologies_loc,
                           sep="\t",header=TRUE);
  
  ## Add a mutation frequency cut-off to separate low and high allelic fraction variants
  
  ## Can input either a file path or a data frame.  If path, read the data in.
  if(!is.data.frame(mutect.file)){
    raw_mutect <- read.table(mutect.file,
                             sep="\t",header=TRUE);
    if (samplename ==""){
      samplename <- gsub(pat,"",basename(mutect.file))
    }
    
  } else {
    
    raw_mutect <- mutect.file
    
    if (samplename ==""){
      if(tolower(input_tool) == "oncotator"){
        samplename <- gsub(pat,"",mutect.file$Tumor_Sample_Barcode[1])
      } else samplename <- gsub(pat,"",mutect.file$tumor_name[1])
    }
  }
  
  ## Can input either a raw mutect file or a mutect file that has been annotated.
  if(tolower(input_tool) == "oncotator"){
    mutect_cols <- c("contig","position","context","ref_allele","alt_allele","tumor_name",
                     "normal_name","score","dbsnp_site","covered","power","tumor_power",
                     "normal_power","total_pairs","improper_pairs","map_Q0_reads","t_lod_fstar","tumor_f",
                     "contaminant_fraction","contaminant_lod","t_ref_count","t_alt_count","t_ref_sum","t_alt_sum",
                     "t_ref_max_mapq","t_alt_max_mapq","t_ins_count","t_del_count","normal_best_gt","init_n_lod",
                     "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum","judgement")
    
    colnames(raw_mutect)[colnames(raw_mutect)=="Chromosome"] <- "contig"
    colnames(raw_mutect)[colnames(raw_mutect)=="Start_position"] <- "position"
    colnames(raw_mutect)[colnames(raw_mutect)=="Reference_Allele"] <- "ref_allele"
    colnames(raw_mutect)[colnames(raw_mutect)=="Tumor_Seq_Allele2"] <- "alt_allele"
    colnames(raw_mutect)[colnames(raw_mutect)=="Tumor_Sample_Barcode"] <- "tumor_name"
    colnames(raw_mutect)[colnames(raw_mutect)=="Matched_Norm_Sample_Barcode"] <- "normal_name"
    
    raw_mutect$contig <- paste("chr",raw_mutect$contig,sep="")
    
    raw_mutect <- raw_mutect[,mutect_cols]
  }
  
  ## remove REJECT samples
  
  raw_mutect <- subset(raw_mutect,judgement=="KEEP")
  
  raw_mutect$cut_off_tumor_f<-rep("< 0.1 VAF",nrow(raw_mutect));
  raw_mutect$cut_off_tumor_f[raw_mutect$tumor_f >= 0.1] <- ">= 0.1 VAF";
  
  ## Check if the number of mutations meets the set threshold
  if (nrow(raw_mutect) >= mut_threshold) {
    new_mutect.file <- paste(output_dir,
                             "tmp_Mutect_file_w_VAF_cutoff.txt",
                             sep="/");
    
    write.table(raw_mutect,
                new_mutect.file,
                sep="\t",
                row.names=FALSE,
                quote=FALSE);
    
    maf <- readMutect(new_mutect.file);
    
    if(hg38){
      sca_motifs <- mutationContext(maf, 
                                    BSgenome.Hsapiens.UCSC.hg38,
                                    k = 3,  
                                    strand = FALSE, 
                                    unify = TRUE, 
                                    check = TRUE);
    }else{
    sca_motifs <- mutationContext(maf, 
                                  BSgenome.Hsapiens.UCSC.hg19,
                                  k = 3,  
                                  strand = FALSE, 
                                  unify = TRUE, 
                                  check = TRUE);
    }
    
    sca_mat <- motifMatrix(sca_motifs, group = "cut_off_tumor_f");
    
    sca_full <- motifMatrix(sca_motifs);
    
    ## Reformat for compatibility with deconstructSigs package
    sig_in <- as.data.frame(t(sca_mat));
    
    ref <- substr(colnames(sig_in),1,1);
    alt <- substr(colnames(sig_in),2,2);
    first <- substr(colnames(sig_in),4,4);
    third <- substr(colnames(sig_in),6,6);
    
    colnames(sig_in) <- paste(first,"[",ref,">",alt,"]",third,sep="");
    
    
    sig_in_full <- as.data.frame(t(sca_full));
    
    ref <- substr(colnames(sig_in_full),1,1);
    alt <- substr(colnames(sig_in_full),2,2);
    first <- substr(colnames(sig_in_full),4,4);
    third <- substr(colnames(sig_in_full),6,6);
    
    colnames(sig_in_full) <- paste(first,"[",ref,">",alt,"]",third,sep="");
    
    
    ## fill in missing contexts
    miss <- setdiff(colnames(randomly.generated.tumors), colnames(sig_in));
    
    add <- matrix(data = 0, 
                  nrow = nrow(sig_in), 
                  ncol = length(miss), 
                  byrow = FALSE,
                  dimnames = list(row.names(sig_in),
                                  miss))
    
    sig_in <- cbind(sig_in,add)
    
    
    miss <- setdiff(colnames(randomly.generated.tumors), colnames(sig_in_full));
    
    add <- matrix(data = 0, 
                  nrow = nrow(sig_in_full), 
                  ncol = length(miss), 
                  byrow = FALSE,
                  dimnames = list(row.names(sig_in_full),
                                  miss))
    
    sig_in_full <- cbind(sig_in_full,add)
    
    ## Calculate signature correlations separately for <0.1 and >=0.1 variant allele fraction and Full compliment of muts
    
    s_full <- whichSignatures(tumor.ref = sig_in_full, sample.id=rownames(sig_in_full),
                              signatures.ref = signatures.nature2013[1:21,], associated = c(),
                              signatures.limit = 10, signature.cutoff = 0.06, contexts.needed = TRUE,
                              tri.counts.method = tri.counts.method);
    
    if(just_pie){
      return(s_full)
    }
    
    
    #print(s_full)
    #print("stop here if just_pie!!!")
    
    #############################################################################################################
    #############################################################################################################
    #############################################################################################################
    
    
    # full_Weights <- s_full$weights
    # 
    # 
    # if(grepl("< 0.1 VAF",row.names(sig_in))[1]) {
    #   s1 <- whichSignatures(tumor.ref = sig_in, sample.id = "< 0.1 VAF",
    #                         signatures.ref = signatures.nature2013[1:21,], associated = c(),
    #                         signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
    #                         tri.counts.method = tri.counts.method);
    #  
    #   VAF_low_Weights <- s1$weights
    #   
    #   no_low <- FALSE
    #   
    # } else no_low <- TRUE
    # 
    # s2 <- whichSignatures(tumor.ref = sig_in, sample.id = ">= 0.1 VAF",
    #                       signatures.ref = signatures.nature2013[1:21,], associated = c(),
    #                       signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
    #                       tri.counts.method = tri.counts.method);
    # 
    # VAF_high_Weights <- s2$weights
    # 
    # ## Only keep signatures with positive values
    # 
    # full_Weights <- full_Weights[,full_Weights[1,]>0, drop=FALSE]
    # 
    # 
    # if(!no_low) VAF_low_Weights <- VAF_low_Weights[,VAF_low_Weights[1,]>0, drop=FALSE]
    # VAF_high_Weights <- VAF_high_Weights[,VAF_high_Weights[1,]>0, drop=FALSE]
    # 
    # if(!no_low) if(class(VAF_low_Weights) == "data.frame") VAF_low_Weights <- VAF_low_Weights[order(VAF_low_Weights[1,], decreasing=TRUE)]
    # if(class(VAF_high_Weights) == "data.frame") VAF_high_Weights <- VAF_high_Weights[order(VAF_high_Weights[1,], decreasing=TRUE)]
    
    ## Output some useful plots
    
    ## original and reconstructed signatures
    
    
    
    # pdf(paste(output_dir,"/", samplename,
    #           "_Reconstructed_signature_ALL_mutations.pdf",sep=""),
    #     height=11,width=8.5)
    # deconstructSigs::plotSignatures(s_full)
    # dev.off()
    # 
    # if(!no_low){
    #   pdf(paste(output_dir,"/", samplename,
    #             "_Reconstructed_signature_low_VAF_mutations.pdf",sep=""),
    #       height=11,width=8.5)
    #   deconstructSigs::plotSignatures(s1)
    #   dev.off()
    # }
    # 
    # pdf(paste(output_dir,"/", samplename,
    #           "_Reconstructed_signature_high_VAF_mutations.pdf",sep=""),
    #     height=11,width=8.5)
    # deconstructSigs::plotSignatures(s2)
    # dev.off()
    # 
    # pdf(paste(output_dir,"/", samplename,
    #           "_Variant_Abundance_Plot.pdf",sep=""),
    #     height=4,width=6)
    # 
    # print(plotVariantAbundance(sca_motifs,"cut_off_tumor_f"))
    # 
    # dev.off()
    # 
    
    ## Include proposed etiologies for matching signatures
    
    # row.names(etiologies) <- etiologies$Signature
    # 
    # if(!no_low) { 
    #   m_etiologies <- etiologies[unique(c(colnames(VAF_high_Weights),
    #                                       colnames(VAF_low_Weights))),]
    # } else {
    #   m_etiologies <- etiologies[unique(c(colnames(VAF_high_Weights))),]
    # }
    # 
    # full_etiologies <- etiologies[unique(c(colnames(full_Weights))),]
    # 
    # 
    # ## Bug where pie cannot be drawn with just one signature
    # if(!no_low) if(sum(s1$weight==1)==1) s1$unknown <- 0.0001;
    # if(sum(s2$weight==1)==1) s2$unknown <- 0.0001;
    
    
    #if(!ppt_pie){  
      
      # pdf(paste(output_dir,"/", samplename,
      #           "_Mutation_Signature_Contributions_ALL_Mutations.pdf",sep=""),
      #     height=8.5,width=11)
      # layout(matrix(c(1,0,1,0), 2))
      # custom_makePie(s_full, radius = pie_rad)
      # pushViewport(viewport(y=.25,height=.5))
      # grid.table(full_etiologies, rows = rep("",nrow(full_etiologies)))
      # dev.off()
      # 
      # pdf(paste(output_dir,"/", samplename,
      #           "_Mutation_Signature_Contributions_VAF_Separated.pdf",sep=""),
      #     height=8.5,width=11)
      # layout(matrix(c(1,0,2,0), 2))
      # if(!no_low) makePie(s1)
      # makePie(s2)
      # pushViewport(viewport(y=.25,height=.5))
      # grid.table(m_etiologies, rows = rep("",nrow(m_etiologies)))
      # 
      # dev.off()
    #} else {
      # Create a new powerpoint document
  #     doc <- pptx()
  #     # Add a new slide into the ppt document 
  #     doc <- addSlide(doc, "Two Content" )
  #     plot_function <- function(){
  #       layout(matrix(c(1,0,2,0), 2))
  #       if(!no_low) makePie(s1)
  #       makePie(s2)
  #       pushViewport(viewport(y=.25,height=.5))
  #       grid.table(m_etiologies, rows = rep("",nrow(m_etiologies)))
  #     }
  #     doc <- addPlot(doc, plot_function, vector.graphic = TRUE )
  #     
  #     writeDoc(doc, file = paste(output_dir,"/", samplename,
  #                                "_Mutation_Signature_Contributions_VAF_Separated.pptx",sep=""))
  #     
  #     
  #     doc <- pptx()
  #     # Add a new slide into the ppt document 
  #     doc <- addSlide(doc, "Two Content" )
  #     plot_function <- function(){
  #       # layout(matrix(c(1,0,1,0), 2))
  #       custom_makePie(s_full, radius = pie_rad)
  #       # pushViewport(viewport(y=.25,height=.5))
  #       # grid.table(full_etiologies, rows = rep("",nrow(full_etiologies)))
  #     }
  #     doc <- addPlot(doc, plot_function, vector.graphic = TRUE )
  #     
  #     writeDoc(doc, file = paste(output_dir,"/", samplename,
  #                                "_Mutation_Signature_Contributions_ALL_Mutations.pptx",sep=""))
  #   #}
  #   
  #   if(!no_low) VAF_low <- s1$weights
  #   VAF_high <- s2$weights
  #   
  #   VAF_all <- s_full$weights
  #   
  #   
  #   if(!no_low) row.names(VAF_low) <- samplename
  #   row.names(VAF_high) <- samplename
  #   
  #   VAF_combined <- list(VAF_low = VAF_high,
  #                        VAF_high = VAF_high,
  #                        VAF_all = VAF_all)
  #   
  # } else {
  #   blank <- signatures.nature2013[1:21,1,drop=FALSE]
  #   colnames(blank) <- samplename
  #   blank[,samplename] <- rep(NA,nrow(blank))
  #   blank <- t(blank)
  #   
  #   VAF_combined <- list(VAF_low = blank,
  #                        VAF_high = VAF_high,
  #                        VAF_all = VAF_all)
  }
  
  #return(VAF_combined);
  
}


msSignature <- function(mutect.files,
                        project = "",
                        pat = ".call_stats.keep",
                        input_tool = "mutect",
                        manual_samplenames = "",
                        mut_threshold = 50,
                        output_dir = "",
                        ppt_pie=FALSE,
                        pie_rad = rep(0.8,length(mutect.files)),
                        tri.counts.method = "V5UTR",
                        just_pie = FALSE) {
  

  
  ## Identify signatures in multiple samples using a list of mutect files
  
  
  if(manual_samplenames==""){
    if(!is.data.frame(mutect.files)){
    manual_samplenames <- basename(mutect.files)
    } 
  }
  
  if (exists("combined_weights_low_vaf")){
    rm(combined_weights_low_vaf)
  }
  if (exists("combined_weights_high_vaf")){
    rm(combined_weights_high_vaf)
  } 
  
  # if(is.data.frame(mutect.files)) {
  #   
  #   if(tolower(input_tool)=="oncotator") {
  #     samples <- unique(mutect.files$Tumor_Sample_Barcode)
  #   } else {
  #     samples <- unique(mutect.files$tumor_name)
  #   }
  #   
  #   for (samp in samples) {
  #     
  #     if(tolower(input_tool)=="oncotator") {
  #       mutect.file <- subset(mutect.files,Tumor_Sample_Barcode==samp)
  #     } else {
  #       mutect.file <- subset(mutect.files,tumor_name==samp)
  #     }
  #     
  #     mut_num <- mut_threshold
  #     if(is.data.frame(mutect.file)) mut_num <- nrow(mutect.file)
  #       
  #     if(mut_num >= mut_threshold){
  #     d <- ssSignature(mutect.file,
  #                      output_dir,
  #                      input_tool = input_tool,
  #                      pat = pat,
  #                      mut_threshold = mut_threshold,
  #                      ppt_pie = ppt_pie,
  #                      pie_rad = pie_rad[ind],
  #                      tri.counts.method = tri.counts.method,
  #                      just_pie = just_pie)
  #     if(just_pie) {
  #       tmp <- list(d)
  #       names(tmp) <- samp
  #       if(samp == samples[1]){
  #         pies <- tmp
  #       } else {
  #         pies <- c(pies,tmp)
  #       }
  #     }
  #     
  #     if((!exists("combined_weights_low_vaf"))) {
  #       combined_weights_low_vaf <- d$VAF_low
  #       combined_weights_high_vaf <- d$VAF_high } else {
  #         
  #         combined_weights_low_vaf <- rbind(combined_weights_low_vaf,d$VAF_low);
  #         combined_weights_high_vaf <- rbind(combined_weights_high_vaf,d$VAF_high);
  #         
  #       }
  #     } else if(just_pie) {
  #       tmp <- list(NA)
  #       
  #       names(tmp) <- samp
  #       if(samp == samples[1]){
  #         pies <- tmp
  #       } else {
  #         pies <- c(pies,tmp)
  #       }
  #     }
  #   }
  #   
  # } else {
    
    ind = 0
    
    for(mutect.file in mutect.files) {
      
      ind = ind+1
      
      sample_name <-  manual_samplenames[ind]
      
      if(is.na(sample_name)) sample_name <- ""
      
      mut_num <- mut_threshold
      if(is.data.frame(mutect.file)) mut_num <- nrow(mutect.file)
      
      if(mut_num >= mut_threshold){
        
      d <- ssSignature(mutect.file,
                       output_dir,
                       samplename = sample_name,
                       input_tool = input_tool,
                       pat = pat,
                       mut_threshold = mut_threshold,
                       ppt_pie = ppt_pie,
                       pie_rad = pie_rad[ind],
                       tri.counts.method = tri.counts.method,
                       just_pie = just_pie)
      
      
      if(just_pie) {
        tmp <- list(d)
        names(tmp) <- sample_name
        if(ind==1){
          pies <- tmp
        } else {
          pies <- c(pies,tmp)
        }
      }
      
      if((!exists("combined_weights_low_vaf"))) {
        combined_weights_low_vaf <- d$VAF_low
        combined_weights_high_vaf <- d$VAF_high 
        combined_weights_all_mut <- d$VAF_all } else {
          
          combined_weights_low_vaf <- rbind(combined_weights_low_vaf,d$VAF_low);
          combined_weights_high_vaf <- rbind(combined_weights_high_vaf,d$VAF_high);
          combined_weights_all_mut <- rbind(combined_weights_all_mut,d$VAF_all); 
        }
      } else if(just_pie) {
        tmp <- list(NA)
        names(tmp) <- sample_name
        if(ind == 1){
          pies <- tmp
        } else {
          pies <- c(pies,tmp)
        }
      }
    }
#  } 
  
  if(just_pie){
    #print(names(pies))
    return(pies)
  } else {
    
    ## Number of times signature was represented
    sample_count <- combined_weights_high_vaf > 0
    sample_sums <- apply(sample_count,2,sum,na.rm=TRUE)
    
    sample_sums <- sample_sums[order(as.numeric(sample_sums),decreasing = TRUE)]
    
    names(sample_sums) <- gsub("Signature.","",names(sample_sums))
    
    
    pdf(paste(output_dir,"/",project,"_Signature_occurence_barplot.pdf",sep=""),
        height=5,width=10) 
    barplot(sample_sums,col="deepskyblue3",ylab="Signature Occurence Frequency")
    dev.off()
    
    ## Number of times signature was the most prominent
    max_sig <- apply(combined_weights_high_vaf, 1, function(x) {
      colnames(combined_weights_high_vaf)[which.max(x)] } )
    
    ## remove missing samples (samples that did not have enough mutations for deconvolve signatures)
    tmp <- as.data.frame(unlist(max_sig))
    
    max_sig_count <- summary(factor(tmp[,1]))
    
    max_sig_count <- max_sig_count[order(as.numeric(max_sig_count),decreasing = TRUE)]
    
    names(max_sig_count) <- gsub("Signature.","",names(max_sig_count))
    
    pdf(paste(output_dir,"/",project,"_Signature_most_represented_barplot.pdf",sep=""),
        height=5,width=10) 
    barplot(max_sig_count,col="deepskyblue3",ylab="Signature Most Represented Frequency")
    dev.off()
    
    ## heatmap
    mat <- as.matrix(combined_weights_high_vaf)
    
    mat <- mat[,order(apply(mat,2,mean,na.rm=TRUE),decreasing = TRUE)]
    mat <- mat[order(mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],
                     mat[,6],mat[,7],mat[,8],mat[,9],mat[,10],
                     decreasing = TRUE),]
    
    
    #  bk1 = seq(min(mat),max(mat),length.out=11)
    
    pdf(paste(output_dir,"/",project,"_Mutation_Signature_Weight_heatmap.pdf",sep=""),
        height=10,width=10) 
    
    heatmap.2(
      t(mat),
      scale = 'none',
      Rowv = NULL,
      Colv = NULL,
      sepcolor = "white",
      sepwidth=c(0.05,0.05),
      margins = c(10,10),
      # breaks=bk1,
      na.color = "gray93",
      colsep=1:ncol(mat),
      rowsep=1:nrow(mat),
      main = paste('Mutation Signatures in',project),
      col=colorpanel(256,"white","lightblue","darkblue"),
      dendrogram = "none",
      # RowSideColors=side.colours,
      density.info="none",
      trace="none"
    )
    
    dev.off()
    
    ## Output the weight data
    write.table(combined_weights_all_mut,
                paste(output_dir,"/",project,"_ALL_mutations_signature_weights.txt",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    write.table(combined_weights_low_vaf,
                paste(output_dir,"/",project,"_Low_VAF_mutation_signature_weights.txt",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    write.table(combined_weights_high_vaf,
                paste(output_dir,"/",project,"_High_VAF_mutation_signature_weights.txt",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    return(combined_weights_high_vaf)
  }
}



custom_makePie <- function (sigs.output, sub = "", add.color = NULL, x_pos, y_pos, radius = 0.8) 
{
  ## Development
  weights <- data.frame(sigs.output[["weights"]])
  unknown <- sigs.output[["unknown"]]
  weights$unknown <- unknown
  
  all.sigs <- c("Signature.1", "Signature.1A", "Signature.1B", 
                "Signature.2", "Signature.3", "Signature.4", "Signature.5", 
                "Signature.6", "Signature.7", "Signature.8", "Signature.9", 
                "Signature.10", "Signature.11", "Signature.12", "Signature.13", 
                "Signature.14", "Signature.15", "Signature.16", "Signature.17", 
                "Signature.18", "Signature.19", "Signature.20", "Signature.21", 
                "Signature.R1", "Signature.R2", "Signature.R3", "Signature.U1", 
                "Signature.U2", "Signature.22", "Signature.23", "Signature.24", 
                "Signature.25", "Signature.26", "Signature.27", "Signature.28", 
                "Signature.29", "Signature.30", "unknown")
  all.colors <- c("#023FA5", "#023FA5", "#7D87B9", "#BEC1D4", 
                  "#D6BCC0", "#BB7784", "gold", "#4A6FE3", "#8595E1", 
                  "#B5BBE3", "#E6AFB9", "#E07B91", "#D33F6A", "#11C638", 
                  "#8DD593", "#C6DEC7", "#EAD3C6", "#F0B98D", "#EF9708", 
                  "#0FCFC0", "#9CDED6", "#D5EAE7", "#F3E1EB", "#F6C4E1", 
                  "#F79CD4", "#866097", "#008941", "#A30059", "#F6C4E1", 
                  "#F79CD4", "#866097", "#008941", "#A30059", "#008080", 
                  "#8B0000", "#F4A460", "#663399", "#706563")
  all.colors <- cbind(as.data.frame(all.sigs), as.data.frame(all.colors))
  colnames(all.colors) <- c("signature", "color")
  ind <- which(weights > 0)
  weights <- weights[, ind]
  missing.colors <- colnames(weights)[!colnames(weights) %in% 
                                        all.colors$signature]
  if (!is.null(add.color)) {
    tmp <- data.frame(missing.colors, add.color)
    colnames(tmp) = c("signature", "color")
    all.colors <- rbind(all.colors, tmp)
  }
  sigs.present <- colnames(weights)
  missing.colors <- sigs.present[!sigs.present %in% all.colors$signature]
  if (length(missing.colors) > 0) {
    warning(paste("No color assigned for: \n", paste(missing.colors, 
                                                     collapse = ", "), ".\nTo assign one, use add.color parameter.", 
                  sep = ""))
  }
  colors.sigs.present = all.colors$color[match(sigs.present, 
                                               all.colors$signature)]
  grDevices::palette(as.character(colors.sigs.present))
  if (sub == "") {
    top.title <- rownames(weights)
  }
  if (sub != "") {
    top.title <- paste(rownames(weights), " -- ", sub, sep = "")
  }
  
  #graphics::pie(t(weights), col = factor(colnames(weights), 
  #                                       levels = unique(colnames(weights))), labels = colnames(weights), 
  #              main = top.title,
  #              radius = radius)
  
  add_pie(t(weights), col = factor(colnames(weights), 
                                   levels = unique(colnames(weights))), labels = colnames(weights), 
          main = top.title,
          x_pos = x_pos,
          y_pos = y_pos,
          radius = radius)
}


add_pie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                       init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                       col = NULL, border = NULL, lty = NULL, main = NULL, x_pos, y_pos, ...) 
{
  
  ## Development ##
  
  # x = t(weights)
  # col = factor(colnames(weights), 
  #              levels = unique(colnames(weights)))
  # levels = unique(colnames(weights))
  # labels = colnames(weights)
  # main = top.title
  # x_pos = x_pos
  # y_pos = y_pos
  # radius = radius
  # 
  # edges = 200
  # clockwise = FALSE
  # init.angle = if (clockwise) 90 else 0
  # density = NULL
  # angle = 45
  # border = NULL
  # lty = NULL
  ########
  
  
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) {
    labels <- as.character(seq_along(x))
  } else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
 # plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) {
    xlim <- (pin[1L]/pin[2L]) * xlim
  } else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  
  #plot.window(xlim, ylim, "", asp = 1)
  
  if (is.null(col)) 
    if (is.null(density)) {
      col <-  c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  } else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  if (clockwise) { 
    twopi <--2 * pi
  } else twopi <-  2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
   # list(x = radius * cos(t2p), y = radius * sin(t2p))
    
    ## new line to move x and y positions 
    list(x = (radius * cos(t2p)) + x_pos, y = (radius * sin(t2p) + y_pos+radius))
    
  }
  for (i in 1L:nx) {
   # n <- max(2, floor(edges * dx[i]))
    
    n <- max(floor(edges * dx[i]))
    

    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    
    # polygon(c(P$x, 0), c(P$y, 0), density = density[i], 
    #         angle = angle[i], border = border[i], col = col[i], 
    #         lty = lty[i])
    
    polygon(c(P$x, x_pos), c(P$y, y_pos+radius), density = density[i], 
            angle = angle[i], border = border[i], col = col[i], 
            lty = lty[i])
    
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    # if (!is.na(lab) && nzchar(lab)) {
    #   lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
    #   
    #    text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
    #         adj = ifelse(P$x < 0, 1, 0), ...)
    # }
  }
 # title(main = main, ...)
#  invisible(NULL)
}
