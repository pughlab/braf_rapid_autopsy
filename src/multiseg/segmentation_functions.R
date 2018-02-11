library(falconx)
centromeres = data.frame(chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                                        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),
                         start.pos = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935, 51644205,
                                       34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129, 13000000, 58632012, 10104553),
                         end.pos   = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679, 42254935, 54644205, 37856694, 19000000, 
                                       19000000, 20000000, 38335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012, 13104553))

falcon.seg.seqz <- function(data.file, chromosome, centromeres, seqz){
  
  require(sequenza)
  seqz <-seqz[seqz$chromosome == chromosome,]
  chrom <- gsub(chromosome, pattern = "chr", replacement = "")
  centromeres$chromosome <- gsub(centromeres$chromosome, pattern = "chr", replacement = "")
  seqz   <- seqz[seqz$zygosity.normal == "het", ]
  get.tauhat <- function(seqz, ...) {
    require(falconx)
    at     <- round(seqz$depth.tumor * seqz$Af, 0)
    bt     <- round(seqz$depth.tumor * seqz$Bf, 0)
    an     <- round(seqz$depth.normal * 0.55, 0)
    bn     <- round(seqz$depth.normal * 0.45, 0)
    
    readmat <- as.data.frame(cbind(an,bn,at,bt))
    biasMatrix <- as.data.frame(cbind(seqz$depth.normal,seqz$depth.tumor))
    colnames(biasMatrix) <- c("sN","sT")
    #getChangepoints.x(at, bt, an, bn, ...)
    getChangepoints.x(readmat, biasMatrix, error = 1e-06,...)
  }
  p <- seqz$position < centromeres$start.pos[centromeres$chromosome == chrom]
  q <- seqz$position > centromeres$end.pos[centromeres$chromosome == chrom]
  pos.p    <- seqz$position[p]
  pos.q    <- seqz$position[q]
  l.p <- length(pos.p)
  l.q <- length(pos.q)
  do.breaks <- function(chrom, tauhat) {
    start.pos <- tauhat
    start.pos[-1] <- tauhat[-1]+1
    data.frame(chrom = chrom,
               start.pos = start.pos[-(length(start.pos))],
               end.pos = tauhat[-1])
  }
  chrom  <- unique(seqz$chromosome)  
  if (l.p > 1 & l.q > 1) {
    tauhat.p <- get.tauhat(seqz[p, ], verbose = FALSE)
    tauhat.p <- c(min(pos.p), pos.p[tauhat.p], max(pos.p))
    tauhat.q <- get.tauhat(seqz[q, ], verbose = FALSE)
    tauhat.q <- c(min(pos.q), pos.q[tauhat.q], max(pos.q))
    breaks.p <- do.breaks(chrom, tauhat.p)
    breaks.q <- do.breaks(chrom, tauhat.q)
    rbind(breaks.p, breaks.q)
  } else if (l.p < 2) {
    tauhat.q <- get.tauhat(seqz[q, ], verbose = FALSE)
    tauhat.q <- c(min(pos.q), pos.q[tauhat.q], max(pos.q))
    do.breaks(chrom, tauhat.q)
  } else if (l.q < 2 ) {
    tauhat.p <- get.tauhat(seqz[p, ], verbose = FALSE)
    tauhat.p <- c(min(pos.p), pos.p[tauhat.p], max(pos.p))
    do.breaks(chrom, tauhat.p)
  } else {
    stop("Segmentation went wrong...")
  }
}



falcon.seg.seqz <- findBreaksFalcon(data.file, chromosome.list ){
  seqz <- read.seqz(data.file)
  brk <- list()
  for (i in chromosome.list){
    brk [[i]] <- falcon.seg.seqz(data.file, chromosome = i, centromeres = centromeres, seqz)
  }
  breaks <- do.call(rbind, brk)
  return(breaks)
}
