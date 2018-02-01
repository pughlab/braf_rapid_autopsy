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