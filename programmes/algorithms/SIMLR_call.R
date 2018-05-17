library(SIMLR)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=3) {
    
    stop("Three arguments must be supplied (input file, output file and number of clusters).n", call.=FALSE)
  
  } else if (length(args)==3) {
    
    filename <- args[1]
    outname <- args[2]
    cluster_num <- as.numeric(args[3])
    
    load(args[1])
    
    set.seed(345454654)
    output <-  SIMLR_Large_Scale(X = counts(sce), c = cluster_num, k = 10, kk = 100)
     
    colData(sce)$SIMLR <- output$y$cluster
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
