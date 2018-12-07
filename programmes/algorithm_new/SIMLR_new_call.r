library(SIMLR)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
      
      stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
      
  } else if (length(args)==2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    set.seed(345454654)
    cluster_num <- SIMLR_Estimate_Number_of_Clusters(X = as.matrix(counts(sce)), NUMC = 2:20, cores.ratio = 1)
    cluster_num <- (min(which.min(cluster_num$K1), which.min(cluster_num$K2))+1)
    output <-  SIMLR_Large_Scale(X = as.matrix(counts(sce)), c = cluster_num, k = 10, kk = 100)
     
    colData(sce)$SIMLR <- output$y$cluster
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
