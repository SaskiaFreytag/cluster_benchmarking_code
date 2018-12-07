library(CountClust)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  
  if (length(args)!=3) {
    
    stop("Three arguments must be supplied (input file, output file, cluster number).n", call.=FALSE)
    
  } else if (length(args)==3) {
    
    filename <- args[1]
    outname <- args[2]
    cluster_number <- as.numeric(args[3])
    
    load(args[1])
    
    input <- as.matrix(counts(sce))
    outs <- FitGoM(t(input), K = cluster_number, tol = 0.1, path_rda = NULL)
  
    
    colData(sce)$countClust <- unlist(apply(outs$fit$omega, 1,  which.max))
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()