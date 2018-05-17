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
    
    input <- counts(sce)
    FitGoM(t(input), K = cluster_number, tol = 0.1, path_rda = "S10X.FitGoM.rda")

    load("S10X.FitGoM.rda")
    
    eval(parse(text=paste0("colData(sce)$countClust <- unlist(apply(Topic_clus_list$clust_", cluster_number, "$omega, 1,  which.max))")))
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()