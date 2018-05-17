library(TSCAN)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
  
  } else if (length(args)==2) {
  
    filename <- args[1]
    outname <- args[2]
  
    load(args[1])
  
    procdata <- preprocess(as.matrix(counts(sce)), minexpr_percent = 0.01)
    set.seed(4343646)
    lpsmclust <- exprmclust(procdata)
    
    colData(sce)$TSCAN <- lpsmclust$clusterid
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
