library(cidr)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    scPBMCs <- scDataConstructor(as.matrix(assay(sce)))
    scPBMCs <- determineDropoutCandidates(scPBMCs)
    scPBMCs <- wThreshold(scPBMCs)
    scPBMCs <- scDissim(scPBMCs)
    scPBMCs <- scPCA(scPBMCs)
    scPBMCs <- nPC(scPBMCs)
    set.seed(89549585)
    scPBMCs <- scCluster(scPBMCs, nPC=5)
    
    colData(sce)$CIDR <- scPBMCs@clusters
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
