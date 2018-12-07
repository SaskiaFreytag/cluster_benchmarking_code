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
    
    scPBMCs <- scDataConstructor(as.matrix(assay(sce)), tagType = "raw")
    scPBMCs <- determineDropoutCandidates(scPBMCs, min1 = 3, min2 = 8,
    N = 2000, alpha = 0.1, fast = TRUE, zerosOnly = FALSE,
    bw_adjust = 1)
    scPBMCs <- wThreshold(scPBMCs, cutoff = 0.5, plotTornado = FALSE)
    scPBMCs <- scDissim(scPBMCs, correction = FALSE, threads = 0,
    useStepFunction = TRUE)
    scPBMCs <- scPCA(scPBMCs, plotPC = FALSE)
    scPBMCs <- nPC(scPBMCs)
    set.seed(89549585)
    scPBMCs <- scCluster(scPBMCs, n = NULL, nCluster = NULL,
    nPC = NULL, cMethod = "ward.D2")
    
    colData(sce)$CIDR <- scPBMCs@clusters
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
