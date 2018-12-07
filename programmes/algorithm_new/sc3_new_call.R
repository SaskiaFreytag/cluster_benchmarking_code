library(SC3)
library(scater)


pre_clean <- function(sce) {
  
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  del_gns <- c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  sce1 <- sce1[-del_gns,]  
  sce1 <- calculateQCMetrics(sce1)
  sce1 <- normalize(sce1)
  
  return(sce1)
}

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {

    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    sce1 <- pre_clean(sce)
    
    set.seed(4758747)
    rowData(sce1)$feature_symbol <- rowData(sce1)$symbol
    sce1 <- sce1[!duplicated(rowData(sce1)$feature_symbol), ]
    isSpike(sce1, "ERCC") <- FALSE
    
    sce1 <- sc3_estimate_k(sce1)
    k_est <- sce1@metadata$sc3$k_estimation
    sce1 <- sc3(sce1, ks = k_est, biology = FALSE, n_cores=10,  k_estimator = FALSE, rand_seed=0)
    
    eval(parse(text=paste0("colData(sce)$SC3 <- colData(sce1)$sc3_", k_est, "_clusters")))
    res <-colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
