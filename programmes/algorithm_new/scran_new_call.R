library(scran)
library(scater)
library(dynamicTreeCut)


pre_clean <- function(sce) {
  
  ave.counts <- rowMeans(as.matrix(counts(sce)))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  del_gns <- c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  sce1 <- sce1[-del_gns,]  
  
  return(sce1)
}

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=3) {
    
    stop("Two arguments must be supplied (input file, output file, marker file).n", call.=FALSE)
    
  } else if (length(args)==3) {
    
    filename <- args[1]
    outname <- args[2]
    makerfile <- args [3]
    
    load(args[1])
    load(args[3])
    
    sce1 <- pre_clean(sce)
    
    index_markers <- pmatch(unlist(marker_genes), rownames(sce1))
    index_markers <- index_markers[-which(is.na(index_markers))]
    
    clusters <- quickCluster(sce1)
    sce1 <- computeSumFactors(sce1, cluster=clusters)
    
    sce1 <- normalize(sce1)

    chosen.exprs <- as.matrix(exprs(sce1)[index_markers, ])
    my.dist <- dist(t(chosen.exprs))
    set.seed(13115645)
    my.tree <- hclust(my.dist, method = "ward.D2")
    my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))

    colData(sce)$scran <- my.clusters
    res <-colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
