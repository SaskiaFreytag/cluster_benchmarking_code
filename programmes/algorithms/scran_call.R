library(scran)
library(scater)
library(dynamicTreeCut)


pre_clean <- function(sce) {
  
  keep_feature <- rowSums(as.matrix(counts(sce) > 0)) > floor(dim(sce)[2]*0.01)
  sce1 <- sce[keep_feature,]
  ave.counts <- rowMeans(as.matrix(counts(sce1)))
  keep <- rowMeans(as.matrix(counts(sce1))) >= 0.05
  sce1 <- sce1[keep,]
  del_gns<-c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  sce1<-sce1[-del_gns,]
  sce1 <- calculateQCMetrics(sce1)
  sce1 <- normalise(sce1)
  
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
    
    index_markers<-pmatch(unlist(marker_genes), rownames(sce1))
    index_markers<-index_markers[-which(is.na(index_markers))]

    high.ab <- calcAverage(sce1) > 1
    clusters <- quickCluster(sce1,  subset.row=high.ab)
    sce1 <- computeSumFactors(sce1, cluster=clusters,  subset.row=high.ab)
    sce1 <- normalise(sce1)

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
