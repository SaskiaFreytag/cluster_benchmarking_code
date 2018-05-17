library(devtools)
install_github("GIS-SP-Group/RCA")
library(RCA)
library(scater)

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
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    sce1 <- pre_clean(sce)
    
    obj <- as.matrix(exprs(sce1))
    rownames(obj) <- rowData(sce1)$symbol
    data_obj <- dataConstruct(obj)
    data_obj <- geneFilt(obj_in = data_obj)
    data_obj <- cellNormalize(data_obj)
    data_obj <- dataTransform(data_obj,"log10")
    data_obj <- featureConstruct(data_obj,method = "GlobalPanel")
    set.seed(20742579)
    data_obj <- cellClust(data_obj)
    
    colData(sce)$RCA <- data_obj$group_labels_color$groupLabel
    res <-colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()