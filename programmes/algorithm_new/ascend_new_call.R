library(BiocParallel)
ncores <- 10
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)
library(ascend)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {
    
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    expression_matrix <- counts(sce)
    rownames(expression_matrix) <- rowData(sce)$id
    
    row_info_data_frame <- as.data.frame(rowData(sce)[, c('id', 'symbol'), drop=FALSE])
    rownames(row_info_data_frame) <- row_info_data_frame$id
    col_info_data_frame <- as.data.frame(colData(sce)[, c('barcode'), drop=FALSE])
    colnames(col_info_data_frame) <- "cell_barcode"
    rownames(col_info_data_frame) <- col_info_data_frame$cell_barcode
    
    control_list <- list(Mt = row_info_data_frame$id[grep("^MT", row_info_data_frame$symbol, ignore.case = TRUE)],
                     Rb = row_info_data_frame$id[grep("^Rps|^Rpl", row_info_data_frame$symbol, ignore.case = TRUE)])
    
    em.set <- newEMSet(assays = list(counts = expression_matrix),
                       colInfo = col_info_data_frame,
                       rowInfo = row_info_data_frame,
                       controls = control_list)
    
    em.set <- filterByControl(em.set, control = 'Mt', pct.threshold = 100)
    em.set <- filterByControl(em.set, control = 'Rb', pct.threshold = 100)
    em.set <- filterLowAbundanceGenes(em.set, pct.threshold = 1)
    em.set <- scranNormalise(em.set, quickCluster = FALSE, min.mean = 1e-05)
    
    set.seed(737368)
    em.set <- runPCA(em.set)
    
    em.set <- runCORE(em.set, dims = 20, nres = 40)
    
    colData(sce)$ascend <- em.set@clusterAnalysis$clusters
    res <-colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()


