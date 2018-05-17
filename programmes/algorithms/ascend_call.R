library(BiocParallel)
ncores <- 10
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE) # setting the number of cores
library(ascend)
library(SingleCellExperiment)

setGeneric(name = "SyncSlots", def = function(object) {
  standardGeneric("SyncSlots")
})

setMethod("SyncSlots", signature("EMSet"), function(object) {
  # Get data frames
  expression.matrix <- object@ExpressionMatrix
  gene.information <- object@GeneInformation
  cell.information <- object@CellInformation
  
  gene.cols <- colnames(gene.information)
  cell.cols <- colnames(cell.information)
  
  # Expression matrix serves as the basis for updating all the other information
  present.cells <- colnames(expression.matrix)
  present.genes <- rownames(expression.matrix)
  
  # Get data frames based on this info
  gene.information <- gene.information[gene.information[, 1] %in% present.genes, ]
  cell.information <- cell.information[cell.information[, 1] %in% present.cells, ]
  
  # Check in case it's not a dataframe
  if (!is.data.frame(gene.information)){
    gene.information <- data.frame(gene.information)
    colnames(gene.information) <- gene.cols
  }
  
  if (!is.data.frame(cell.information)){
    cell.information <- data.frame(cell.information)
    colnames(cell.information) <- cell.cols
  }
  
  # Put back into the slots
  object@GeneInformation <- gene.information
  object@CellInformation <- cell.information
  
  return(object)
})

FilterByExpressedGenesPerCell <- function (object, pct.value = 1)
{
  filtered.object <- object
  expression.matrix <- filtered.object@ExpressionMatrix
  cells.per.gene <- Matrix::rowSums(expression.matrix)
  remove.genes <- cells.per.gene < (ncol(expression.matrix) *
                                      (pct.value/100))
  if (any(remove.genes)) {
    removed.genes <- names(which(remove.genes))
    remove.idx <- which(rownames(expression.matrix) %in%
                          removed.genes)
    filtered.matrix <- expression.matrix[-remove.idx, ]
  }
  else {
    removed.genes <- list()
    filtered.matrix <- expression.matrix
  }
  filtered.object <- ReplaceExpressionMatrix(filtered.object,
                                             filtered.matrix)
  filtered.object@Log$FilterByExpressedGenesPerCell <- removed.genes
  if (is.null(filtered.object@Log$FilteringLog)) {
    filtered.df <- data.frame()
  }
  else {
    filtered.df <- filtered.object@Log$FilteringLog
  }
  output.list <- list(FilterByExpressedGenesPerCell = length(removed.genes))
  filtered.df <- list(filtered.df, output.list) #changed this as it was throwing errors
  filtered.object@Log$FilteringLog <- filtered.df
  filtered.object@Log$FilterByExpressedGenesPerCell <- list(RemovedGenes = removed.genes)
  filtered.object <- SyncSlots(filtered.object)
  return(filtered.object)
}

pre_clean <- function (sce) {
  
  GeneInformation <- as(rowData(sce), "data.frame")
  rownames(GeneInformation) <- rownames(sce)
  CellInformation <- as(colData(sce), "data.frame")
  rownames(CellInformation) <- colnames(sce)
  emset <- NewEMSet(ExpressionMatrix = as.matrix(counts(sce)) , GeneInformation = GeneInformation, CellInformation = CellInformation)
  controls <- list(Mt = rowData(sce)$id[rowData(sce)$Mito], Rb = rowData(sce)$id[rowData(sce)$ProteinRibo])
  emset <- UpdateControls(emset, controls)
  emset <- FilterByExpressedGenesPerCell(emset, 1)
  normset <- NormaliseByRLE(emset)
  normset <- ExcludeControl(normset, "Mt")
  normset <- ExcludeControl(normset, "Rb")
  
  return(normset)
  
}

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {
    

    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    normset <- pre_clean(sce)
    
    set.seed(737368)
    pcaset <- RunPCA(normset)
    pcaset <- ReduceDimensions(pcaset, n = 20)
    clustered_set <- RunCORE(pcaset)
    
    colData(sce)$ascend <- GetCellInfo(clustered_set)$cluster
    res <-colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()

