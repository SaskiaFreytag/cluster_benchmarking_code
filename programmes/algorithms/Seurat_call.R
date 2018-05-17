library(Seurat)
library(SingleCellExperiment)

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
  
  } else if (length(args)==2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    pbmc <-  CreateSeuratObject(raw.data = as.matrix(counts(sce)), min.cells = 3
    , min.genes = 0,  project = "10X_PBMC")
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
    percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
    pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
    set.seed(635465465)
    pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
    pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
    set.seed(44487152)
    pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:12, resolution = 0.6, print.output = 0, save.SNN = TRUE)
    colData(sce)$Seurat <- pbmc@meta.data$res.0.6
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
