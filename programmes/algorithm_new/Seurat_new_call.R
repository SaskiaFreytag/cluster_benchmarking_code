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
    
    pbmc <-  CreateSeuratObject(raw.data = as.matrix(counts(sce)),
                                min.cells = 3, min.genes = 0,  project = "10X_PBMC",
                                is.expr = 0, normalization.method = NULL,
                                scale.factor = 10000, do.scale = FALSE, do.center = FALSE,
                                names.field = 1, names.delim = "_", meta.data = NULL,
                                display.progress = FALSE)
    mito.genes <- grep(pattern = "^MT", x = rownames(x = pbmc@data), value = TRUE)
    percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
    pbmc <- AddMetaData(pbmc, percent.mito, col.name = "percent.mito")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000,
                          assay.type = "RNA", display.progress = FALSE)
    pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,
                              y.high.cutoff = Inf, num.bin = 20, binning.method = "equal_width",
                              selection.method = "mean.var.plot", top.genes = 1000, do.recalc = TRUE,
                              sort.results = TRUE, do.cpp = TRUE, display.progress = FALSE)
    pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"), 
                      model.use = "linear", use.umi = FALSE, do.scale = TRUE,
                      do.center = TRUE, scale.max = 10, block.size = 1000,
                      min.cells.to.block = 3000, display.progress = FALSE, assay.type = "RNA",
                      do.cpp = TRUE, check.for.norm = TRUE, do.par = FALSE, num.cores = 1)
    set.seed(635465465)
    pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, reduction.name = "pca",
                   reduction.key = "PC", assay.type = "RNA", seed.use = 42)
    pbmc <- JackStraw(pbmc, num.pc = 20, num.replicate = 100, prop.freq = 0.01,
                      display.progress = FALSE, do.par = FALSE, num.cores = 1, maxit = 1000)
    set.seed(44487152)
    pbmc <- FindClusters(object = pbmc, reduction.type = "pca", 
                         dims.use = 1:12,  print.output = FALSE, save.SNN = FALSE,
                         genes.use = NULL,
                         k.param = 30, plot.SNN = FALSE, prune.SNN = 1/15,
                         distance.matrix = NULL, 
                         reuse.SNN = FALSE, force.recalc = FALSE, nn.eps = 0,
                         modularity.fxn = 1, resolution = 0.8, algorithm = 1, n.start = 100,
                         n.iter = 10, random.seed = 0)
    colData(sce)$Seurat <- pbmc@meta.data$res.0.8
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()
