library(SingleCellExperiment)
library(RaceID)


main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)<2) {
    
    stop("Two argument must be supplied (input file, output file).n", call.=FALSE)
    
  } else if (length(args)>=2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(filename)
    
    sc <- SCseq(as.data.frame(as.matrix(counts(sce))))
    sc <- filterdata(sc, mintotal=1, minexpr = 5, minnumber = 5,
                     LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                     bmode = "RaceID")
    sc <- compdist(sc, metric="pearson", FSelect = TRUE, knn = NULL)
    sc <- clustexp(sc, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
                   bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
    sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                       outdistquant = 0.95)
    
    res <- colData(sce)
    res$RaceID3 <- sc@cpart
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()