library(devtools)
library(SingleCellExperiment)
source_url("https://raw.githubusercontent.com/dgrun/StemID/master/RaceID2_StemID_class.R")

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=2) {
    
    stop("Two arguments must be supplied (input file and output file).n", call.=FALSE)
    
  } else if (length(args)==2) {
    
    filename <- args[1]
    outname <- args[2]
    
    load(args[1])
    
    sc <- SCseq(as.data.frame(as.matrix(counts(sce))))
    sc <- filterdata(sc, mintotal=1, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
    sc <- clustexp(sc, clustnr=30, bootnr=50, metric="pearson", do.gap=FALSE,sat=TRUE, SE.method="Tibs2001SEmax", 
                   SE.factor=.25, B.gap=50, cln=0, rseed=17000, FUNcluster="kmedoids")
    sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)
    
    colData(sce)$RaceID2 <- sc@cpart
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
  }
  
}

main()