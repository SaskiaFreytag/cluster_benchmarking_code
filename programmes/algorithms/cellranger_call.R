library(SingleCellExperiment)


main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)<3) {
    
    stop("Three argument must be supplied (input file, output file, genes).n", call.=FALSE)
    
  } else if (length(args)>=3) {
    
    filename <- args[1]
    outname <- args[2]
    third_parameter <- args[3]
    if(length(args)==4){
        fourth_parameter <- args[4]
    } else {
        fourth_parameter <- ""
    }
    
    system(paste0("module add cellranger \n", "cellranger reanalyze --id=sample_tmp --matrix=sample333/outs/filtered_gene_bc_matrices_h5.h5 ",
    third_parameter, " ", fourth_parameter))
    
    load(filename)

    colData(sce)$`Cell Ranger` <- read.csv("sample_tmp/outs/analysis/clustering/graphclust/clusters.csv", header=T)[,2]
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
    system("/bin/rm -f -r sample_tmp")
  }
  
}

main()


