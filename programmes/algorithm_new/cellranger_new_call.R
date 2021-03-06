library(SingleCellExperiment)


main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)<2) {
    
    stop("Two argument must be supplied (input file, output file).n", call.=FALSE)
    
  } else if (length(args)>=2) {
    
    filename <- args[1]
    outname <- args[2]

    dataset_name <- gsub(".RData", "", strsplit(filename, "_")[[1]][2])
    
    system(paste0("module add cellranger \n", "cellranger reanalyze --id=sample_tmp_", dataset_name, 
                  " --matrix=Cellranger_Clustering/", dataset_name, "_analysis/outs/raw_gene_bc_matrices_h5.h5 ",
                  "--barcodes=cleaned_barcodes_", tolower(dataset_name), ".csv --localcores=10 --localmem=400"))
    
    load(filename)

    colData(sce)$`Cell Ranger` <- read.csv(paste0("sample_tmp_", dataset_name,
                                "/outs/analysis/clustering/graphclust/clusters.csv"), header=T)[,2]
  
    res <- colData(sce)
    
    write.table(res, file=outname, quote=FALSE, row.names = FALSE, col.names=TRUE)
    system(paste0("/bin/rm -f -r sample_tmp_", dataset_name))
  }
  
}

main()


