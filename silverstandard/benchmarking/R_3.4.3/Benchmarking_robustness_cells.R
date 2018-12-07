library(BiocParallel)
library(SingleCellExperiment)

set.seed(162622829)
subsampled <- sapply(1:5, function(x) sample(colData(sce_general)$barcode, 3000))

for (i in 1:5) {
  
  load("Sce_Dataset3.RData")
  
  output_link <- paste0("Clustering_Results/Stability/Cells/Iteration_",i, "/")
  system(paste0("mkdir ", output_link))
  algorithm_link <- "../../programmes/algorithms/"
  input_data <- "Sce_Dataset3_Cells.RData"
  final_output <- "Clustering_Result_Stability.RData"
  
  sce <- sce[,subsampled[,i]]
  save(sce, file=input_data)
  
  ## Save list of genes for reanalysis with cellranger
  barcode<-subsampled[,i]
  barcode<-as.matrix(barcode)
  colnames(barcode)<-"Barcode"
  write.csv(barcode, file=paste0("cleaned_barcodes_supsampled.csv"), row.names = F, quote=F)
  
  third_param <- c("", "", 8, "", "", "", "", "../prior_info/BulkRNASeq_Marker_Genes.RData",
                   "", 8, "", "--barcode=cleaned_barcodes_supsampled.csv")
  algorithm <- c("ascend", "cidr", "CountClust", "RaceID", "RaceID2", "RCA", "sc3", "scran", "Seurat", "SIMLR", "TSCAN", "cellranger")
  
  arguments <- list(output_link = output_link, algorithm_link = algorithm_link, algorithm = algorithm, input_data = input_data,
                    third_param = third_param)
  
  call_algo <- function(x, arguments) {
    
    system(paste0( "Rscript --vanilla " , arguments$algorithm_link ,
                                arguments$algorithm[x], "_call.R ", arguments$input_data,
                                " ", arguments$output_link, arguments$algorithm[x], ".txt ", 
                                arguments$third_param[x]))
  }
  
  multicoreParam <- MulticoreParam(workers = 13)
  bplapply(1:12, call_algo, arguments=arguments,  BPPARAM = multicoreParam)
  
  
  system(paste0("/bin/rm -f ", input_data))
  system(paste0("/bin/rm -f cleaned_barcodes_supsampled.csv"))
}

res <- list(ascend=list(), cidr=list(), CountClust=list(), RaceID=list(), RaceID2=list(), RCA=list(), sc3=list(), scran=list(), Seurat=list(), SIMLR=list(),
TSCAN=list(), cellranger=list())

for(i in 1:5){
    output_link <- paste0("Clustering_Results/Stability/Cells/Iteration_",i, "/")
    
    for (x in algorithm){
        par_res <- read.table(paste0(output_link, x, ".txt"), header=TRUE)
        par_res <- par_res[,dim(par_res)[2]]
        eval(parse(text=paste0("res$", x, "[[", i, "]]", "<- par_res")))
    }
    
}

save(res, file=final_output)


