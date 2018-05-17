library(BiocParallel)
library(SingleCellExperiment)

final_output <- "Clustering_Robustness_Genes.RData"

set.seed(12345)
random_gene_sets <- replicate(10, sample(1:dim(sce)[1], floor(dim(sce)[1]/2)))

for (i in 1:10) {

  load("Sce_Dataset1.RData")
    
  output_link <- paste0("Clustering_Results/Stability/Genes/Iteration_",i, "/")
  system(paste0("mkdir ", output_link))
  algorithm_link <- "../../programmes/algorithms/"
  input_data <- "Sce_Dataset1_Half.RData"
  
  sce <- sce[random_gene_sets[,i], ]
  save(sce, file=input_data)
  
  ## Save list of genes for reanalysis with cellranger
  gene<-rowData(sce)$id
  gene<-as.matrix(gene)
  colnames(gene)<-"Gene"
  write.csv(gene, file="cleaned_genes_half.csv", row.names = F, quote=F)

  third_param <- c("", "", 8, "", "", "", "", "../prior_info/BulkRNASeq_Marker_Genes.RData",
                   "", 8, "", "--barcode=cleaned_barcodes_daatset1.csv --genes=cleaned_genes_half.csv")
  algorithm <- c("ascend", "cidr", "CountClust", "RaceID", "RaceID2", "RCA", "sc3", "scran", "Seurat", "SIMLR", "TSCAN", "cellranger")

  arguments <- list(output_link = output_link, algorithm_link = algorithm_link, algorithm = algorithm, input_data = input_data,
                  third_param = third_param)

  call_algo <- function(x, arguments) {
  
    system.time(system( paste0( "Rscript --vanilla " , arguments$algorithm_link , 
                                arguments$algorithm[x], "_call.R ", arguments$input_data,
                                " ", arguments$output_link, arguments$algorithm[x], ".txt ", 
                                 arguments$third_param[x]) ) )
  }

  multicoreParam <- MulticoreParam(workers = 13)
  time_taken <- bplapply(1:12, call_algo, arguments=arguments,  BPPARAM = multicoreParam)
  
  
  names(time_taken) <- algorithm
  
  save(time_taken, file=paste0(output_link, "Time_Taken_Half_Iteration_", i, ".RData"))
  system(paste0("/bin/rm -f ", input_data))
  system(paste0("/bin/rm -f cleaned_genes_half.csv"))
  }

 res <- list(ascend=list(), cidr=list(), CountClust=list(), RaceID=list(), RaceID2=list(), RCA=list(), sc3=list(), scran=list(), Seurat=list(), SIMLR=list(), 
             TSCAN=list(), cellranger=list())

for(i in 1:10){
  output_link <- paste0("Clustering_Results/Stability/Genes/Iteration_",i, "/")
  
  for (x in algorithm){
    par_res <- read.table(paste0(output_link, x, ".txt"), header=TRUE)
    par_res <- par_res[,dim(par_res)[2]]
    eval(parse(text=paste0("res$", x, "[[", i, "]]", "<- par_res")))                      
  }
  
} 
 
save(res, file=final_output) 
 
