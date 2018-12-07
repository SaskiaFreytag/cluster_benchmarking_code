library(BiocParallel)
library(SingleCellExperiment)

load("Sce_Dataset4.RData")
sce_general <- sce

rank <- order(rowData(sce)$mean_counts, decreasing = TRUE)
nn <- dim(sce)[1]


for (i in 1:5) {
  
  output_link <- paste0("Clustering_Results/Stability/Genes/Iteration_",i, "/")
  system(paste0("mkdir ", output_link))
  algorithm_link <- "../programmes/algorithm_new/"
  input_data <- "Sce_Dataset4_Cells.RData"
  final_output <- "Clustering_Result_Genes_Stability.RData"
  
  cutoff <- (i*0.1)*nn
  sub_genes <- which(rank <= cutoff)
  
  sce <- sce_general[sub_genes,]
  save(sce, file=input_data)
        
  ## Save list of genes for reanalysis with cellranger
  gene<-rowData(sce)$id
  gene<-as.matrix(gene)
  colnames(gene)<-"Gene"
  write.csv(gene, file="cleaned_genes_half.csv", row.names = F, quote=F)
  
  third_param <- c("", "--barcode=cleaned_barcodes_dataset4.csv --genes=cleaned_genes_half.csv", "", " 20", "", "", "", "", " ~/10xGenomics/Clustering/Reference/BulkRNASeq_Marker_Genes.RData", "", "")
  algorithm <- c("ascend", "cellrangerS", "cidr", "CountClust", "RaceID", "RaceID2", "RCA", "sc3", "scran", "Seurat", "TSCAN")
  
  arguments <- list(output_link = output_link, algorithm_link = algorithm_link, algorithm = algorithm, input_data = input_data, third_param = third_param)
  
  call_algo <- function(x, arguments) {
    
    system.time(system(paste0( "Rscript --vanilla " , arguments$algorithm_link ,
                                arguments$algorithm[x], "_new_call.R ", arguments$input_data,
                                " ", arguments$output_link, arguments$algorithm[x], ".txt ", 
                                arguments$third_param[x])))
  }
  
  multicoreParam <- MulticoreParam(workers = 13)
  time_taken<-bplapply(1:11, call_algo, arguments=arguments,  BPPARAM = multicoreParam)
  
  names(time_taken) <- algorithm
  
  save(time_taken, file=paste0(output_link, "Time_Taken_Half_Iteration_", i, ".RData"))
  system(paste0("/bin/rm -f ", input_data))
  system(paste0("/bin/rm -f cleaned_genes_half.csv"))
}

res <- list(ascend=list(), cellranger=list(), cidr=list(), CountClust=list(), RaceID=list(), RaceID2=list(), RCA=list(), sc3=list(), scran=list(), Seurat=list(), TSCAN=list())

for(i in seq(1,5)){
    output_link <- paste0("Clustering_Results/Stability/Genes/Iteration_",i, "/")
    
    for (x in algorithm){
        try({par_res <- read.table(paste0(output_link, x, ".txt"), header=TRUE)
        par_res <- par_res[,dim(par_res)[2]]
        eval(parse(text=paste0("res$", x, "[[", i, "]]", "<- par_res")))})
    }
    
}

save(res, file=final_output)


