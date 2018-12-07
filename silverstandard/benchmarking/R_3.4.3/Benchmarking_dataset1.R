library(BiocParallel)

output_link <- "Clustering_Results/Dataset1/"
algorithm_link <- "../../programmes/algorithms/"
input_data <- "Sce_Dataset1.RData"
final_output <- "Clustering_Result_Dataset1.RData"


third_param <- c("", "", 8, "", "", "", "", "../prior_info/BulkRNASeq_Marker_Genes.RData", "", 8, "", "")
algorithm <- c("ascend", "cidr", "CountClust", "RaceID", "RaceID2", "RCA", "sc3", "scran", "Seurat", "SIMLR", "TSCAN", "cellranger")

arguments <- list(output_link = output_link, algorithm_link = algorithm_link, algorithm = algorithm, input_data = input_data,
                  third_param = third_param)

call_algo <- function(x, arguments) {
  system( paste0( "Rscript --vanilla " , arguments$algorithm_link , arguments$algorithm[x], "_call.R ", arguments$input_data,
                  " ", arguments$output_link, arguments$algorithm[x], ".txt", arguments$third_param[x]) ) 
}

multicoreParam <- MulticoreParam(workers = 12)
bplapply(1:12, call_algo, arguments=arguments,  BPPARAM = multicoreParam)

files <- list.files(output_link)

res <- read.table(paste0(output_link, files[1]), header = TRUE)
res <- res [, -dim(res)[2]]

for (i in files){
  tmp <- read.table(paste0(output_link, i), header = TRUE)
  eval(parse( text = paste0( "res$", colnames(tmp)[dim(tmp)[2]], "<- tmp[, dim(tmp)[2]]")))
}

save(res, file=final_output)

