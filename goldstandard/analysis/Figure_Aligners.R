###########################################################################################
#                                                                                         #
#                                                                                         #
#                       ASSESSING WITHIN CLUSTERING RESULTS                              #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(mclust)
library(NMI)
library(ggplot2)
library(reshape2)


## Load results cellranger
load("Clustering_Result_CellRanger.RData")
res_cellranger <- res

## Load results scpipe subread
load("Clustering_Result_ScPipe_SubRead.RData")
res_scpipe <- res

## Load results scpipe star
load("Clustering_Result_ScPipe_STAR.RData")
res_star <- res

programs <- c("CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran", "Seurat", "TSCAN")      


## Subset to barcodes that are common and make sure they are in the same order

res_cellranger$barcode <- gsub("-1", "", as.character(res_cellranger$barcode))
all_barcodes <- intersect(res_star$barcode, intersect(as.character(res_scpipe$barcode), as.character(res_cellranger$barcode)))

res_cellranger<-res_cellranger[pmatch(all_barcodes, as.character(res_cellranger$barcode)),]
res_scpipe<-res_scpipe[pmatch(all_barcodes, as.character(res_scpipe$barcode)),]
res_star<-res_star[pmatch(all_barcodes, as.character(res_star$barcode)),]

## Function to compare clustering results
comp_clusters<-function(res_clusters1, res_clusters2){
  
  a<-dim(res_clusters1)[2]
  nc<-dim(res_clusters1)[1]
  comp<-sapply(1:a, function(x)  adjustedRandIndex(res_clusters1[, x], res_clusters2[, x]))
  
  return(comp) 
}

res_comp <- cbind(comp_clusters(res_cellranger[,programs], res_scpipe[,programs]),
      comp_clusters(res_cellranger[,programs], res_star[,programs]),
      comp_clusters(res_star[,programs], res_scpipe[,programs]))

colnames(res_comp) <- c("Cell Ranger- ScPipe Subread", "Cell Ranger - ScPipe STAR", "ScPipe Subread - ScPipe Star")
rownames(res_comp) <- programs

res_comp <- melt(res_comp)
colnames(res_comp) <- c("Method", "Comparison","Adjusted Rand Index")

res_comp <- as.data.frame(res_comp)
res_comp$`Adjusted Rand Index`<-as.numeric(as.character(res_comp$`Adjusted Rand Index`))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(12)[c(3,4,5,6,7,8,9,10,12)]

gg1 <- ggplot(res_comp, aes(x=Method, y=`Adjusted Rand Index`)) + geom_col(aes(fill=Method)) + theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~ Comparison)+scale_fill_manual(values=colors_gg)
gg1

ggsave(gg1, file="Figure_Aligners.pdf", width=12, height=5)


