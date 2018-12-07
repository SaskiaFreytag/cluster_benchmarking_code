###########################################################################################
#                                                                                         #
#                                                                                         #
#                             Produce Supplementary Figure 2                              #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(scater)
library(data.table)
library(ggplot2)
library(mclust)
library(ggrepel)
library(gridExtra)
library(reshape2)
library(cowplot)

## Read in clustering results
load("Clustering_Result_Dataset1.RData")
res_1 <- res

load("Assigned_Cell_Types_Info_Dataset1.RData")
res_1 <- cbind(res_csl, csl_id1)
colnames(res_1)[(dim(res_1)[2]-1):(dim(res_1)[2])] <- c("Assigned Cell Type", "Correlation")

## Read in clustering results
load("Clustering_Result_Dataset2.RData")
res_2 <- res

load("Assigned_Cell_Types_Info_Dataset2.RData")
res_2 <- cbind(res_3, csl_id1)
colnames(res_2)[(dim(res_2)[2]-1):(dim(res_2)[2])] <- c("Assigned Cell Type", "Correlation")

## Read in clustering results
load("Clustering_Result_Dataset3.RData")
res_3 <- res

load("Assigned_Cell_Types_Info_Dataset3.RData")
res_3 <- cbind(res_3, csl_id1)
colnames(res_3)[(dim(res_3)[2]-1):(dim(res_3)[2])] <- c("Assigned Cell Type", "Correlation")

res_all <- list(`dataset 1`=res_1 , `dataset 2`=res_2, `dataset 3`= res_3)

# Find index of clustering results
programs <- c( "ascend", "Cell.Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat", "SIMLR"
               ,"TSCAN")         
index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))


rr <- list(seq(0.5, 0.70, 0.025),seq(0.1, 0.3, 0.025), seq(0.55, 0.75, 0.025))
rr_comp <- list() 
for (i in 1:9){
  print(i)
  res_tmp <- lapply(1:length(res_all), function(x) 
    res_all[[x]][which(as.numeric(as.character(res_all[[x]]$Correlation))>=rr[[x]][i]),])
  rr_comp[[i]] <- lapply(1:length(index_all), function(x) apply(res_tmp[[x]][, index_all[[x]]], 2, function(y) 
    adjustedRandIndex(y, res_tmp[[x]]$`Assigned Cell Type`)))
}

rr_comp <-rbind(cbind(Reduce(c, lapply(rr_comp, function(x) x[[1]])), rep(rr[[1]], each=length(programs)), "dataset 1"),
      cbind(Reduce(c, lapply(rr_comp, function(x) x[[2]])), rep(rr[[2]], each=length(programs)), "dataset 2"),
      cbind(Reduce(c, lapply(rr_comp, function(x) x[[3]])), rep(rr[[3]], each=length(programs)), "dataset 3"))

rr_comp <- as.data.frame(rr_comp)
rr_comp$Method <- rownames(rr_comp)
rownames(rr_comp) <- NULL
colnames(rr_comp) <- c("ARI_truth", "Minimum Correlation", "Dataset", "Method")
rr_comp$Method[which(rr_comp$Method=="Cell.Ranger")] <- "Cell Ranger"
rr_comp$ARI_truth <- as.numeric(as.character(rr_comp$ARI_trut))
rr_comp$`Minimum Correlation` <- as.numeric(as.character(rr_comp$`Minimum Correlation`))

g2 <- ggplot(rr_comp, aes(x=`Minimum Correlation`, y=`ARI_truth`, color=Method)) + theme_bw(base_size = 11, base_family = "")  + geom_line(aes(`Minimum Correlation`, `ARI_truth`, color=Method)) + facet_grid(.~Dataset, scales = "free")


ggsave("Figure_Supp5.pdf", plot=g2, height=5, width=12)
