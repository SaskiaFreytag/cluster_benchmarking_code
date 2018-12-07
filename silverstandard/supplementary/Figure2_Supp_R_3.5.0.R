###########################################################################################
#                                                                                         #
#                                                                                         #
#                             Produce Figure 2 Paper                                      #
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


### Read in clustering results
load("Clustering_Result_Dataset2a.RData")
res_dataset2a <- res
colnames(res_dataset2a)[77] <- "Cell Ranger"

load("Assigned_Cell_Types_Info_Dataset2a.RData")
res_dataset2a <- cbind(res_dataset2a, csl_id1)
colnames(res_dataset2a)[(dim(res_dataset2a)[2]-1):(dim(res_dataset2a)[2])] <- c("Assigned Cell Type", "Correlation")

#

load("Clustering_Result_Dataset3a.RData")
res_dataset3a <- res
colnames(res_dataset3a)[77] <- "Cell Ranger"

load("Assigned_Cell_Types_Info_Dataset3a.RData")
res_dataset3a <- cbind(res_dataset3a, csl_id1)
colnames(res_dataset3a)[(dim(res_dataset3a)[2]-1):(dim(res_dataset3a)[2])] <- c("Assigned Cell Type", "Correlation")

#

load("Clustering_Result_Dataset4.RData")
res_dataset4 <- res
colnames(res_dataset4)[77] <- "Cell Ranger"

load("Assigned_Cell_Types_Info_Dataset4.RData")
res_dataset4 <- cbind(res_dataset4, csl_id1)
colnames(res_dataset4)[(dim(res_dataset4)[2]-1):(dim(res_dataset4)[2])] <- c("Assigned Cell Type", "Correlation")

#

load("Clustering_Result_Dataset5.RData")
res_dataset5 <- res
colnames(res_dataset5)[77] <- "Cell Ranger"

load("Assigned_Cell_Types_Info_Dataset5.RData")
res_dataset5 <- cbind(res_dataset5, csl_id1)
colnames(res_dataset5)[(dim(res_dataset5)[2]-1):(dim(res_dataset5)[2])] <- c("Assigned Cell Type", "Correlation")

res_all <- list(`dataset 2a`=res_dataset2a,  `dataset 3a`=res_dataset3a , `dataset 4`= res_dataset4, `dataset 5`= res_dataset5)

# Find index of clustering results
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat"
               ,"TSCAN")         
index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))


rr <- list(seq(0.1, 0.3, 0.025), seq(0.5, 0.70, 0.025), seq(0.15, 0.35, 0.025), seq(0.45, 0.65, 0.025))
rr_comp <- list() 
for (i in 1:9){
  print(i)
  res_tmp <- lapply(1:length(res_all), function(x) 
    res_all[[x]][which(as.numeric(as.character(res_all[[x]]$Correlation))>=rr[[x]][i]),])
  rr_comp[[i]] <- lapply(1:length(index_all), function(x) apply(res_tmp[[x]][, index_all[[x]]], 2, function(y) 
    adjustedRandIndex(y, res_tmp[[x]]$`Assigned Cell Type`)))
}

rr_comp <-rbind(cbind(Reduce(c, lapply(rr_comp, function(x) x[[1]])), rep(rr[[1]], each=length(programs)), "dataset 2a"),
      cbind(Reduce(c, lapply(rr_comp, function(x) x[[2]])), rep(rr[[2]], each=length(programs)), "dataset 3a"),
      cbind(Reduce(c, lapply(rr_comp, function(x) x[[3]])), rep(rr[[3]], each=length(programs)), "dataset 4"),
      cbind(Reduce(c, lapply(rr_comp, function(x) x[[4]])), rep(rr[[4]], each=length(programs)), "dataset 5"))

rr_comp <- as.data.frame(rr_comp)
rr_comp$Method <- programs
rownames(rr_comp) <- NULL
colnames(rr_comp) <- c("ARI_truth", "Minimum Correlation", "Dataset", "Method")
rr_comp$Method[which(rr_comp$Method=="Cell.Ranger")] <- "Cell Ranger"
rr_comp$ARI_truth <- as.numeric(as.character(rr_comp$ARI_trut))
rr_comp$`Minimum Correlation` <- as.numeric(as.character(rr_comp$`Minimum Correlation`))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(12)[-11]

g2 <- ggplot(rr_comp, aes(x=`Minimum Correlation`, y=`ARI_truth`, color=Method, group=Method)) + theme_bw(base_size = 11, base_family = "")  + geom_line(aes(color=Method)) + facet_grid(.~Dataset, scales = "free") + scale_color_manual(values=colors_gg)


ggsave("Figure_Supp5bb.pdf", plot=g2, height=5, width=12)
