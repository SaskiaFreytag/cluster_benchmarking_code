###########################################################################################
#                                                                                         #
#                                                                                         #
#                     Produce Plot of Gene Filtering Stability                            #
#                                                                                         #
#                                                                                         #
###########################################################################################


library(cowplot)

# Load results from stability analysis
load("Clustering_Result_Genes_Stability.RData")
res_gen <- res

res_gen [[2]] <- res_gen[[12]]
res_gen <- res_gen[-12]

load("Assigned_Cell_Types_Dataset4.RData")

# Find intersection of all subsampled datasets
load("Clustering_Result_Dataset4.RData")

library(mclust)
library(NMI)

colnames(res)[76:86] <- c("ascend", "CellRanger", "CIDR", "countClust", "RaceID2", "RaceID", "RCA", "SC3", "scran",       "Seurat","TSCAN")

algorithms <- tolower(names(res_gen))
index <- match(algorithms, tolower(colnames(res)))

comp_cluster<-function(tmp, tmp1){
  
  comp<-mclust::adjustedRandIndex(tmp, tmp1)
  
  return(comp)
}

res_comp<-matrix(NA, ncol=length(res_gen), nrow=5)
colnames(res_comp)<-names(res_gen)


for(i in 1:length(res_gen)){
  
  tmp1 <- res[, index[i]]
  res_tmp <- sapply(res_gen[[i]], function(x) try(comp_cluster(x, tmp1)))
  try(res_comp[,i]<-res_tmp)
}

res_comp <- as.data.frame(res_comp)
ind <- apply(res_comp, 1, function(x) grep("Error", x))
ind <- sapply(1:3, function(x) cbind(x, ind[[x]]))
ind <- Reduce(rbind, ind)

res_comp[ind] <- NA
rownames(res_comp) <- c("10%", "20%", "30%", "40%", "50%")
res_comp <- apply(res_comp, 2, as.numeric)

library(reshape2)

res_comp <- melt(res_comp)
colnames(res_comp) <- c("Percentage", "Method", "ARI_comp")
res_comp$Percentage <- res_comp$Percentage*10
res_comp$Percentage <- as.factor(res_comp$Percentage)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(12)[-11]

library(ggplot2)
gg1 <- ggplot(res_comp, aes(x=`Percentage`, y=`ARI_comp`, group=Method)) + geom_point(aes(color=Method)) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_manual(values=colors_gg) +   geom_line(aes(color=Method)) + guides(color=FALSE)
gg1 <- ggdraw(gg1) +  draw_plot_label("a") 

## compare to truth


comp_truth<-function(tmp, assigned_cell_types){
  
  comp <- mclust::adjustedRandIndex(tmp,  assigned_cell_types)
  
  return(comp)
}

res_truth <- matrix(NA, ncol=length(res_gen), nrow=5)
colnames(res_truth)<-names(res_gen)

for(i in 1:length(res_gen)){
  
  res_tmp <- sapply(res_gen[[i]], function(x) try(comp_truth(x, assigned_cell_types$Assigned_CellType)))
  try(res_truth[,i]<-res_tmp)
}
  
ind <- apply(res_truth, 1, function(x) grep("Error", x))
ind <- sapply(1:3, function(x) cbind(x, ind[[x]]))
ind <- Reduce(rbind, ind)

res_truth[ind] <- NA
res_truth <- apply(res_truth, 2, as.numeric)

library(reshape2)

res_truth <- melt(res_truth)
colnames(res_truth) <- c("Percentage", "Method", "ARI_truth")
res_truth$Percentage <- res_truth$Percentage*10
res_truth$Percentage <- as.factor(res_truth$Percentage)


gg2 <- ggplot(res_truth, aes(x=`Percentage`, y=`ARI_truth`, group=Method)) + geom_point(aes(color=Method)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_color_manual(values=colors_gg) + geom_line(aes(color=Method)) 
gg2 <- ggdraw(gg2) +  draw_plot_label("b") 

gg <- grid.arrange(gg1, gg2, ncol=2, widths=4:5)

ggsave(gg, file="Figure_Gene_Stability.pdf", width=12, height=5)
