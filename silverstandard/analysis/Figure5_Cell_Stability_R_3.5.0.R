###########################################################################################
#                                                                                         #
#                                                                                         #
#                     Produce Plot of Gene Filtering Stability                            #
#                                                                                         #
#                                                                                         #
###########################################################################################


library(cowplot)

# Load results from stability analysis
load("Clustering_Result_Stability.RData")

load("Assigned_Cell_Types_Dataset5.RData")

# Find intersection of all subsampled datasets
load("Sce_Dataset5.RData")
set.seed(162622829)
subsampled <- lapply(1:5, function(x) sample(colData(sce)$barcode, 4000))

inter_subsampled <- Reduce(intersect, subsampled)
barcodes_indices <- lapply(subsampled, function(x) pmatch(inter_subsampled, x))

library(mclust)
library(NMI)


comp_cluster<-function(tmp){
  
  a<-length(tmp)-1
  b<-length(tmp)
  nc<-length(tmp[[1]])
  
  comp<-sapply(1:a, function(x) sapply((x+1):b, function(y) mclust::adjustedRandIndex(tmp[[x]], tmp[[y]])))
  
  return(unlist(comp))
}

res_comp<-matrix(NA, ncol=length(res), nrow=10)
colnames(res_comp)<-names(res)


for(i in 1:length(res)){
  
  tmp<-lapply(1:5, function(x) res[[i]][[x]][barcodes_indices[[x]]])
  
  res_comp[,i]<-comp_cluster(tmp)
}

res_comp <- res_comp[,order(colnames(res_comp), decreasing=F)]
res_comp <- as.data.frame(res_comp)

library(reshape2)

res_comp <- melt(res_comp)
colnames(res_comp) <- c("Method", "ARI_comp")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(12)[-11]

library(ggplot2)
gg1 <- ggplot(res_comp, aes(x=`Method`, y=`ARI_comp`)) + geom_boxplot(aes(color=Method)) + theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_manual(values=colors_gg) + guides(color=FALSE)
gg1 <- ggdraw(gg1) +  draw_plot_label("a") 

## compare to truth

indices <- lapply(subsampled, function(x) match(x, sce$barcode))

comp_truth<-function(tmp, index){
  
  b<-length(tmp)
  
  comp<-sapply(1:b, function(y)
  {assigned_cell_types <- assigned_cell_types$Assigned_CellType[index[[y]]]
    mclust::adjustedRandIndex(tmp[[y]],  assigned_cell_types)})
  
  return(unlist(comp))
}

  
res_truth<-lapply(1:length(res), function(x) comp_truth(res[[x]], indices))
res_truth<-Reduce(rbind, res_truth)
res_truth <- t(res_truth)
colnames(res_truth)<-names(res)  
res_truth <- res_truth[, order(colnames(res_truth), decreasing=F)]

res_truth <- melt(res_truth)
colnames(res_truth) <- c("Number", "Method", "ARI_truth")


gg2 <- ggplot(res_truth, aes(x=`Method`, y=`ARI_truth`)) + geom_boxplot(aes(color=Method)) + theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_color_manual(values=colors_gg)
gg2 <- ggdraw(gg2) +  draw_plot_label("b") 


gg <- grid.arrange(gg1, gg2, ncol=2, widths=4:5)

ggsave(gg, file="Figure_Cell_Stability.pdf", width=12, height=5)
