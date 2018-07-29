###########################################################################################
#                                                                                         #
#                                                                                         #
#                     Produce Plot of Gene Filtering Stability                            #
#                                                                                         #
#                                                                                         #
###########################################################################################


library(cowplot)
library(mclust)
library(NMI)
library(ggplot2)
library(reshape2)

# Load results from stability analysis
load("Clustering_Robustness_Cells.RData")


# Find intersection of all subsampled datasets
load("Sce_Dataset3.RData")
set.seed(162622829)
subsampled <- lapply(1:5, function(x) sample(colData(sce)$barcode, 3000))

inter_subsampled <- Reduce(intersect, subsampled)
barcodes_indices <- lapply(subsampled, function(x) pmatch(inter_subsampled, x))

## Figure out ARI for each combination of clustering solutions in the same method
comp_cluster<-function(tmp){
  
  a<-length(tmp)-1
  b<-length(tmp)
  nc<-length(tmp[[1]])
  
  comp<-sapply(1:a, function(x) sapply((x+1):b, function(y) adjustedRandIndex(tmp[[x]], tmp[[y]])))
  
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


res_comp <- melt(res_comp)
colnames(res_comp) <- c("Method", "ARI_comp")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(length(res))

gg1 <- ggplot(res_comp, aes(x=`Method`, y=`ARI_comp`)) + geom_boxplot(aes(color=Method))+ theme_minimal() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
gg1 <- ggdraw(gg1) +  draw_plot_label("a") 

ggsave(gg1, file="Figure_Cell_Stability.pdf", width=8, height=5)
