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
load("Clustering_Robustness_Genes.RData")

## Figure out ARI for each combination of clustering solutions in the same method
comp_cluster<-function(tmp){
  
  a<-length(tmp)-1
  b<-length(tmp)
  nc<-length(tmp[[1]])
  
  comp<-sapply(1:a, function(x) sapply((x+1):b, function(y) adjustedRandIndex(tmp[[x]], tmp[[y]])))
  
  return(unlist(comp))
}

res_comp<-matrix(NA, ncol=length(res), nrow=45)
colnames(res_comp)<-names(res)


for(i in 1:length(res)){
  
  tmp<-lapply(1:10, function(x) res[[i]][[x]])
  
  res_comp[,i]<-comp_cluster(tmp)
}

res_comp <- res_comp[,order(colnames(res_comp), decreasing=F)]
res_comp <- as.data.frame(res_comp)


res_comp <- melt(res_comp)
colnames(res_comp) <- c("Method", "Adjusted Rand Index")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(length(res))

gg1 <- ggplot(res_comp, aes(x=`Method`, y=`Adjusted Rand Index`)) + geom_boxplot(aes(color=Method))+ theme_minimal() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + annotate("text", label = "*", x = 5, y = 1, size = 8, colour = colors_gg[5])
gg1 <- ggdraw(gg1) + draw_plot_label("b")

ggsave(gg1, file="Figure_Gene_Stability.pdf", width=8, height=5)


