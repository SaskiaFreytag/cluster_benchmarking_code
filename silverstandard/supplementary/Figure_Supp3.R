library(cowplot)
library(scater)
library(mclust)
library(NMI)
library(reshape2)
library(dplyr)
library(ggplot2)


comp_cluster<-function(tmp){
  
  a<-length(tmp)-1
  b<-length(tmp)
  nc<-length(tmp[[1]])
  
  comp<-sapply(1:a, function(x) sapply((x+1):b, function(y) adjustedRandIndex(tmp[[x]], tmp[[y]])))
  
  return(unlist(comp))
}

## Load results for gene stability analysis
load("Clustering_Robustness_Genes.RData")

res_comp<-matrix(NA, ncol=length(res), nrow=45)
colnames(res_comp)<-names(res)

for(i in 1:length(res)){
  
  tmp<-lapply(1:10, function(x) res[[i]][[x]])
  
  res_comp_gene[,i]<-comp_cluster(tmp)
}

res_comp_gene <- res_comp_gene[,order(colnames(res_comp_gene), decreasing=F)]
res_comp_gene <- as.data.frame(res_comp_gene)

# Load results from cell stability analysis
load("Clustering_Robustness_Cells.RData")

# Find intersection of all subsampled datasets
load("Sce_Dataset3.RData")
set.seed(162622829)
subsampled <- lapply(1:5, function(x) sample(colData(sce)$barcode, 3000))

inter_subsampled <- Reduce(intersect, subsampled)
barcodes_indices <- lapply(subsampled, function(x) pmatch(inter_subsampled, x))

res_comp_cell<-matrix(NA, ncol=length(res), nrow=10)
colnames(res_comp_cell)<-names(res)


for(i in 1:length(res)){
    
    tmp<-lapply(1:5, function(x) res[[i]][[x]][barcodes_indices[[x]]])
    
    res_comp[,i]<-comp_cluster(tmp)
}

res_comp_cell <- res_comp[,order(colnames(res_comp), decreasing=F)]
res_comp_cell <- as.data.frame(res_comp_cell)

## Coerce data, work out median and join

res_comp_cell <- melt(res_comp_cell)
colnames(res_comp_cell) <- c("Method", "ARI_comp")

res_comp_gene <- melt(res_comp_gene)
colnames(res_comp_gene) <- c("Method", "ARI_comp")

res_comp_gene <- res_comp_gene %>% group_by(Method) %>% summarize(`Median ARI_comp Gene`=median(ARI_comp)) %>% ungroup()
res_comp_cell <- res_comp_cell %>% group_by(Method) %>% summarize(`Median ARI_comp Cell`=median(ARI_comp)) %>% ungroup()
res_comp_gene$Method <- res_comp_cell$Method

res <- left_join(res_comp_cell, res_comp_gene)
res$Method <- as.character(res$Method)
res$Method[res$Method == "RaceID"] <-"RaceID*"


# Plot

gg <- ggplot(res, aes(`Median ARI_comp Cell`, `Median ARI_comp Gene`)) +  
  geom_point(aes(`Median ARI_comp Cell`, `Median ARI_comp Gene`, color = Method)) + 
  geom_text_repel(aes(`Median ARI_comp Cell`, `Median ARI_comp Gene`, label = Method, color = Method)) +
  theme_minimal() +  guides(color = FALSE)

ggsave(gg, file="Supplementary_Stability.pdf")
