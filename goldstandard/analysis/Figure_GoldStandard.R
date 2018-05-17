###########################################################################################
#                                                                                         #
#                                                                                         #
#                             Produce Figure 1 Paper                                      #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(scater)
library(data.table)
library(ggplot2)
library(mclust)
library(ggrepel)
library(gridExtra)
library(cowplot)


## Read in clustering results
load("Clustering_Result_CellRanger.RData")
res_cellranger <- res
res_cellranger$barcode <- gsub("-1", "", res_cellranger$barcode)

## Load demuxlet results
demuxlet <- fread("demuxlet_calls.best")

## Add demuxlet results
res_cellranger$demuxlet <- demuxlet$SNG.1ST[pmatch(res_cellranger$barcode, demuxlet$BARCODE)]

# Find index of clustering results
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran", "SIMLR"
               , "Seurat","TSCAN")         
index_cellranger <- pmatch(programs, colnames(res_cellranger))
index_cellranger<- index_cellranger[!is.na(index_cellranger)]


##  Figure out the inaccuracy of each method
inaccuracy <- function(tab) {
  all_inacc <- apply(tab, 1, function(x) {
    max_val <- which.max(x)
    sum(x[-max_val])/sum(x)}
  )
  sum(all_inacc)/dim(tab)[1]
}

tab_cellranger <- lapply(index_cellranger, function(x) table(res_cellranger[, c(x, which(colnames(res_cellranger) == "demuxlet"))]))
names(tab_cellranger) <- colnames(res_cellranger)[index_cellranger]
accuracy_cellranger <- cbind(programs, sapply(tab_cellranger, function(x) 1-inaccuracy(x)))

colnames(accuracy_cellranger) <- c("Method", "Accuracy") 
accuracy_cellranger <- as.data.frame(accuracy_cellranger)

accuracy_cellranger$Accuracy<-as.numeric(as.character(accuracy_cellranger$Accuracy))

## Plot accuracy
g1 <- ggplot(accuracy_cellranger, aes(x=Method, y=Accuracy)) + geom_col(aes(fill=Method)) + 
  theme_minimal() + guides(fill=FALSE) +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
g1 <- ggdraw(g1)+ draw_plot_label("b") 

##  Figure out the ARI given demuxlet of each method and compare to number of clusters
comp<-sapply(index_cellranger, function(x) adjustedRandIndex(res_cellranger[, x], res_cellranger$demuxlet))
cluster_number <- sapply(index_cellranger, function(x) length(unique(res_cellranger[, x])))
comp <- data.frame(Method = programs, `Adjusted Rand Index` = comp, `Number of Clusters` = cluster_number, check.names=FALSE)

# Plot ARI versus number of clusters
g2 <- ggplot(comp, aes(x=`Number of Clusters`, y=`Adjusted Rand Index`, color=Method)) + theme_bw(base_size = 11, base_family = "") + guides(color=FALSE) + expand_limits(x=c(-5, 120), y=c(0,1)) + scale_x_log10(breaks=c(5, 10, 20, 40, 80, 160)) + geom_point(aes(`Number of Clusters`, `Adjusted Rand Index`, color=Method)) + geom_text_repel(aes(`Number of Clusters`, `Adjusted Rand Index`, label = Method, color=Method)) + geom_vline(xintercept=3, linetype = 2)   
  
g2 <-  ggdraw(g2) + draw_plot_label("a")

gg <- grid.arrange(g2, g1, nrow=1)

ggsave("Figure_GoldStandard.pdf", plot=gg, width=12, height=5)




