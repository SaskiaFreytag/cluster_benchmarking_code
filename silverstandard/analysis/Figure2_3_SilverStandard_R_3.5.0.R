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
library(infotheo)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colours <- gg_color_hue(12)[-11]

## Read in clustering results
load("Clustering_Result_Dataset2b.RData")
res_dataset2a <- res
colnames(res_dataset2a)[77] <- "Cell Ranger"
colnames(res_dataset2a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset1a.RData")
res_dataset2a$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

#

load("Clustering_Result_Dataset3a.RData")
res_dataset3a <- reb
colnames(res_dataset3a)[77] <- "Cell Ranger"
colnames(res_dataset3a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset3b.RData")
res_dataset3a$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

#

load("Clustering_Result_Dataset4.RData")
res_dataset4 <- res
colnames(res_dataset4)[77] <- "Cell Ranger"
colnames(res_dataset4)[81] <- "RaceID"


load("Assigned_Cell_Types_Dataset4.RData")
res_dataset4$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType


#

load("Clustering_Result_Dataset5.RData")
res_dataset5 <- res
colnames(res_dataset5)[77] <- "Cell Ranger"
colnames(res_dataset5)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset5.RData")
res_dataset5$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

#

res_all <- list(`dataset 2a`=res_dataset2a,  `dataset 3a`=res_dataset3a , `dataset 4`= res_dataset4, `dataset 5`= res_dataset5)

# Find index of clustering results
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat"
               ,"TSCAN")         
index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))


homogeneity_score <- function(a, b) {
  mi <- mutinformation(a, b)
  entropy.a <- entropy(a)
  entropy.b <- entropy(b)
  if (entropy.a == 0.0) {
    homogeneity <- 1.0
  } else {
    homogeneity <- mi / entropy.a
  }
  homogeneity
}

homo_all <- sapply(1:length(index_all), function(x) apply(res_all[[x]][, index_all[[x]]], 2, function(y) homogeneity_score(res_all[[x]]$Assigned_CellTypes, y)))
colnames(homo_all) <- names(res_all)
homo_all <- melt(homo_all)

colnames(homo_all) <- c("Method", "Dataset", "Homogeneity") 


(g1 <- ggplot(homo_all, aes(x=Method, y=Homogeneity)) + 
  geom_col(aes(fill=Method, group=Dataset), position = "dodge") + 
  theme_minimal() + 
  guides(fill=FALSE) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=colours) +
  facet_grid(~Dataset))

(g1 <- ggdraw(g1)+ draw_plot_label("b"))

##  Figure out the ARI given labeling

comp <- sapply(1:length(index_all), function(x) apply(res_all[[x]][, index_all[[x]]], 2, function(y) mclust::adjustedRandIndex(y, res_all[[x]]$Assigned_CellTypes)))
cluster_number <- sapply(1:length(index_all), function(x) apply(res_all[[x]][, index_all[[x]]], 2, function(y) length(unique(y))))

colnames(comp) <- names(res_all)
colnames(cluster_number) <- names(res_all)

comp <- melt(comp)
cluster_number <- melt(cluster_number)

colnames(comp) <- c("Method", "Dataset", "ARI_truth")
colnames(cluster_number) <- c("Method", "Dataset", "Number of Clusters")

comp2 <- merge(comp, cluster_number, by=c("Method", "Dataset"))


g2 <- ggplot(comp2, aes(x=`Number of Clusters`, y=`ARI_truth`, color=Method)) + theme_bw(base_size = 11, base_family = "") + guides(color=FALSE) + 
  expand_limits(x=c(-5, 120), y=c(0,1)) + scale_x_log10(breaks=c(5, 10, 20, 40, 80, 160)) + 
  geom_point(aes(`Number of Clusters`, `ARI_truth`, color=Method, shape=`Dataset`)) + 
  geom_text_repel(aes(`Number of Clusters`, `ARI_truth`, label = Method)) + 
  scale_color_manual(values=colours) +
  geom_vline(xintercept=11, linetype = 2)
(g2 <- ggdraw(g2) + draw_plot_label("b"))



ggsave("Figure_SilverStandard_Homo_2.pdf", plot=g1, height=5, width=12)
ggsave("Figure_SilverStandard_ARI_2.pdf", plot=g2, height=5, width=6)


