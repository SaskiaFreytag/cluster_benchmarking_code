library(mclust)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(infotheo)
library(reshape2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


setwd("~/10xGenomics/Clustering/")

load("Clustering_Result_Dataset2b.RData")
res_2a <- res
colnames(res_2a)[77] <- "Cell Ranger"
colnames(res_2a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset2b.RData")
res_2a$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset3b.RData")
res_3a <- res
colnames(res_3a)[77] <- "Cell Ranger"
colnames(res_3a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset3b.RData")
res_3a$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset2.RData")
res_2 <- res

load("Assigned_Cell_Types_Dataset2.RData")
res_2$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType


load("Clustering_Result_Dataset3.RData")
res_3 <- res

load("Assigned_Cell_Types_Datset3.RData")
res_3$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType


res_2 <- list(`dataset 2`=res_2, `dataset 2a`=res_2a)
res_3 <- list(`dataset 3`=res_3, `dataset 3a`=res_3a)

## Figure out ARI of the solutions

programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat","TSCAN")         

index_2 <- lapply(res_2, function(x) pmatch(programs, colnames(x)))
index_3 <- lapply(res_3, function(x) pmatch(programs, colnames(x)))

## Adjusted Rand Index

comp_2 <- lapply(1:length(index_2), function(x) apply(res_2[[x]][, index_2[[x]]], 2, function(y) mclust::adjustedRandIndex(y, res_2[[x]]$Assigned_CellTypes)))
comp_3 <- lapply(1:length(index_3), function(x) apply(res_3[[x]][, index_3[[x]]], 2, function(y) mclust::adjustedRandIndex(y, res_3[[x]]$Assigned_CellTypes)))

comp_2 <- Reduce(rbind, comp_2)
rownames(comp_2) <- c("ARI_truth Evaluation 1", "ARI_truth Evaluation 2")
comp_3 <- Reduce(rbind, comp_3)
rownames(comp_3) <- c("ARI_truth Evaluation 1", "ARI_truth Evaluation 2")

comp_3 <- data.frame (t(comp_3), Dataset="Dataset 3/3a", Method = colnames(comp_3), check.names = FALSE)  
comp_2 <- data.frame (t(comp_2), Dataset="Dataset 2/2a", Method = colnames(comp_2), check.names = FALSE)  

comp <- rbind(comp_2, comp_3)

colours <- gg_color_hue(12)

a_plot <- ggplot(comp, aes(x=`ARI_truth Evaluation 1`, y=`ARI_truth Evaluation 2`)) + geom_point(aes(col=Method, shape=Dataset), size=2) + scale_color_manual(values=colours[-11]) +  geom_abline(slope=1, linetype="dotted") +
   guides(color=FALSE) + guides(shape=FALSE) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
(a_plot <- ggdraw(a_plot) + draw_plot_label("a"))

## Homogeneity

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

homo_2 <- lapply(1:length(index_2), function(x) apply(res_2[[x]][, index_2[[x]]], 2, function(y) homogeneity_score(res_2[[x]]$Assigned_CellTypes, y)))
homo_3 <- lapply(1:length(index_3), function(x) apply(res_3[[x]][, index_3[[x]]], 2, function(y) homogeneity_score(res_3[[x]]$Assigned_CellTypes, y)))

homo_2 <- Reduce(rbind, homo_2)
rownames(homo_2) <- c("Homogeneity Evaluation 1", "Homogeneity Evaluation 2")
homo_3 <- Reduce(rbind,homo_3)
rownames(homo_3) <- c("Homogeneity Evaluation 1", "Homogeneity Evaluation 2")

homo_3 <- data.frame (t(homo_3), Dataset="Dataset 3/3a", Method = colnames(homo_3), check.names = FALSE)  
homo_2 <- data.frame (t(homo_2), Dataset="Dataset 2/2a", Method = colnames(homo_2), check.names = FALSE)  

homo <- rbind(homo_2, homo_3)

b_plot <- ggplot(homo, aes(x=`Homogeneity Evaluation 1`, y=`Homogeneity Evaluation 2`)) + geom_point(aes(col=Method, shape=Dataset), size=2) + scale_color_manual(values=colours[-11]) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_abline(slope=1, linetype="dotted")
(b_plot <- ggdraw(b_plot) + draw_plot_label("b"))

gg <- grid.arrange(a_plot, b_plot, widths = 4:5, ncol=2)
ggsave("Figure_Compare_R.pdf", plot=gg, height=5, width=12)
