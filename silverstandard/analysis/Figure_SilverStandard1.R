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

## Read in clustering results
load("Clustering_Result_Dataset1.RData")
res_csl <- res
load("Assigned_Cell_Types_Dataset1.RData")
res_csl$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset3.RData")
res_pbmc4 <- res
load("Assigned_Cell_Types_Dataset3.RData")
res_pbmc4$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType


load("Clustering_Result_Dataset2.RData")
res_pbmc3 <- res
load("Assigned_Cell_Types_Dataset2.RData")
res_pbmc3$Assigned_CellTypes <- assigned_cell_types$Assigned_CellType

res_all <- list(`dataset 1`=res_csl,  `dataset 2`=res_pbmc3, `dataset 3`= res_pbmc4)

# Find index of clustering results
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat", "SIMLR"
               ,"TSCAN")         
index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))

## Figure out homogeneity with regards to truth

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

homo_all <- melt(homo_all)

colnames(homo_all) <- c("Method", "Dataset", "Homogeneity")
homo_all$Dataset <- paste("dataset", homo_all$Dataset)


## Plot homogeneity
g1 <- ggplot(homo_all, aes(x=Method, y=Homogeneity)) + geom_col(aes(fill=Method, group=Dataset), position = "dodge") +
geom_col(aes(alpha=Dataset), position = "dodge", fill="black") +
scale_alpha_manual(values=c(0.05, 0.2, 0.5)) +
theme_minimal() + guides(fill=FALSE) +  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g1 <- ggdraw(g1)+ draw_plot_label("b")

## Figure out ARI and number of clusters with respect to the truth
comp <- sapply(1:length(index_all), function(x) apply(res_all[[x]][, index_all[[x]]], 2, function(y) adjustedRandIndex(y, res_all[[x]]$Assigned_CellTypes)))
cluster_number <- sapply(1:length(index_all), function(x) apply(res_all[[x]][, index_all[[x]]], 2, function(y) length(unique(y))))

colnames(comp) <- names(res_all)
colnames(cluster_number) <- names(res_all)

comp <- melt(comp)
cluster_number <- melt(cluster_number)

colnames(comp) <- c("Method", "Dataset", "Adjusted Rand Index")
colnames(cluster_number) <- c("Method", "Dataset", "Number of Clusters")

comp <- merge(comp, cluster_number)

## Plot ARI
g2 <- ggplot(comp, aes(x=`Number of Clusters`, y=`Adjusted Rand Index`, color=Method)) + theme_bw(base_size = 11, base_family = "") + guides(color=FALSE) + expand_limits(x=c(-5, 120), y=c(0,1)) + scale_x_log10(breaks=c(5, 10, 20, 40, 80, 160)) + geom_point(aes(`Number of Clusters`, `Adjusted Rand Index`, color=Method, shape=`Dataset`)) + geom_text_repel(aes(`Number of Clusters`, `Adjusted Rand Index`, label = Method, color=Method)) + geom_vline(xintercept=11, linetype = 2)
g2 <- ggdraw(g2) + draw_plot_label("a")

gg <- grid.arrange(g2, g1, nrow=1)

ggsave("Figure_SilverStandard1.pdf", plot=gg, height=5, width=12)
