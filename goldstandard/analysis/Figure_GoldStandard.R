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

## Figure out homogeneity for each method given the demuxlet results
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

homo_all <- apply(res_cellranger[, index_cellranger], 2, function(y) homogeneity_score(res_cellranger$demuxlet, y))
homo_all <- melt(homo_all)
homo_all$Method <- rownames(homo_all)

colnames(homo_all) <- c("Homogeneity", "Method")

## Plot homogeneity
g1 <- ggplot(homo_all, aes(x=Method, y=Homogeneity)) + geom_col(aes(fill=Method)) +
theme_minimal() + guides(fill=FALSE) +  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g1 <- ggdraw(g1)+ draw_plot_label("b")


##  Figure out the ARI given demuxlet of each method and compare to number of clusters

comp<-sapply(index_cellranger, function(x) adjustedRandIndex(res_cellranger[, x], res_cellranger$demuxlet))
cluster_number <- sapply(index_cellranger, function(x) length(unique(res_cellranger[, x])))
comp <- data.frame(Method = programs, `ARI_truth` = comp, `Number of Clusters` = cluster_number, check.names=FALSE)


g2 <- ggplot(comp, aes(x=`Number of Clusters`, y=`ARI_truth`, color=Method)) +
    theme_bw(base_size = 11, base_family = "") + guides(color=FALSE) +
    expand_limits(x=c(-5, 120), y=c(0,1)) +
    scale_x_log10(breaks=c(5, 10, 20, 40, 80, 160)) +
    geom_point(aes(`Number of Clusters`, `ARI_truth`, color=Method)) +
    geom_text_repel(aes(`Number of Clusters`, `ARI_truth`, label = Method, color=Method)) +
    geom_vline(xintercept=3, linetype = 2)

g2 <-  ggdraw(g2) + draw_plot_label("a")

gg <- grid.arrange(g2, g1, nrow=1)

ggsave("Figure_GoldStandard.pdf", plot=gg, width=12, height=5)
