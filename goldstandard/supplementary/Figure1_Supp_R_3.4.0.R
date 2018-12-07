###########################################################################################
#                                                                                         #
#                                                                                         #
#                        MAKE TSNE PLOTS FOR SUPPLEMENTARY                                #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(scater)
library(data.table)
library(ggplot2)
library(Rtsne)
library(gridExtra)
library(cowplot)


## Read in 10x Genomics data
load("Sce_CellRanger.RData")

## Read in clustering results
load("Clustering_Result.RData")
colData(sce) <- res
colData(sce)$barcode <- gsub("-1", "", colData(sce)$barcode) # fix barcodes

## Read in demuxlet calls
demuxlet <- fread("demuxlet_calls.best")

## Match barcodes and add to objects
colData(sce)$demuxlet <- demuxlet$SNG.1ST[pmatch(colData(sce)$barcode, demuxlet$BARCODE)]

## Normalize and drop genes 
sce <- normalize(sce)
drop_genes <- apply(exprs(sce), 1, function(x) {var(x) == 0})
sce <- sce[!drop_genes, ]

## Plotting t-SNE
set.seed(1)
tsne <- Rtsne(t(exprs(sce)), perplexity = 12)

dr <- data.frame(`t-SNE Dim 1`=tsne$Y[,1], `t-SNE Dim 2`=tsne$Y[,2], `RaceID2`=as.factor(colData(sce)$RaceID2), `Seurat`=colData(sce)$Seurat, `SC3`=as.factor(colData(sce)$SC3), `demuxlet`=colData(sce)$demuxlet, check.names=FALSE)

g1 <- ggplot(dr, aes(y=`t-SNE Dim 2`, x=`t-SNE Dim 1`)) + geom_point(aes(shape=demuxlet, col=Seurat)) + theme_minimal()
g2 <- ggplot(dr, aes(y=`t-SNE Dim 2`, x=`t-SNE Dim 1`)) + geom_point(aes(shape=demuxlet, col=SC3)) + theme_minimal ()
g3 <- ggplot(dr, aes(y=`t-SNE Dim 2`, x=`t-SNE Dim 1`)) + geom_point(aes(shape=demuxlet, col=RaceID2)) +theme_minimal ()


g1 <- ggdraw(g1) + draw_figure_label("c")
g2 <- ggdraw(g2) + draw_figure_label("b")
g3 <- ggdraw(g3) + draw_figure_label("a")

ggsave(g1, file="TSNE_Seurat.pdf")
ggsave(g2, file="TSNE_SC3.pdf")
ggsave(g3, file="TSNE_RaceID2.pdf")

