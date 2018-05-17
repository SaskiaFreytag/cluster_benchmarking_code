###########################################################################################
#                                                                                         #
#                                                                                         #
#                             SUPPLEMENTARY TABLE 1                                       #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(xtable)

# Set working directory
setwd("~/10xGenomics/Clustering/")

load("PBMC3K/Assigned_Cell_Types.RData")
ass_pbmc3k <- assigned_cell_types

load("PBMC4K/Assigned_Cell_Types.RData")
ass_pbmc4k <- assigned_cell_types

load("CSL/Assigned_Cell_Types.RData")
ass_csl <- assigned_cell_types

res_all <- sapply(list(`dataset 1`=ass_csl$Assigned_CellType, `dataset 2`=ass_pbmc3k$Assigned_CellType, `dataset 3`=ass_pbmc4k$Assigned_CellType), function(x) table(x)/length(x))

xtable(res_all)

library(ggplot2)
library(reshape2)
library(RColorBrewer)

res_all <-  melt(res_all)
colnames(res_all) <- c("Cell type", "Dataset", "Proportion") 

colorRampPalette(brewer.pal(9,"Reds"))(11)


gg1<-ggplot(res_all, aes(x=Dataset, y=Proportion)) + geom_col(aes(fill=`Cell type`)) + theme_minimal() + scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Reds"))(11))

ggsave(gg1, file="Proportions.pdf")
                                                                                                                      