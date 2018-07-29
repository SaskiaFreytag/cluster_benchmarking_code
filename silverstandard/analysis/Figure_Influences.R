###########################################################################################
#                                                                                         #
#                                                                                         #
#                              FIGURE INFLUENCE OF FACTORS                                #
#                                                                                         #
#                                                                                         #
###########################################################################################


library(dplyr)
library(reshape2)
library(ggplot2)

## Load results
load("Clustering_Result_Dataset1.RData")
res_dataset1 <- res
load("Assigned_Cell_Types_Dataset1.RData")
res_dataset1$Truth <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset2.RData")
res_dataset2 <- res
load("Assigned_Cell_Types_Dataset3.RData")
res_dataset2$Truth <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset3.RData")
res_dataset3 <- res
load("Assigned_Cell_Types_Dataset3.RData")
res_dataset3$Truth <- assigned_cell_types$Assigned_CellType

programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran", "Seurat","SIMLR", "TSCAN", "Truth")

index_dataset1 <- pmatch(programs, colnames(res_dataset1))
index_dataset2 <- pmatch(programs, colnames(res_dataset2))
index_dataset3 <- pmatch(programs, colnames(res_dataset3))


## Function to find adjusted R-square for cell trait
find_adRsq <- function(res, clusters, attributes=c("total_counts", "total_features", "pct_counts_mt", "pct_counts_rbr", "pct_counts_rbp"),
                       better_lables=c("Total counts", "Total features", "Mitochondrial (%)", 
                                       "Ribosomal RNA (%)", "Ribosomal Protein (%)")){
  indices <- pmatch(attributes, colnames(res))
  
  res_mat <- matrix(NA, ncol=length(indices), nrow=length(clusters))
  colnames(res_mat) <- better_lables
  rownames(res_mat) <- colnames(res)[clusters]
  
  for (i in 1:length(indices)){
    res_mat[,i] <- sapply(clusters, function(x) summary(lm(res[,indices[i]]~as.factor(res[,x])))$adj.r.squared)
  }
  
  
  influences <- melt(res_mat)
  colnames(influences)<-c("Method", "Cell Feature", "value")

  
  return(influences)
}

influences_dataset1 <- find_adRsq(res_dataset1, index_dataset1)
influences_dataset2 <- find_adRsq(res_dataset2, index_dataset2)
influences_dataset3 <- find_adRsq(res_dataset3, index_dataset3)

influences_dataset1$Dataset <- "Dataset 1"  
influences_dataset2$Dataset <- "Dataset 2"  
influences_dataset3$Dataset <- "Dataset 3"  

influences <- rbind(influences_dataset1, influences_dataset2, influences_dataset3)
influences_mean <- group_by(influences, Method ,`Cell Feature`) %>% summarize (. , `value` = mean(`value`))

influences_mean_max <- influences_mean %>% group_by(`Cell Feature`) %>% filter(value == max(value))

influences_mean <- list(influences_mean,  influences_mean_max)

# functions for radar plot
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


coord_radar <- function (theta = "x", start = 0, direction = 1)
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# plot function for spider plot
plot_spider <- function(influences){
    
    a<- ggplot(influences[[1]], aes(x = `Cell Feature`, y = value)) +
    geom_polygon(aes(group = Method, color = Method), fill = NA, size = 1) +
    geom_line(aes(group = Method, color = Method), size = 1) +
    scale_color_manual(values=c(gg_color_hue(12), "black")) +
    annotate("text", x = influences[[2]]$`Cell Feature`, y = (influences[[2]]$value + 0.05),
    label = round(influences[[2]]$value, 2), size = 4) +
    theme(strip.text.x = element_text(size = rel(1), color="black"),
    axis.text.x = element_text(size = rel(1.5), color="black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_line(color="grey70"),
    panel.grid.minor =element_line(color="grey70"),
    legend.text=element_text(size=rel(1)),
    legend.position="bottom") +
    xlab("") + ylab("") + theme(axis.line = element_blank()) +
    guides(color = guide_legend() )  +
    coord_radar() + theme(panel.background = element_rect(fill = "transparent", colour = "transparent"))
    
    return(a)
}

gg <- plot_spider(influences_mean )

ggsave("Figure_Influences.pdf", plot=gg, height=10, width=12)
