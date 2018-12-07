###########################################################################################
#                                                                                         #
#                                                                                         #
#                              FIGURE INFLUENCE OF FACTORS                                #
#                                                                                         #
#                                                                                         #
###########################################################################################



## Load results
load("Clustering_Result_Dataset1b.RData")
res_dataset2a <- res
colnames(res_dataset2a)[77] <- "Cell Ranger"
colnames(res_dataset2a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset1b.RData")
res_dataset2a$Truth <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset3b.RData")
res_dataset3a <- res
colnames(res_dataset3a)[77] <- "Cell Ranger"
colnames(res_dataset3a)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset3b.RData")
res_dataset3a$Truth <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset4.RData")
res_dataset4 <- res
colnames(res_dataset4)[77] <- "Cell Ranger"
colnames(res_dataset4)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset4.RData")
res_dataset4$Truth <- assigned_cell_types$Assigned_CellType

load("Clustering_Result_Dataset5.RData")
res_dataset5 <- res
colnames(res_dataset5)[77] <- "Cell Ranger"
colnames(res_dataset5)[81] <- "RaceID"

load("Assigned_Cell_Types_Dataset5.RData")
res_dataset5$Truth <- assigned_cell_types$Assigned_CellType

res_all <- list(`dataset 2a`=res_dataset2a,  `dataset 3a`=res_dataset3a , `dataset 4`= res_dataset4, `dataset 5`= res_dataset5)

programs <- c("ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran", "Seurat", "TSCAN", "Truth")         

index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))


## Function to find adjusted R-square for cell trait
library(dplyr)
library(reshape2)

find_adRsq <- function(res, clusters, attributes=c("total_counts", "total_features", "pct_counts_mt", "pct_counts_rbp"),
                       better_lables=c("Total counts", "Total genes detected", "Mitochondrial (%)", "Ribosomal Protein (%)")){
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

influences_all <- lapply(1:length(index_all), function(x) find_adRsq(res_all[[x]], index_all[[x]]))



influences <- Reduce(rbind, influences_all)

library(dplyr)

influences_mean <- group_by(influences, Method ,`Cell Feature`) %>% summarize (. , `value` = mean(`value`))


influences_mean_max <- influences_mean %>% group_by(`Cell Feature`) %>% dplyr::filter(value == max(value))

influences_mean <- list(influences_mean,  influences_mean_max)

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

plot_spider <- function(influences){
  
  a<- ggplot(influences[[1]], aes(x = `Cell Feature`, y = value)) + 
    geom_polygon(aes(group = Method, color = Method), fill = NA, size = 1) +
    geom_line(aes(group = Method, color = Method), size = 1) + 
    scale_color_manual(values=c(gg_color_hue(12)[-11], "black")) +
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

ggsave("Figure_Influences_R35.pdf", plot=gg, height=10, width=12)
