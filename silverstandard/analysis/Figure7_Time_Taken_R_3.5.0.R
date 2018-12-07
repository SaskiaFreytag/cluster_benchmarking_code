#######################################################################################
#                                                                                     #
#                                                                                     #
#                           MAKE TIME PLOT FOR PAPER                                  #
#                                                                                     #
#                                                                                     #
#######################################################################################


algorithm  <- c( "ascend", "CellRanger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran","Seurat"
               ,"TSCAN")   

output_link_1 <- "Clustering_Results/Stability/Cells"
all_data <- list.dirs(output_link_1)
all_data <- all_data[-c(1)]

res_mat <- matrix(NA, nrow=length(all_data), ncol=length(algorithm))
colnames(res_mat) <- algorithm


for (i in 1:length(all_data)){
  
  load(paste0(all_data[i], "/Time_Taken_Half_Iteration_", i, ".RData" ))
  res_mat[i, ] <- sapply(time_taken, function(x) x[4])
}

library(reshape2)

res_mat<-melt(res_mat)
colnames(res_mat) <- c("Iteration", "Method", "Running Time (sec)")

res_mat <- as.data.frame(res_mat)
res_mat$`log10 Running Time (sec)` <- log10(res_mat$`Running Time (sec)`)

library(dplyr)

res_mat_grouped <- group_by(res_mat, Method)
res_mat_grouped <- summarize(res_mat_grouped, `Mean log10 Running Time (sec)` = mean(`log10 Running Time (sec)`), `Mean Running Time (sec)` = mean(`Running Time (sec)`))


res_mat_grouped

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors_gg<- gg_color_hue(12)[-11]

res_mat_grouped <- res_mat_grouped[order(as.character(res_mat_grouped$Method)),]
res_mat_grouped$Method <- as.factor(as.character(res_mat_grouped$Method))

gg<-ggplot(res_mat_grouped, aes(y=`Mean log10 Running Time (sec)`, x=Method)) + geom_col(aes(fill = `Method`)) + theme_minimal(base_size = 11, base_family = "") + scale_fill_manual(values=colors_gg)

ggsave(file="Figure_Time_Taken.pdf", plot=gg, height=5, width=10)
