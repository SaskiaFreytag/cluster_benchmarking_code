###########################################################################################
#                                                                                         #
#                                                                                         #
#                           Make Figure Silver Standard 4b                                #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(corrplot)
library(mclust)
library(NMI)


## Load results
load("Clustering_Result_Dataset2b.RData")
res_dataset2a <- res
colnames(res_dataset2a)[77] <- "Cell Ranger"
colnames(res_dataset2a)[81] <- "RaceID"

load("Clustering_Result_Dataset3b.RData")
res_dataset3a <- res
colnames(res_dataset3a)[77] <- "Cell Ranger"
colnames(res_dataset3a)[81] <- "RaceID"

load("Clustering_Result_Dataset4.RData")
res_dataset4 <- res
colnames(res_dataset4)[77] <- "Cell Ranger"
colnames(res_dataset4)[81] <- "RaceID"

load("Clustering_Result_Dataset5.RData")
res_dataset5 <- res
colnames(res_dataset5)[77] <- "Cell Ranger"
colnames(res_dataset5)[81] <- "RaceID"


res_all <- list(`dataset 2a`=res_dataset2a,  `dataset 3a`=res_dataset3a , `dataset 4`= res_dataset4, `dataset 5`= res_dataset5)

# Find index of clustering results
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "scran",  "Seurat"
               ,"TSCAN")         
index_all <- lapply(res_all, function(x) pmatch(programs, colnames(x)))


# Function to compare clustering results
comp_cluster_mat<-function(res_clusters1){
  
  a<-dim(res_clusters1)[2]-1
  b<-dim(res_clusters1)[2]
  nc<-dim(res_clusters1)[1]
  comp<-sapply(1:a, function(x) sapply((x+1):b, function(y) adjustedRandIndex(res_clusters1[, x], res_clusters1[, y])))
  comp2<-sapply(1:a, function(x) sapply((x+1):b, function(y) NMI(cbind(1:nc, res_clusters1[, x]), cbind(1:nc, res_clusters1[, y]))$value))
  
  comp_mat <- matrix(0, nrow=b, ncol=b)
  comp_mat2<-comp_mat
  comp_mat[lower.tri(comp_mat)] <- unlist(comp)
  
  comp_mat2[lower.tri(comp_mat2)]<-unlist(comp2)
  comp_mat2<-t(comp_mat2)
  
  comp_mat<-comp_mat+comp_mat2
  rownames(comp_mat)<-colnames(res_clusters1)
  colnames(comp_mat)<-colnames(res_clusters1)
  return(comp_mat)
}

## Function for plotting
plot_corr_sim <- function(mat1, clust_num){
  mat_sim <- mat1
  mat_sim[upper.tri(mat_sim)] <- mat1[lower.tri(mat_sim)]
  order_methods <- corrplot::corrplot(mat_sim, is.corr=TRUE, order="FPC")
  order_methods <- pmatch(rownames(order_methods), rownames(mat1))
  mat1 <- mat1[order_methods, order_methods]
  clust_num <- clust_num[order_methods]
  
  cols1<-colorRampPalette(c("white","yellow","purple"))(200)
  return({
    corrplot::corrplot(mat1, cl.lim=c(0,1), col=cols1, diag=FALSE, method="color", is.corr=FALSE, tl.col="black", 
                       addgrid.col="white", mar = c(0, 1, 1, 0))
    
    for(i in 1:length(clust_num)){
      text(i, (dim(mat1)[1]-i+1), labels=clust_num[i])
    }
    mtext(side=2,line=2, text="ARI_comp", col="darkgrey", cex=1.2)
    mtext(side=3,line=2, text="NMI", col="darkgrey", cex=1.2)
  })
}


## Calculate similarity, order them according to ARI and NMI and plot 

mat_all <- lapply(1:length(res_all), function(x) comp_cluster_mat(res_all[[x]][,index_all[[x]]]))
mat_all1 <- Reduce("+", mat_all)/4 

clust_num_1<- lapply(1:length(res_all), function(x1) apply(res_all[[x1]][,index_all[[x1]]], 2, function(x) length(unique(x))))

clust_num_1 <- Reduce("+", clust_num_1)/4

pdf("Plot_SilverStandard_Compare_3.5.0.pdf", width=10.45, height=10.45)
plot_corr_sim(mat_all1, clust_num_1)
text("b", x=-1, y=12.5, cex=2.5)
dev.off()

