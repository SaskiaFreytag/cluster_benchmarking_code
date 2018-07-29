###########################################################################################
#                                                                                         #
#                                                                                         #
#                           Make Figure Silver Standard 2b                                #
#                                                                                         #
#                                                                                         #
###########################################################################################

library(corrplot)
library(mclust)
library(NMI)

## Load results
load("Clustering_Result_Dataset2.RData")
res_cellranger <- res

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
programs <- c( "ascend", "Cell Ranger", "CIDR", "countClust", "RaceID", "RaceID2", "RCA", "SC3", "SIMLR", "scran", "Seurat","TSCAN")         

index_cellranger <- pmatch(programs, colnames(res_cellranger))
mat1_cellranger <- comp_cluster_mat(res_cellranger[,index_cellranger])

# Calculate number of clusters in each method
clust_num_cellranger <- apply(res_cellranger[,index_cellranger], 2, function(x) length(unique(x)))

pdf("Plot_SilverStandard_2b.pdf", width=10.45, height=10.45)
plot_corr_sim(mat1_cellranger, clust_num_cellranger)
text("b", x=-1, y=12.5, cex=2.5)
dev.off()

