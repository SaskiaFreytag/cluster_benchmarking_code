###########################################################################################
#                                                                                         #
#                                                                                         #
#                     ASSIGNING CLUSTER NAMES VIA 10X METHOD DATASET3                     #
#                                                                                         #
#                                                                                         #
###########################################################################################


## Load required packages and code 
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)
library(plyr)
library(data.table)
library(pheatmap)
library(devtools)
source_url("https://raw.githubusercontent.com/10XGenomics/single-cell-3prime-paper/master/pbmc68k_analysis/util.R")

## Load dataset and reference datasets
load("Sce_Dataset3.RData")

# Download data from https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
pure_11 <-  readRDS("../Reference/all_pure_select_11types.rds")
pbmc_68k <- readRDS('../Reference/all_pure_pbmc_data.rds')

purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$all_data$`15852`$hg19$genes)
m<-t(counts(sce))
l<-.normalize_by_umi(m)   
m_n<-l$m
df<-.get_variable_gene(m_n) 

# do for 1000 most variable genes
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off

set.seed(0)
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)

m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-rowData(sce)$symbol[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-rowData(sce)$id[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11) 
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])

colData(sce)$Assigned_CellType<-cls_id

assigned_cell_types<-colData(sce)
save(assigned_cell_types, file="Assigned_Cell_Types_Dataset3.RData")

