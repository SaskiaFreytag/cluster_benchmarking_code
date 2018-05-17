###########################################################################################
#                                                                                         #
#                                                                                         #  
#                                 PREPROCESSING DATASET1                                  #
#                                                                                         #
#                                                                                         #
###########################################################################################


## Load data via cellranger package
library(cellrangerRkit)
gene_bc_matrix <- load_cellranger_matrix("dataset1") 

## Make scater object
library(scater)

sce <- SingleCellExperiment(
  assays = list(counts = as(exprs(gene_bc_matrix), "matrix")), colData = as(phenoData(gene_bc_matrix), 'data.frame'), rowData=as(featureData(gene_bc_matrix), 'data.frame'))


# find rRNAs ENS-IDs
library(AnnotationHub)
library(GenomicRanges)
edb <- AnnotationHub()[["AH57757"]]

anno<-as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype", "seq_name"), filter=GeneIdFilter(rownames(sce))), "data.frame")
anno<-anno[pmatch(rownames(sce),anno$gene_id), c(6:7)]
colnames(anno) <- c("id", "symbol")
anno<-as(anno, "DataFrame")
rowData(sce)<-anno
rownames(rowData(sce))<-rownames(assay(sce))

## Define mitochondrial, ribosomal RNA and protein genes

# find rRNAs ENS-IDs
rRNA<-as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype", "seq_name"), filter=GeneBiotypeFilter("rRNA")), "data.frame")


rowData(sce)$Mito <- grepl("^MT-", rowData(sce)$symbol)
rowData(sce)$ProteinRibo <- grepl("^RP", rowData(sce)$symbol)
rowData(sce)$RNARibo <- rowData(sce)$id %in% rRNA$gene_id


## Calculate quality statistics 
sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo, rbr=rowData(sce)$RNARibo))

hist(colData(sce)$pct_counts_mt)
hist(colData(sce)$pct_counts_rbr)
hist(colData(sce)$pct_counts_rbp)

## Delete bad cells 
libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
mito.drop <- isOutlier(sce$pct_counts_mt, nmads = 3, type = "higher", log = TRUE)
sce <- sce[,!(libsize.drop | feature.drop |mito.drop)]


sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo, rbr=rowData(sce)$RNARibo))


## Save list of barcodes for reanalysis
barcode<-colData(sce)$barcode
barcode<-as.matrix(barcode)
colnames(barcode)<-"Barcode"
write.csv(barcode, file="cleaned_barcodes_dataset1.csv", row.names = F, quote=F)

## Save scater object for further analysis
save(sce, file="Sce_Dataset1.RData")
