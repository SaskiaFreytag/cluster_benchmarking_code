###########################################################################################
#                                                                                         #
#                                                                                         #  
#                                PRE-PROCESSING DATASET2b                                 #
#                                                                                         #
#                                                                                         #
###########################################################################################


## Load data via TENxPBMCData package
library(TENxPBMCData)
sce <- TENxPBMCData(dataset = "pbmc3k")

## Load scater
library(scater)

## Define mitochondrial, ribosomal protein genes 

rowData(sce)$Mito <- grepl("^MT", rowData(sce)$Symbol)
rowData(sce)$ProteinRibo <- grepl("^RP", rowData(sce)$Symbol)

## Calculate quality statistics 
sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo))

hist(colData(sce)$pct_counts_mt)
hist(colData(sce)$pct_counts_rbp)

## Delete bad cells 
libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
mito.drop <- isOutlier(sce$pct_counts_mt, nmads = 3, type = "higher", log = TRUE)
sce <- sce[,!(libsize.drop | feature.drop |mito.drop)]

sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo))

rowData(sce)$symbol <- rowData(sce)$Symbol_TENx
rowData(sce)$id <- rowData(sce)$ENSEMBL_ID
colData(sce)$barcode <- colData(sce)$Barcode
counts(sce) <- as.matrix(counts(sce))
colnames(sce) <- colData(sce)$barcode
rownames(counts(sce)) <- rowData(sce)$id
sce@rowRanges <- as(sce@rowRanges, "GRangesList")


## Save list of barcodes for reanalysis
barcode<-colData(sce)$Barcode
barcode<-as.matrix(barcode)
colnames(barcode)<-"Barcode"
write.csv(barcode, file="cleaned_barcodes_dataset2b.csv", row.names = F, quote=F)

## Save scater object for further analysis
save(sce, file="Sce_Dataset2b.RData")
