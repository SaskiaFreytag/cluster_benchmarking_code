###########################################################################################
#                                                                                         #
#                                                                                         #
#                          PREPROCESSING DATA OF SCPIPE OUTPUT                            #
#                                                                                         #
#                                                                                         #
###########################################################################################


## Change functions in scpipe to avoid biomaRt
library(Homo.sapiens)
library(biomaRt)
library(scPipe)

get_genes_by_go_alt <- function (go = NULL) 
{
  if (is.null(go)) {
    stop("must provide GO term. (i.e go=c('GO:0005739'))")
  }

  G_list <- select(
    Homo.sapiens,
    key = go,
    keytype = "GO",
    columns = c("ENSEMBL", "GENENAME")
  )
  return(G_list[, 4])
}


## Load data via csv 
library(scPipe)
library(SingleCellExperiment)
sce = create_sce_by_dir("scpipe/STAR")
dim(sce)
gene_number <- colSums(assay(sce, "counts") > 0)
QC_metrics(sce)$number_of_genes = gene_number
QC_metrics(sce)$total_count_per_cell = colSums(assay(sce, 
                                                     "counts"))
mt_genes <- get_genes_by_go_alt(go = c("GO:0005739"))
mt_count <- colSums(assay(sce, "counts")[rownames(sce) %in% 
                                           mt_genes, ])
QC_metrics(sce)$non_mt_percent <- (colSums(assay(sce, 
                                                 "counts")) - mt_count)/(colSums(assay(sce, 
                                                                                       "counts")) + 1e-05)
ribo_genes <- get_genes_by_go_alt (go = c("GO:0005840"))
ribo_count <- colSums(assay(sce, "counts")[rownames(sce) %in% 
                                             ribo_genes, ])
QC_metrics(sce)$non_ribo_percent = (colSums(assay(sce, 
                                                  "counts")) - ribo_count)/(colSums(assay(sce, 
                                                                                          "counts")) + 0.01)

plot_demultiplex(sce)
plot_UMI_dup(sce)


#sce = calculate_QC_metrics(sce)
sce = detect_outlier(sce)

plot_mapping(sce, percentage = TRUE, dataname = "sc_sample_data")
plot_QC_pairs(sce)

sce = remove_outliers(sce)

library(AnnotationHub)
library(GenomicRanges)
#query(AnnotationHub(), "EnsDb.Hsapiens.")
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

# find nuclear mitochondrial ENS-IDs 
nuclear_mito<-read.table("~/10xGenomics/Clustering/Reference/Mito_Genes_Nuclear.txt", sep="\t", fill=T, quote="")

rowData(sce)$Mito <- grepl("MT", rowData(sce)$symbol)
rowData(sce)$ProteinRibo <- grepl("^RP", rowData(sce)$symbol)
rowData(sce)$RNARibo <- rowData(sce)$id %in% rRNA$gene_id


## Calculate quality statistics 
library(scater)
sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo, rbr=rowData(sce)$RNARibo))

hist(colData(sce)$pct_counts_mt)
hist(colData(sce)$pct_counts_rbr)
hist(colData(sce)$pct_counts_rbp)

## Delete bad cells 
libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
mito.drop <- isOutlier(sce$pct_counts_mt, nmads = 3, type = "higher", log = TRUE) 
sce <- sce[,!(libsize.drop | feature.drop| mito.drop)]


## Save barcodes
#Load barcodes
all_barcodes <- read.csv("barcode_annotation_scpipe_star.csv", header=F)
#Match barcodes
barcode <- all_barcodes$V2[pmatch(colnames(sce), all_barcodes$V1)]
colData(sce)$barcode <- barcode
barcode<-as.matrix(barcode)
colnames(barcode)<-"Barcode"
write.csv(barcode, file="cleaned_barcodes_scpipe_star.csv", row.names = F, quote=F)

sce <- calculateQCMetrics(sce, feature_controls = list(mt=rowData(sce)$Mito , rbp=rowData(sce)$ProteinRibo, rbr=rowData(sce)$RNARibo))


## Save sce object for further analysis
save(sce, file="Sce_scPipe_STAR.RData")



