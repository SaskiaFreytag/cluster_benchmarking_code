###########################################################################################
#                                                                                         #
#                                                                                         #  
#                     COMPARE CELLRANGER AND SCPIPE OUTCOMES                              #
#                                                                                         #
#                                                                                         #
###########################################################################################

## Load libraries
library(SingleCellExperiment)
library(UpSetR)
library(scater)
library(ggplot2)

## Load scPipe object and resave
load("Sce_ScPipe_SubRead.RData")
sce_scpipe <- sce
remove(sce)

## Load cellranger object and resave
load("Sce_CellRanger.RData")
sce_cellranger <- sce
colData(sce_cellranger)$barcode<-gsub("-1", "", colData(sce_cellranger)$barcode)
remove(sce)

## Load scPipe_STAR object and resave
load("Sce_ScPipe_STAR.RData")
sce_STARpipe <- sce
remove(sce)


## Plot Venn diagram of the barcodes
all_barcodes <- union(as.character(colData(sce_cellranger)$barcode), union(as.character(colData(sce_scpipe)$barcode), as.character(colData(sce_STARpipe)$barcode)))
barcode_matrix <- cbind(as.numeric(all_barcodes %in% as.character(colData(sce_cellranger)$barcode)), 
                        as.numeric(all_barcodes %in% as.character(colData(sce_scpipe)$barcode)),
                        as.numeric(all_barcodes %in% as.character(colData(sce_STARpipe)$barcode)))
rownames(barcode_matrix) <- all_barcodes
colnames(barcode_matrix) <- c("Cell Ranger", "ScPipe Subread", "ScPipe STAR")

pdf("Barcodes.pdf", height=5, width=6)
upset(as.data.frame(barcode_matrix))
dev.off()

## Plot Venn diagram of the genes
all_genes <- union(rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0], union(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0]))
gene_matrix <- cbind(as.numeric(all_genes %in% rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0]), as.numeric(all_genes %in% rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0]), as.numeric(all_genes %in% rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0]))
rownames(gene_matrix) <- all_genes
colnames(gene_matrix) <- c("Cell Ranger", "ScPipe Subread", "ScPipe STAR")

pdf("Genes.pdf", height=5, width=6)
upset(as.data.frame(gene_matrix))
dev.off()

## Find out what these genes are predominatley
genes_unique_cellranger<- setdiff(rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0], union(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0]))
genes_unique_scpipe<- setdiff(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], union(rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0], rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0]))
genes_unique_STARpipe<- setdiff(rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0], union(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0]))
genes_unique_pipe<- setdiff(union(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0]), rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0])
genes_unique_cell_STAR<- setdiff(union(rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0], rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0]), rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0])
genes_unique_cell_scpipe<- setdiff(union(rownames(sce_scpipe)[rowData(sce_scpipe)$total_counts>0], rownames(sce_cellranger)[rowData(sce_cellranger)$total_counts>0]), rownames(sce_STARpipe)[rowData(sce_STARpipe)$total_counts>0])

library(AnnotationHub)
library(GenomicRanges)
edb <- AnnotationHub()[["AH57757"]]

anno_unique_cellranger <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_cellranger)), "data.frame")
anno_unique_scpipe <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_scpipe)), "data.frame")
anno_unique_STARpipe <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_STARpipe)), "data.frame")
anno_unique_pipe <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_pipe)), "data.frame")
anno_unique_cell_STAR <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_cell_STAR)), "data.frame")
anno_unique_cell_scpipe <- as(genes(edb, columns = c("gene_id","symbol","gene_name", "gene_biotype"), filter=GeneIdFilter(genes_unique_cell_scpipe)), "data.frame")

anno_unique_genes <- rbind(cbind(anno_unique_cellranger$gene_biotype, "Cell Ranger"), cbind(anno_unique_scpipe$gene_biotype, "ScPipe Subread"), cbind(anno_unique_STARpipe$gene_biotype, "ScPipe STAR"), cbind(anno_unique_pipe$gene_biotype, "ScPipe"), cbind(anno_unique_cell_STAR$gene_biotype, "ScPipe STAR + Cell Ranger"), cbind(anno_unique_cell_scpipe$gene_biotype, "ScPipe Subread + Cell Ranger"))
anno_unique_genes <- as.data.frame(anno_unique_genes)
colnames(anno_unique_genes) <- c("Biotype", "Program")
anno_unique_genes$Biotype<-as.factor(anno_unique_genes$Biotype)
anno_unique_genes$Program<-as.factor(anno_unique_genes$Program)

ggsave("Classes_Unique_Genes.pdf")
ggplot(anno_unique_genes, aes(Biotype)) + geom_bar(aes(fill=Program), position = position_stack(reverse = TRUE)) + coord_flip() + ggtitle("Unique Genes")
dev.off()

## Print number of reads
intersect_barcodes <- intersect(as.character(colData(sce_cellranger)$barcode), intersect(as.character(colData(sce_scpipe)$barcode), as.character(colData(sce_STARpipe)$barcode)))
all_counts <- cbind(colData(sce_cellranger)$total_counts[pmatch(intersect_barcodes, colData(sce_cellranger)$barcode)], colData(sce_scpipe)$total_counts[pmatch(intersect_barcodes, colData(sce_scpipe)$barcode)],
                    colData(sce_STARpipe)$total_counts[pmatch(intersect_barcodes, colData(sce_STARpipe)$barcode)])
rownames(all_counts) <- intersect_barcodes
colnames(all_counts) <- c("Cell Ranger", "ScPipe Subread", "ScPipe STAR")
all_counts <- all_counts[order(all_counts[,1], decreasing=FALSE),]
all_counts <- as.data.frame(all_counts)


g1 <- ggdraw(ggplot(all_counts, aes(x=`Cell Ranger`, y=`ScPipe Subread`)) + geom_point()  + theme_minimal()) +draw_figure_label("a")
g2 <- ggdraw(ggplot(all_counts, aes(x=`Cell Ranger`, y=`ScPipe STAR`)) + geom_point()  + theme_minimal()) + draw_figure_label("b")
g3 <- ggdraw(ggplot(all_counts, aes(x=`ScPipe STAR`, y=`ScPipe Subread`)) + geom_point() + theme_minimal() ) + draw_figure_label("c")

ggsave(g1, file="Comparison_Total_Counts_1.pdf")
ggsave(g2, file="Comparison_Total_Counts_2.pdf")
ggsave(g3, file="Comparison_Total_Counts_3.pdf")



## Correlations for every sample
sub_cellranger <- counts(sce_cellranger[(rowData(sce_cellranger)$total_counts>0), pmatch(intersect_barcodes, colData(sce_cellranger)$barcode)])
sub_scpipe <- counts(sce_scpipe[(rowData(sce_scpipe)$total_counts>0), pmatch(intersect_barcodes, colData(sce_scpipe)$barcode)])
sub_STARpipe <- counts(sce_STARpipe[(rowData(sce_STARpipe)$total_counts>0), pmatch(intersect_barcodes, colData(sce_STARpipe)$barcode)])

genes <- intersect(rownames(sub_cellranger), intersect(rownames(sub_scpipe), rownames(sub_STARpipe)))
sub_cellranger <- sub_cellranger[pmatch(genes, rownames(sub_cellranger)),]
sub_scpipe <- sub_scpipe[pmatch(genes, rownames(sub_scpipe)),]
sub_STARpipe <- sub_STARpipe[pmatch(genes, rownames(sub_STARpipe)),]

corrs <- cbind(sapply(1:dim(sub_scpipe)[1], function(x) cor(sub_scpipe[x,], sub_cellranger[x,], method="spearman")),
               sapply(1:dim(sub_STARpipe)[1], function(x) cor(sub_STARpipe[x,], sub_cellranger[x,], method="spearman")),
               sapply(1:dim(sub_STARpipe)[1], function(x) cor(sub_STARpipe[x,], sub_scpipe[x,], method="spearman")))
colnames(corrs) <- c("Cell Ranger ScPipe Subread", "Cell Ranger ScPipe STAR", "ScPipe Subread ScPipe STAR")
rownames(corrs) <- genes
summary(corrs)

## Split correlations by gene groups
corrs <- as.data.frame(corrs)
corrs$Gene <- rownames(corrs)
library(reshape2)
corrs_all <- melt(corrs)
colnames(corrs_all) <- c("Gene", "Comparison", "Correlation")


anno_corrs <- as(genes(edb, columns = c("gene_biotype"), filter=GeneIdFilter(corrs_all$Gene)), "data.frame")
corrs_all$Category <- anno_corrs$gene_biotype[match(corrs_all$Gene, anno_corrs$gene_id)]

corrs_all_filt <- corrs_all[which(corrs_all$Category == "unprocesses_pseudogene"| corrs_all$Category == "protein_coding" | corrs_all$Category == "processed_pseudogene" | corrs_all$Category == "lincRNA" |
  corrs_all$Category == "antisense_RNA"),]



gg1 <- ggplot(corrs_all_filt, aes(x=Category, y=Correlation)) + geom_boxplot(aes(fill=Comparison)) + theme_minimal()
ggsave(gg1, file="Correlations.pdf", width=12, height=5)





