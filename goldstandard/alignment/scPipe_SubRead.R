# run scPipe on 10X pbmc sample
# you can find the fastq on 10X webpage: https://support.10xgenomics.com/single-cell-gene-expression/datasets
library(scPipe)
library(Rsubread)

# annotation file path:
## exon annotation, we dont need ERCC for 10X
exon_annotation = "Homo_sapiens.GRCh38.90.gff3" 
# you can find the annotation on: ftp://ftp.ensembl.org/pub/release-90/gff3/

## cell barcode annotation
barcode_annotation_file = "barcode_annotation_scpipe_subread.csv"
## genome index, depends on which alinger you use. 
genome_index = "Homo_sapiens_GR38_90_primary_index"

# fastq file. fq_R1 should contain the transcript

fq_R1 = "Undetermined_R2.fastq.gz"
fq_R2 = "Undetermined_R1.fastq.gz"

# output folder
data_dir = "scpipe/subread"
system("mkdir scpipe/subread")
## output files:
combined_fastq = file.path(data_dir, "combined.fastq")
aligned_bam = file.path(data_dir, "out.aligned.bam")
mapped_bam = file.path(data_dir, "out.aligned.mapped.bam")

# trim barcode and move them to fastq read name.
## the read structure for 10X pbmc sample is:
## no index in read1, read2 starts with 16bp index2, then 10bp UMI. 
tenX_read_structure = list(bs1=-1, bl1=0, bs2=0, bl2=16, us=16, ul=10) # 0-index
sc_trim_barcode(combined_fastq,
               fq_R1,
               fq_R2,
               read_structure = tenX_read_structure) 

# alignment
Rsubread::align(index=genome_index,
               readfile1=combined_fastq,
               output_file=aligned_bam)

# map reads to genes
sc_exon_mapping(inbam=aligned_bam,
                outbam=mapped_bam,
                annofn=c(exon_annotation),
                bc_len = 16, UMI_len = 10) # the UMI length is 10 bp, the total index length is 16 bp.

# # detect cell barcodes
# ## the cell barcode in drop-seq like protocol is random, so we need to detect barcode before demultiplex
sc_detect_bc(infq=combined_fastq,
             outcsv=barcode_annotation_file, # output the cell barcode annotation
             bc_len=16, # 16bp index
             max_reads=5000000, # only read first 5 million reads
             min_count=100) # discard the cell barcode if it has reads lower than 100

# barcode demultiplex and UMI deduplication
sc_demultiplex(inbam=mapped_bam, outdir=data_dir, bc_anno=barcode_annotation_file)

sc_gene_counting(outdir=data_dir, 
                 bc_anno=barcode_annotation_file)

print("Done")


