README Single Cell Clustering Benchmarking Project

goldstandard: Includes all code for the benchmarking analysis of gold standard dataset:
GSE111108
programmes: contains folder algorithms with code for all clustering tools used
silverstandard: Includes all code for the benchmarking analysis of silver standard 
datasets: GSE115189, 10xGenomics datasets 
(https://support.10xgenomics.com/single-cell-gene-expression/datasets)

The following versions of each software were used in this project:

R version 3.4.3 (2017-11-30)
Base packages: base, datasets, graphics, grDevices, methods, parallel, stats, stats4,
utils
Other packages: abind~1.4-5, acepack~1.4.1, ascend~0.4.0, ade4~1.7-11, 
AnnotationDbi~1.40.0, ape~5.1, ascend~0.1.1, assertthat~0.2.0, backports~1.1.2, 
base64enc~0.1-3, beeswarm~0.2.3, bindr~0.1.1, bindr~0.1.1, bindrcpp~0.2.2, Biobase~2.38.0, 
Biobase~2.38.0, BiocGenerics~0.24.0, BiocInstaller~1.28.0, BiocParallel~1.12.0, 
BiocParallel~1.12.0, biomaRt~2.34.2, biomaRt~2.34.2, Biostrings~2.46.0, bit~1.1-12, 
bit64~0.9-7, bit64~0.9-7, bitops~1.0-6, bitops~1.0-6, blob~1.1.1, broom~0.4.4, 
caret~6.0-79, caTools~1.17.1, cellrangerRkit~1.1.0, checkmate~1.8.5, cidr~0.1.5, 
class~7.3-14, cluster~2.0.7-1, clusterCrit~1.2.7, codetools~0.2-15, colorspace~1.3-2, 
colorspace~1.3-2, combinat~0.0-8, compiler~3.4.3, corrplot~0.84, CountClust~1.4.1, 
cowplot~0.9.2, curl~3.2, CVST~0.2-1, data.table~1.10.4-3, data.table~1.10.4-3, DBI~0.8, 
ddalpha~1.3.2, DelayedArray~0.4.1, DEoptimR~1.0-8, DEoptimR~1.0-8, devtools~1.13.5, 
diffusionMap~1.1-0, digest~0.6.15, digest~0.6.15, dimRed~0.1.0, diptest~0.75-7, 
doParallel~1.0.11, doRNG~1.6.6, doSNOW~1.0.16, dplyr~0.7.4, dplyr~0.7.4, DRR~0.0.3, 
DT~0.4, dtw~1.18-1, dynamicTreeCut~1.63-1, e1071~1.6-8, edgeR~3.20.9, fastICA~1.2-1, 
fitdistrplus~1.0-9, flexmix~2.3-14, FNN~1.1, FNN~1.1, foreach~1.4.4, foreach~1.4.4, 
foreign~0.8-69, Formula~1.2-2, fpc~2.1-11, gdata~2.18.0, GenomeInfoDb~1.14.0, 
GenomeInfoDb~1.14.0, GenomeInfoDbData~1.0.0, GenomeInfoDbData~1.0.0, 
GenomicAlignments~1.14.2, GenomicFeatures~1.30.3, GenomicRanges~1.30.3, geometry~0.3-6, 
GGally~1.3.2, ggbeeswarm~0.6.0, ggplot2~2.2.1, ggplot2~2.2.1, ggrepel~0.7.0, 
ggridges~0.5.0, git2r~0.21.0, glue~1.2.0, glue~1.2.0, GO.db~3.5.0, gower~0.1.2, 
gplots~3.0.1, graph~1.56.0, grid~3.4.3, grid~3.4.3, gridExtra~2.3, gridExtra~2.3, 
gtable~0.2.0, gtools~3.5.0, Hmisc~4.1-1, Homo.sapiens~1.3.1, htmlTable~1.11.2, 
htmltools~0.3.6, htmlwidgets~1.0, htmlwidgets~1.0, httpuv~1.3.6.2, httr~1.3.1, 
ica~1.0-1, igraph~1.2.1, infotheo~1.2.0, ipred~0.9-6, IRanges~2.12.0, IRanges~2.12.0, 
irlba~2.3.2, iterators~1.0.9, iterators~1.0.9, kernlab~0.9-25, KernSmooth~2.23-15, 
knitr~1.20, lars~1.2, lattice~0.20-35, lattice~0.20-35, latticeExtra~0.6-28, 
lava~1.6.1, lazyeval~0.2.1, limma~3.34.9, lmtest~0.9-36, locfit~1.5-9.1, lubridate~1.7.4, 
magic~1.5-8, magrittr~1.5, maptpx~1.9-2, MASS~7.3-49, MASS~7.3-49, Matrix~1.2-14, 
matrixStats~0.53.1, matrixStats~0.53.1, mclust~5.4, memoise~1.1.0, metap~0.8, 
mgcv~1.8-23, mime~0.5, minpack.lm~1.2-1, mixtools~1.1.0, mnormt~1.5-5, ModelMetrics~1.1.0, 
modeltools~0.2-21, munsell~0.4.3, mvtnorm~1.0-7, mvtnorm~1.0-7, NbClust~3.0, nlme~3.1-137, 
NMI~2.0, nnet~7.3-12, nnet~7.3-12, numDeriv~2016.8-1, org.Hs.eg.db~3.5.0, 
org.Hs.eg.db~3.5.0, OrganismDbi~1.20.0, pbapply~1.3-4, pcaPP~1.9-73, permute~0.9-4, 
pheatmap~1.0.8, picante~1.6-2, pillar~1.2.1, pillar~1.2.1, pkgconfig~2.0.1, pkgmaker~0.22, 
plyr~1.8.4, png~0.1-7, prabclus~2.2-6, pracma~2.1.4, prettyunits~1.0.2, prodlim~1.6.1, 
progress~1.1.2, proxy~0.4-22, psych~1.8.3.3, purrr~0.2.4, R.methodsS3~1.7.1, R.oo~1.21.0, 
R.utils~2.6.0, R6~2.2.2, RaceID ~ https://github.com/dgrun/RaceID, RaceID2 ~https://github.com/dgrun/StemID, ranger~0.9.0, RANN~2.5.1, RBGL~1.54.0, RBGL~1.54.0, RCA~1.0, 
RColorBrewer~1.1-2, Rcpp~0.12.16, RcppAnnoy~0.0.10, RcppEigen~0.3.3.4.0, 
RcppParallel~4.4.0, RcppRoll~0.2.2, RCurl~1.95-4.10, RCurl~1.95-4.10, recipes~0.1.2, 
registry~0.5, reshape~0.8.7, reshape2~1.4.3, rhdf5~2.22.0, rhdf5~2.22.0, Rhtslib~1.10.0, 
rjson~0.2.15, rlang~0.2.0, Rmisc~1.5, RMySQL~0.10.14, rngtools~1.2.4, robustbase~0.92-8, 
robustbase~0.92-8, ROCR~1.0-7, rpart~4.1-13, rrcov~1.4-3, Rsamtools~1.30.0, 
RSpectra~0.12-0, RSQLite~2.1.0, Rsubread~1.28.1, rstudioapi~0.7, rtracklayer~1.38.3, Rtsne~0.13, 
S4Vectors~0.16.0, SC3~1.7.7, scales~0.5.0, scales~0.5.0, scater~1.6.3, 
scatterplot3d~0.3-41, scPipe~1.0.6, scran~1.6.9, SDMTools~1.1-221, segmented~0.5-3.0, 
Seurat~2.3.0, sfsmisc~1.1-2, shiny~1.0.5, shinydashboard~0.7.0, SIMLR~1.4.1, 
SingleCellExperiment~1.0.0, SingleCellExperiment~1.0.0, slam~0.1-42, sn~1.5-1, snow~0.4-2, 
splines~3.4.3, statmod~1.4.30, stringi~1.1.7, stringi~1.1.7, stringr~1.3.0, 
SummarizedExperiment~1.8.1, survival~2.41-3, tclust~1.3-1, tibble~1.4.2, tidyr~0.8.0, 
tidyselect~0.2.4, timeDate~3043.102, tools~3.4.3, tools~3.4.3, trimcluster~0.1-2, 
TSCAN~1.16.0, tsne~0.1-3, TxDb.Hsapiens.UCSC.hg19.knownGene~3.2.2, tximport~1.6.0, 
vegan~2.5-1, VGAM~1.0-5, vipor~0.4.5, viridis~0.5.1, viridisLite~0.3.0, withr~2.1.2, 
WriteXLS~4.0.0, XML~3.98-1.10, xtable~1.8-2, XVector~0.18.0, XVector~0.18.0, yaml~2.1.18, 
zlibbioc~1.24.0, zoo~1.8-1
Other software
demuxlet version 0.0.1
CellRanger version 2.0.0
STAR version 020201
Subread version 1.5.2)
