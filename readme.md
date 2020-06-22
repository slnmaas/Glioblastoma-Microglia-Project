# Glioblastoma Microglia Project

This repository contains all R files used to generate the figures for the manuscript “Glioblastoma Hijacks Microglial Transcriptional Networks to ” by Sybren L.N. Maas, Erik R. Abels, Lieke L. Van De Haar, Xuan Zhang, Liza Morsett, Srinjoy Sil, Joana Guedes, Pritha Sen, Shilpa Prabhakar, Suzanne E. Hickman, Charles P. Lai David T. Ting, Xandra O. Breakefield, Marike L.D. Broekman, and Joseph El Khoury.

This manuscript was published on June 18th 2020 in <a href="https://jneuroinflammation.biomedcentral.com/articles/10.1186/s12974-020-01797-2" target="_blank">Journal of Neuroinflammation</a>. 

Using these scripts for other scientific projects is encouraged. However, if these scripts help the generation of a new manuscript please cite the publication listed above.

## Getting Started

To run these files R is required. For information on the R programming language see: https://www.r-project.org/
For this project we used RStudio to run the R files. See https://www.rstudio.com/ for more information.

Data can be downloaded from the GEO repositories listed in the manuscript. 

### Prerequisites

The following packages are required:

```
install.packages("DESeq2")
install.packages("gplots")
install.packages("limma")
install.packages("VennDiagram")
install.packages("UpSetR")
```

### sessionInfo()

For compatibility information, the scripts were executed with the following versions:

```
> sessionInfo()
R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.12.6 (unknown)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] UpSetR_1.3.3               VennDiagram_1.6.16         futile.logger_1.4.1        limma_3.26.3              
 [5] gplots_2.17.0              DESeq2_1.10.0              RcppArmadillo_0.6.300.2.0  Rcpp_0.12.12              
 [9] SummarizedExperiment_1.0.1 Biobase_2.30.0             GenomicRanges_1.22.2       GenomeInfoDb_1.6.1        
[13] IRanges_2.4.6              S4Vectors_0.8.5            BiocGenerics_0.16.1       

loaded via a namespace (and not attached):

 [1] RColorBrewer_1.1-2   plyr_1.8.3           XVector_0.10.0       bitops_1.0-6         futile.options_1.0.0
 [6] tools_3.2.3          zlibbioc_1.16.0      rpart_4.1-10         RSQLite_1.0.0        annotate_1.48.0     
[11] gtable_0.1.2         lattice_0.20-33      DBI_0.3.1            proto_0.3-10         gridExtra_2.0.0     
[16] genefilter_1.52.0    cluster_2.0.3        caTools_1.17.1       gtools_3.5.0         locfit_1.5-9.1      
[21] nnet_7.3-11          AnnotationDbi_1.32.2 XML_3.98-1.3         survival_2.38-3      BiocParallel_1.4.0  
[26] foreign_0.8-66       gdata_2.17.0         latticeExtra_0.6-26  Formula_1.2-1        geneplotter_1.48.0  
[31] ggplot2_2.0.0        lambda.r_1.1.7       scales_0.4.1         Hmisc_3.17-0         splines_3.2.3       
[36] colorspace_1.2-6     xtable_1.8-0         KernSmooth_2.23-15   acepack_1.3-3.3      munsell_0.4.2  
```
