##
##    Microglia_sample_comparison.R
##
##

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the microglia dataframe. 
source("R_scripts/shared/Datasets_microglia.R")

## PART 3: Setup DEseq2 analysis for all sample groups. Comparing to microglia from controls (CP6)

### DESeq2 Normalization and DE of P5 v P6
colDataD56 <- data.frame(condition=rep(NA,length(colnames(cntNiek_microglia))), row.names=colnames(cntNiek_microglia))
colDataD56$condition[grep("T.P5",rownames(colDataD56))] <- "TP5"
colDataD56$condition[grep("T..P5",rownames(colDataD56))] <- "TP5"
colDataD56$condition[grep("C.P6",rownames(colDataD56))] <- "CP6"
colDataD56$condition[grep("C..P6",rownames(colDataD56))] <- "CP6"
colDataD56$condition[grep("T.P6",rownames(colDataD56))] <- "TP6"
colDataD56$condition[grep("T..P6",rownames(colDataD56))] <- "TP6"
dds <- DESeqDataSetFromMatrix(countData = cntNiek_microglia,colData = colDataD56,design = ~condition)
dds$condition <- factor(dds$condition,levels=c("CP6",
                                               "TP6",
                                               "TP5"))

dds <- DESeq(dds)

rm(colDataD56, cntNiek_microglia)

# For comparative and descriptive analysis we use the rlog (see manual for descriptions). blind=FALSE to keep sample group info for calculating normalizations. 
rld56 <- rlog(dds, blind = FALSE)

# Extract the rld transformed values
rld56_values <- assay(rld56)

## PART 4: UNSUPERVISED CLUSTERING

# Calculate standard deviations per gene (all 9 samples) and order based on std dev (highest at top)
geneSDdds <- rep(NA, nrow(rld56_values))
for (i in 1:nrow(rld56_values)) {geneSDdds[i] <- sd(rld56_values[i,])}
rld56_values_sd_ord <- rld56_values[order(-geneSDdds),]

my_palette <- colorRampPalette(c("blue4","white", "red"))(n = 1200)

col_breaks = c(seq(6,9,length=40), # for blue
               seq(9.007519,12,length=40), # for white
               seq(12.007519,15,length=40)) # for red

pdf(file="output/pdf/microglia_unsupervised/clustering-P5-P6_170905_triplicates_750.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_sd_ord[1:750,]),
           #main="P5 and P6 Unsupervised Clustering (log2 of DESeq2 counts)", # heat map title
           col=my_palette,
           #breaks=col_breaks,
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=0.44,
           cexRow=0.44,
           keysize=0.8,
           labRow = "",
           #scale = c("row"),
           Colv=T,
           Rowv=T,
           dendrogram="column",
           colsep=1:ncol(rld56_values_sd_ord), sepcolor='black', sepwidth=0.01,
           lmat=rbind( c(0, 3), c(2,1), c(0,4) ),
           lhei=c(0.5, 16, 0.25 ) # dendrogram height, row height, key height?
)
dev.off()
