##
##    DE_analysis_microglia_V2.R
##

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the dataframe. We use the minus get microglia datasets file
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

## PART 4: Do differential expression analysis for groups TP6 vs CP6, TP5 vs CP6, TP5 vs TP6

# Do differential expression analysis

# Group 1: compare Tumor microglia GFP- (TP6), to Control microglia GFP- (CP6)
dds_resT6C6 <- results(dds, c("condition","TP6","CP6"))
DESeq2::plotMA(dds_resT6C6,  alpha = 0.05, main="MA plot TP6 vs CP6", ylim=c(-5,5))

# Obtain the three groups that have to be plotted: padj = NA, padj >= 0.05 and padj < 0.05
dds_resT6C6_NA <- dds_resT6C6[is.na(dds_resT6C6$padj),]
dds_resT6C6_nonNA <- dds_resT6C6[!is.na(dds_resT6C6$padj),]
dds_resT6C6_sig <- dds_resT6C6_nonNA[dds_resT6C6_nonNA$padj<0.05,]
dds_resT6C6_nonSig <- dds_resT6C6_nonNA[dds_resT6C6_nonNA$padj>=0.05,]

# Write tables
write.table(dds_resT6C6_NA,file="output/tables/MA_plot_data/microglia/resT6C6_MA_plot_data_NA_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT6C6_sig,file="output/tables/MA_plot_data/microglia/resT6C6_MA_plot_data_sig_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT6C6_nonSig,file="output/tables/MA_plot_data/microglia/resT6C6_MA_plot_data_nonSig_170215.txt",col.names=T,quote=F,sep="\t")


# Group 2: compare Tumor microglia GFP+ (TP5), to Control microglia GFP- (CP6)
dds_resT5C6 <- results(dds, c("condition","TP5","CP6"))
DESeq2::plotMA(dds_resT5C6,  alpha = 0.05, main="MA plot TP5 vs CP6", ylim=c(-5,5))

# Obtain the three groups that have to be plotted: padj = NA, padj >= 0.05 and padj < 0.05
dds_resT5C6_NA <- dds_resT5C6[is.na(dds_resT5C6$padj),]
dds_resT5C6_nonNA <- dds_resT5C6[!is.na(dds_resT5C6$padj),]
dds_resT5C6_sig <- dds_resT5C6_nonNA[dds_resT5C6_nonNA$padj<0.05,]
dds_resT5C6_nonSig <- dds_resT5C6_nonNA[dds_resT5C6_nonNA$padj>=0.05,]

# Write tables
write.table(dds_resT5C6_NA,file="output/tables/MA_plot_data/microglia/resT5C6_MA_plot_data_NA_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT5C6_sig,file="output/tables/MA_plot_data/microglia/resT5C6_MA_plot_data_sig_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT5C6_nonSig,file="output/tables/MA_plot_data/microglia/resT5C6_MA_plot_data_nonSig_170215.txt",col.names=T,quote=F,sep="\t")


# Group 3: compare Tumor microglia GFP+ (TP5), to Tumor microglia GFP- (TP6)
dds_resT5T6 <- results(dds, c("condition","TP5","TP6"))
DESeq2::plotMA(dds_resT5T6,  alpha = 0.05, main="MA plot TP5 vs TP6", ylim=c(-5,5))

# Obtain the three groups that have to be plotted: padj = NA, padj >= 0.05 and padj < 0.05
dds_resT5T6_NA <- dds_resT5T6[is.na(dds_resT5T6$padj),]
dds_resT5T6_nonNA <- dds_resT5T6[!is.na(dds_resT5T6$padj),]
dds_resT5T6_sig <- dds_resT5T6_nonNA[dds_resT5T6_nonNA$padj<0.05,]
dds_resT5T6_nonSig <- dds_resT5T6_nonNA[dds_resT5T6_nonNA$padj>=0.05,]

# Write tables
write.table(dds_resT5T6_NA,file="output/tables/MA_plot_data/microglia/resT5T6_MA_plot_data_NA_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT5T6_sig,file="output/tables/MA_plot_data/microglia/resT5T6_MA_plot_data_sig_170215.txt",col.names=T,quote=F,sep="\t")
write.table(dds_resT5T6_nonSig,file="output/tables/MA_plot_data/microglia/resT5T6_MA_plot_data_nonSig_170215.txt",col.names=T,quote=F,sep="\t")
