##
##    XY_DE_microglia_V2.R
##
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

## PART 4: Do differential expression analysis for all possible comparisons

# Group 1: compare Tumor microglia GFP- (TP6), to Control microglia GFP- (CP6)
dds_resT6C6 <- results(dds, c("condition","TP6","CP6"))
#DESeq2::plotMA(dds_resT6C6,  alpha = 0.05, main="MA plot TP6 vs CP6", ylim=c(-5,5))

# Group 2: compare Tumor microglia GFP+ (TP5), to Control microglia GFP- (CP6)
dds_resT5C6 <- results(dds, c("condition","TP5","CP6"))

# Group 3: compare Tumor microglia GFP+ (TP5), to Tumor microglia GFP- (TP6)
dds_resT5T6 <- results(dds, c("condition","TP5","TP6"))

# Group 4: compare Tumor microglia GFP- (TP6), to Tumor microglia GFP+ (TP5)
dds_resT6T5 <- results(dds, c("condition","TP6","TP5"))
#DESeq2::plotMA(dds_resT6T5,  alpha = 0.05, main="MA plot TP6 vs TP5", ylim=c(-5,5))

# First order the datasets on row names so that output is always alphabetically
dds_resT6C6Ord = dds_resT6C6[order(row.names(dds_resT6C6)),]
dds_resT5C6Ord = dds_resT5C6[order(row.names(dds_resT5C6)),]
dds_resT5T6Ord = dds_resT5T6[order(row.names(dds_resT5T6)),]
dds_resT6T5Ord = dds_resT6T5[order(row.names(dds_resT6T5)),]

rm(dds_resT5T6, dds_resT5C6, dds_resT6C6, dds_resT6T5, colDataD56, cntNiek_microglia, dds)



# COMPARISON 1: First work on XY plot for TP5 vs CP6 and TP6 vs CP6.
All_TP5_CP6_names = row.names(dds_resT5C6Ord)
dds_resT5C6OrdNN <- dds_resT5C6Ord[!is.na(dds_resT5C6Ord$padj),]
Sig_TP5_CP6 = row.names(dds_resT5C6OrdNN[dds_resT5C6OrdNN$padj<0.05,])

All_TP6_CP6_names = row.names(dds_resT6C6Ord)
dds_resT6C6OrdNN <- dds_resT6C6Ord[!is.na(dds_resT6C6Ord$padj),]
Sig_TP6_CP6 = row.names(dds_resT6C6OrdNN[dds_resT6C6OrdNN$padj<0.05,])

# Shared genes is equal to the union of the sets
Union_Sig_TP5CP6_TP6CP6 = intersect(Sig_TP5_CP6, Sig_TP6_CP6)

# Extract genes that are only significant in the specific analyses
Only_Sig_TP5CP6 = setdiff(Sig_TP5_CP6, Union_Sig_TP5CP6_TP6CP6)
Only_Sig_TP6CP6 = setdiff(Sig_TP6_CP6, Union_Sig_TP5CP6_TP6CP6)

All_sig_genes_TP5CP6_TP6CP6 = unique(c(Union_Sig_TP5CP6_TP6CP6, Only_Sig_TP5CP6, Only_Sig_TP6CP6))
All_non_sig_genes_TP5CP6 = setdiff(All_TP5_CP6_names, All_sig_genes_TP5CP6_TP6CP6)

# Now write all needed .txt files
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP5CP6_All_non_sig_genes_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Only_Sig_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP5CP6_Only_Sig_TP5CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Only_Sig_TP6CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP5CP6_Only_Sig_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Union_Sig_TP5CP6_TP6CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP5CP6_Union_Sig_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")

write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP6CP6_All_non_sig_genes_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Only_Sig_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP6CP6_Only_Sig_TP5CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Only_Sig_TP6CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP6CP6_Only_Sig_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Union_Sig_TP5CP6_TP6CP6), 
            file="output/tables/XY_plots/TP5CP6_TP6CP6/TP6CP6_Union_Sig_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")





# COMPARISON 2: First work on XY plot for TP5 vs CP6 and TP5 vs TP6.
All_TP5_CP6_names = row.names(dds_resT5C6Ord)
dds_resT5C6OrdNN <- dds_resT5C6Ord[!is.na(dds_resT5C6Ord$padj),]
Sig_TP5_CP6 = row.names(dds_resT5C6OrdNN[dds_resT5C6OrdNN$padj<0.05,])

All_TP5_TP6_names = row.names(dds_resT5T6Ord)
dds_resT5T6OrdNN <- dds_resT5T6Ord[!is.na(dds_resT5T6Ord$padj),]
Sig_TP5_TP6 = row.names(dds_resT5T6OrdNN[dds_resT5T6OrdNN$padj<0.05,])

Union_Sig_TP5CP6_TP5TP6 = intersect(Sig_TP5_CP6, Sig_TP5_TP6)

Only_Sig_TP5CP6 = setdiff(Sig_TP5_CP6, Union_Sig_TP5CP6_TP5TP6)
Only_Sig_TP5TP6 = setdiff(Sig_TP5_TP6, Union_Sig_TP5CP6_TP5TP6)

All_sig_genes_TP5CP6_TP5TP6 = unique(c(Union_Sig_TP5CP6_TP5TP6, Only_Sig_TP5CP6, Only_Sig_TP5TP6))
All_non_sig_genes_TP5CP6 = setdiff(All_TP5_CP6_names, All_sig_genes_TP5CP6_TP5TP6)

# Now write all needed .txt files
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5CP6_All_non_sig_genes_TP5CP6_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Only_Sig_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5CP6_Only_Sig_TP5CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Only_Sig_TP5TP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5CP6_Only_Sig_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% Union_Sig_TP5CP6_TP5TP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5CP6_Union_Sig_TP5CP6_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")

write.table(subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5TP6_All_non_sig_genes_TP5CP6_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% Only_Sig_TP5CP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5TP6_Only_Sig_TP5CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% Only_Sig_TP5TP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5TP6_Only_Sig_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% Union_Sig_TP5CP6_TP5TP6), 
            file="output/tables/XY_plots/TP5CP6_TP5TP6/TP5TP6_Union_Sig_TP5CP6_TP5TP6_170805.txt",col.names=T,quote=F,sep="\t")





# COMPARISON 3: First work on XY plot for TP6 vs CP6 and TP6 vs TP5.
All_TP6_CP6_names = row.names(dds_resT6C6Ord)
dds_resT6C6OrdNN <- dds_resT6C6Ord[!is.na(dds_resT6C6Ord$padj),]
Sig_TP6_CP6 = row.names(dds_resT6C6OrdNN[dds_resT6C6OrdNN$padj<0.05,])

All_TP6_TP5_names = row.names(dds_resT6T5Ord)
dds_resT6T5OrdNN <- dds_resT6T5Ord[!is.na(dds_resT6T5Ord$padj),]
Sig_TP6_TP5 = row.names(dds_resT6T5OrdNN[dds_resT6T5OrdNN$padj<0.05,])

Union_Sig_TP6CP6_TP6TP6 = intersect(Sig_TP6_CP6, Sig_TP6_TP5)

Only_Sig_TP6CP6 = setdiff(Sig_TP6_CP6, Union_Sig_TP6CP6_TP6TP6)
Only_Sig_TP6TP5 = setdiff(Sig_TP6_TP5, Union_Sig_TP6CP6_TP6TP6)

All_sig_genes_TP5CP6_TP6CP6 = unique(c(Union_Sig_TP6CP6_TP6TP6, Only_Sig_TP6CP6, Only_Sig_TP6TP5))
All_non_sig_genes_TP5CP6 = setdiff(All_TP6_CP6_names, All_sig_genes_TP5CP6_TP6CP6)

# Now write all needed .txt files
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6CP6_All_non_sig_genes_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Only_Sig_TP6CP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6CP6_Only_Sig_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Only_Sig_TP6TP5), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6CP6_Only_Sig_TP6TP5_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% Union_Sig_TP6CP6_TP6TP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6CP6_Union_Sig_TP5CP6_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")

write.table(subset(dds_resT6T5Ord, rownames(dds_resT6T5Ord) %in% All_non_sig_genes_TP5CP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6TP5_All_non_sig_genes_TP6CP6_TP6TP5_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6T5Ord, rownames(dds_resT6T5Ord) %in% Only_Sig_TP6CP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6TP5_Only_Sig_TP6CP6_170805.txt",col.names=T,quote=F,sep="\t")            
write.table(subset(dds_resT6T5Ord, rownames(dds_resT6T5Ord) %in% Only_Sig_TP6TP5), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6TP5_Only_Sig_TP6TP5_170805.txt",col.names=T,quote=F,sep="\t")
write.table(subset(dds_resT6T5Ord, rownames(dds_resT6T5Ord) %in% Union_Sig_TP6CP6_TP6TP6), 
            file="output/tables/XY_plots/TP6CP6_TP6TP5/TP6TP5_Union_Sig_TP6CP6_TP6TP5_170805.txt",col.names=T,quote=F,sep="\t")
