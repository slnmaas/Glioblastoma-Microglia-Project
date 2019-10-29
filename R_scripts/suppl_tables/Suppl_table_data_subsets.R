##
##    Suppl_table_data_all_genes.R
##    

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the dataframe. We use the minus get microglia datasets file
source("R_scripts/shared/Datasets_microglia.R")

## PART 3: Setup the DEseq2 analysis to obtain the normalized counts

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

## For the supplementary tables we work with the normalized gene counts
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- normalized_counts[,c("MC.P6", "FC1.P6", "FC2.P6", "MT.P6", "FT.P6", "FT1.P6", "MT.P5", "FT.P5", "FT1.P5")]

CP6_normalized_counts <- normalized_counts[, c("MC.P6", "FC1.P6", "FC2.P6")]
TP6_normalized_counts <- normalized_counts[, c("MT.P6", "FT.P6", "FT1.P6")]
TP5_normalized_counts <- normalized_counts[, c("MT.P5", "FT.P5", "FT1.P5")]

rm(normalized_counts)

##  PART 4: Calculate the normalized read counts, calculate the averages per group and write files

CP6_avg_ncounts <- as.data.frame(CP6_normalized_counts)
TP6_avg_ncounts <- as.data.frame(TP6_normalized_counts)
TP5_avg_ncounts <- as.data.frame(TP5_normalized_counts)

CP6_avg_ncounts$mean_ncounts = rowMeans(CP6_normalized_counts[,c("MC.P6", "FC1.P6", "FC2.P6")], na.rm=FALSE)
TP6_avg_ncounts$mean_ncounts = rowMeans(TP6_normalized_counts[,c("MT.P6", "FT.P6", "FT1.P6")], na.rm=FALSE)
TP5_avg_ncounts$mean_ncounts = rowMeans(TP5_normalized_counts[,c("MT.P5", "FT.P5", "FT1.P5")], na.rm=FALSE)

rm(CP6_normalized_counts, TP6_normalized_counts, TP5_normalized_counts)

# Only keep the column we need for the table
CP6_avg_ncounts = CP6_avg_ncounts[c("mean_ncounts")]
TP6_avg_ncounts = TP6_avg_ncounts[c("mean_ncounts")]
TP5_avg_ncounts = TP5_avg_ncounts[c("mean_ncounts")]

# Obtain the subgroups and filter data
## PART 2: Load GSEA gene sets
source("R_scripts/shared/Groups_for_analysis.R")

# The sets that we will be working with are: GSEA_mouse_TGFbeta1, microglia_sensome_updated

# GSEA_mouse_TGFbeta1
CP6_GSEA_mouse_TGFbeta1_ncounts = subset(CP6_avg_ncounts, rownames(CP6_avg_ncounts) %in% GSEA_mouse_TGFbeta1)
write.table(CP6_GSEA_mouse_TGFbeta1_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/GSEA_mouse_TGFbeta1/CP6_Avg_ncounts_GSEA_mouse_TGFbeta1_170730_ncounts.txt",col.names=T,quote=F,sep="\t")
TP6_GSEA_mouse_TGFbeta1_ncounts = subset(TP6_avg_ncounts, rownames(TP6_avg_ncounts) %in% GSEA_mouse_TGFbeta1)
write.table(TP6_GSEA_mouse_TGFbeta1_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/GSEA_mouse_TGFbeta1/TP6_Avg_ncounts_GSEA_mouse_TGFbeta1_170730_ncounts.txt",col.names=T,quote=F,sep="\t")
TP5_GSEA_mouse_TGFbeta1_ncounts = subset(TP5_avg_ncounts, rownames(TP5_avg_ncounts) %in% GSEA_mouse_TGFbeta1)
write.table(TP5_GSEA_mouse_TGFbeta1_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/GSEA_mouse_TGFbeta1/TP5_Avg_ncounts_GSEA_mouse_TGFbeta1_170730_ncounts.txt",col.names=T,quote=F,sep="\t")

# microglia_sensome_updated
CP6_microglia_sensome_updated_ncounts = subset(CP6_avg_ncounts, rownames(CP6_avg_ncounts) %in% microglia_sensome_updated)
write.table(CP6_microglia_sensome_updated_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/microglia_sensome_updated/CP6_Avg_ncounts_microglia_sensome_updated_170730_ncounts.txt",col.names=T,quote=F,sep="\t")
TP6_microglia_sensome_updated_ncounts = subset(TP6_avg_ncounts, rownames(TP6_avg_ncounts) %in% microglia_sensome_updated)
write.table(TP6_microglia_sensome_updated_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/microglia_sensome_updated/TP6_Avg_ncounts_microglia_sensome_updated_170730_ncounts.txt",col.names=T,quote=F,sep="\t")
TP5_microglia_sensome_updated_ncounts = subset(TP5_avg_ncounts, rownames(TP5_avg_ncounts) %in% microglia_sensome_updated)
write.table(TP5_microglia_sensome_updated_ncounts,file="output/tables/Suppl_tables/suppl_table2/avg_ncounts_subgroups/microglia_sensome_updated/TP5_Avg_ncounts_microglia_sensome_updated_170730_ncounts.txt",col.names=T,quote=F,sep="\t")



##  PART 5: Setup DEseq2 DE comparisons for all groups

# Set 1: DE comparing expression in CP6 to other groups
dds_resCP6_to_TP6 <- results(dds, c("condition","CP6","TP6"))
dds_resCP6_to_TP5 <- results(dds, c("condition","CP6","TP5"))

# Set 2: DE comparing expression in TP6 to other groups
dds_resTP6_to_CP6 <- results(dds, c("condition","TP6","CP6"))
dds_resTP6_to_TP5 <- results(dds, c("condition","TP6","TP5"))

# Set 3: DE comparing expression in TP5 to other groups
dds_resTP5_to_CP6 <- results(dds, c("condition","TP5","CP6"))
dds_resTP5_to_TP6 <- results(dds, c("condition","TP5","TP6"))


##  PART 6: Write all files for the different files to disk
# The sets that we will be working with are: GSEA_mouse_TGFbeta1, microglia_sensome_updated

# GSEA_mouse_TGFbeta1
dds_resCP6_to_TP6_GSEA_mouse_TGFbeta1 = subset(dds_resCP6_to_TP6, rownames(dds_resCP6_to_TP6) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resCP6_to_TP6_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resCP6_to_TP6_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")
dds_resCP6_to_TP5_GSEA_mouse_TGFbeta1 = subset(dds_resCP6_to_TP5, rownames(dds_resCP6_to_TP5) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resCP6_to_TP5_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resCP6_to_TP5_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")

dds_resTP6_to_CP6_GSEA_mouse_TGFbeta1 = subset(dds_resTP6_to_CP6, rownames(dds_resTP6_to_CP6) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resTP6_to_CP6_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resTP6_to_CP6_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")
dds_resTP6_to_TP5_GSEA_mouse_TGFbeta1 = subset(dds_resTP6_to_TP5, rownames(dds_resTP6_to_TP5) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resTP6_to_TP5_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resTP6_to_TP5_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")

dds_resTP5_to_CP6_GSEA_mouse_TGFbeta1 = subset(dds_resTP5_to_CP6, rownames(dds_resTP5_to_CP6) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resTP5_to_CP6_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resTP5_to_CP6_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")
dds_resTP5_to_TP6_GSEA_mouse_TGFbeta1 = subset(dds_resTP5_to_TP6, rownames(dds_resTP5_to_TP6) %in% GSEA_mouse_TGFbeta1)
write.table(dds_resTP5_to_TP6_GSEA_mouse_TGFbeta1,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/GSEA_mouse_TGFbeta1/dds_resTP5_to_TP6_GSEA_mouse_TGFbeta1_170730.txt",col.names=T,quote=F,sep="\t")

# microglia_sensome_updated
dds_resCP6_to_TP6_microglia_sensome_updated = subset(dds_resCP6_to_TP6, rownames(dds_resCP6_to_TP6) %in% microglia_sensome_updated)
write.table(dds_resCP6_to_TP6_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resCP6_to_TP6_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")
dds_resCP6_to_TP5_microglia_sensome_updated = subset(dds_resCP6_to_TP5, rownames(dds_resCP6_to_TP5) %in% microglia_sensome_updated)
write.table(dds_resCP6_to_TP5_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resCP6_to_TP5_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")

dds_resTP6_to_CP6_microglia_sensome_updated = subset(dds_resTP6_to_CP6, rownames(dds_resTP6_to_CP6) %in% microglia_sensome_updated)
write.table(dds_resTP6_to_CP6_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resTP6_to_CP6_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")
dds_resTP6_to_TP5_microglia_sensome_updated = subset(dds_resTP6_to_TP5, rownames(dds_resTP6_to_TP5) %in% microglia_sensome_updated)
write.table(dds_resTP6_to_TP5_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resTP6_to_TP5_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")

dds_resTP5_to_CP6_microglia_sensome_updated = subset(dds_resTP5_to_CP6, rownames(dds_resTP5_to_CP6) %in% microglia_sensome_updated)
write.table(dds_resTP5_to_CP6_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resTP5_to_CP6_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")
dds_resTP5_to_TP6_microglia_sensome_updated = subset(dds_resTP5_to_TP6, rownames(dds_resTP5_to_TP6) %in% microglia_sensome_updated)
write.table(dds_resTP5_to_TP6_microglia_sensome_updated,file="output/tables/Suppl_tables/suppl_table2/DE_subgroups/microglia_sensome_updated/dds_resTP5_to_TP6_microglia_sensome_updated_170730.txt",col.names=T,quote=F,sep="\t")
