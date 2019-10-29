##
##    Heatmaps_rlog_values_microglia.R
##
##

## PART 1: Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the dataframe. We use the minus get microglia datasets file
source("R_scripts/shared/Datasets_microglia.R")

## PART 3: Setup DEseq2 analysis for all sample groups. 

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


# PART 4: extract rlog values and calculate averages per sample groups

# For comparative and descriptive analysis we use the rlog (see manual for descriptions). blind=FALSE to keep sample group info for calculating normalizations. 
rld56 <- rlog(dds, blind = FALSE)
rld56_values <- assay(rld56)

rld56_values <- rld56_values[,c("MC.P6", "FC1.P6", "FC2.P6", "MT.P6", "FT.P6", "FT1.P6", "MT.P5", "FT.P5", "FT1.P5")]

CP6_rlog_values <- rld56_values[, c("MC.P6", "FC1.P6", "FC2.P6")]
TP6_rlog_values <- rld56_values[, c("MT.P6", "FT.P6", "FT1.P6")]
TP5_rlog_values <- rld56_values[, c("MT.P5", "FT.P5", "FT1.P5")]

rm(colDataD56, cntNiek_microglia)

CP6_avg_rlog <- as.data.frame(CP6_rlog_values)
TP6_avg_rlog <- as.data.frame(TP6_rlog_values)
TP5_avg_rlog <- as.data.frame(TP5_rlog_values)

CP6_avg_rlog$mean_rlog = rowMeans(CP6_rlog_values[,c("MC.P6", "FC1.P6", "FC2.P6")], na.rm=FALSE)
TP6_avg_rlog$mean_rlog = rowMeans(TP6_rlog_values[,c("MT.P6", "FT.P6", "FT1.P6")], na.rm=FALSE)
TP5_avg_rlog$mean_rlog = rowMeans(TP5_rlog_values[,c("MT.P5", "FT.P5", "FT1.P5")], na.rm=FALSE)

# Only keep the column we need for the table
CP6_avg_rlog = CP6_avg_rlog[c("mean_rlog")]
TP6_avg_rlog = TP6_avg_rlog[c("mean_rlog")]
TP5_avg_rlog = TP5_avg_rlog[c("mean_rlog")]

rm(CP6_rlog_values, TP6_rlog_values, TP5_rlog_values, rld56)


## PART 5: Extract data for specific gene groups from the rlog values data obtained in part 4

# load groups of genes as manually currated in the Groups_for_analysis.R file
source("R_scripts/shared/Groups_for_analysis.R")

# mouse_sensome
microglia_sensome_updated_CP6 = subset(CP6_avg_rlog, rownames(CP6_avg_rlog) %in% microglia_sensome_updated)
microglia_sensome_updated_TP6 = subset(TP6_avg_rlog, rownames(TP6_avg_rlog) %in% microglia_sensome_updated)
microglia_sensome_updated_TP5 = subset(TP5_avg_rlog, rownames(TP5_avg_rlog) %in% microglia_sensome_updated)

# extract the average rlog values
CP6_rlog_mouse_sensome = microglia_sensome_updated_CP6$mean_rlog
TP6_rlog_mouse_sensome = microglia_sensome_updated_TP6$mean_rlog
TP5_rlog_mouse_sensome = microglia_sensome_updated_TP5$mean_rlog

mouse_sensome = data.frame(CP6_rlog_mouse_sensome, TP6_rlog_mouse_sensome, TP5_rlog_mouse_sensome)
row.names(mouse_sensome) = row.names(microglia_sensome_updated_CP6)

# Now obtain the right DE order for EV-GFPpos vs WT microglia
dds_resT5C6 <- results(dds, c("condition","TP5","CP6"))
dds_resT5C6Ord_log2fold <- dds_resT5C6[order(-dds_resT5C6$log2FoldChange),]

#Only the sensome
dds_resT5C6Ord_log2fold_sensome = subset(dds_resT5C6Ord_log2fold, rownames(dds_resT5C6Ord_log2fold) %in% microglia_sensome_updated)
sensome_order_T5C6 = rownames(dds_resT5C6Ord_log2fold_sensome)

mouse_sensome = mouse_sensome[match(sensome_order_T5C6, rownames(mouse_sensome)),]


## PART 7: Perform the actual plotting and saving of the files
rld56_values_sensome = subset(rld56_values, rownames(rld56_values) %in% microglia_sensome_updated)
rld56_values_sensome = as.data.frame(rld56_values_sensome)

rld56_values_sensome = rld56_values_sensome[match(sensome_order_T5C6, rownames(rld56_values_sensome)),]


## PART 6: Setup plotting variables for the heatmaps
my_palette  = colorRampPalette(c("blue4","white", "red"))(n = 32)
col_breaks  = c(seq(2,14, length=33)) 
colsep      = c(0,3,6)
sepwidth    = c(0.05,0.05)
sepcolor    = "white"
cexCol      = 1.0
cexRow      = 0.8
lmat        = rbind(c(0,3),c(2,1),c(0,4))
lwid        = c(1.5,1) # second position for column width
lhei        = c(1.5,2,0.3) # second position for row height


## PART 7: Perform the actual plotting and saving of the files
pdf(file="output/pdf/rlog_heatmaps/170823_mouse_sensome_sort_ALL_TP5CP6_just.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_sensome),
           main="mouse_sensome", # heat map title
           col=my_palette,
           #breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rld56_values_sensome)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           #keysize=0.05,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=-7.4,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()


# GSEA_mouse_TGFbeta1

dds_resT5C6Ord_log2fold_TGFb = subset(dds_resT5C6Ord_log2fold, rownames(dds_resT5C6Ord_log2fold) %in% GSEA_mouse_TGFbeta1)
TGFb_order_T5C6 = rownames(dds_resT5C6Ord_log2fold_TGFb)

rld56_values_TGFb = subset(rld56_values, rownames(rld56_values) %in% GSEA_mouse_TGFbeta1)
rld56_values_TGFb = as.data.frame(rld56_values_TGFb)

rld56_values_TGFb = rld56_values_TGFb[match(TGFb_order_T5C6, rownames(rld56_values_TGFb)),]

## Perform the actual plotting and saving of the files
pdf(file="output/pdf/rlog_heatmaps/170823_mouse_TGFb_sort_ALL_TP5CP6.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_TGFb),
           main="TGFb", # heat map title
           col=my_palette,
           #breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rld56_values_TGFb)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           #keysize=0.05,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=-7.4,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()