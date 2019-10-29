library(limma)
library(DESeq2)
library(gplots)
library(limma)
library(RColorBrewer)
library(pheatmap)

Szulzewsky_human_sets_C_T <- read.table(file="input/Szulzewsky_sets/Szulzewsky_human_sets_C_T.txt")

# Remove features (certain rownames) we no longer need
htseq_drop_rows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique","ERCC-")
drop <- rowSums(sapply(htseq_drop_rows, grepl, rownames(Szulzewsky_human_sets_C_T)))>0
Szulzewsky_human_sets_C_T <- Szulzewsky_human_sets_C_T[!drop,]

# Print dimensions of 
dim(Szulzewsky_human_sets_C_T)

# Keep genes with more than 1 read in at least 2 samples
keep <- rowSums(Szulzewsky_human_sets_C_T>=1)>=2
table(keep)
Szulzewsky_human_sets_C_T <- Szulzewsky_human_sets_C_T[keep,]
dim(Szulzewsky_human_sets_C_T)

# Drop low read count (less than 5 reads) samples. First show these
hist(colSums(Szulzewsky_human_sets_C_T>5), breaks=20)

# Print colSums output to spot the samples with less counts
colSums(Szulzewsky_human_sets_C_T>5)

# Drop samples where less than 6000 genes were detected with more than 5 reads
drop <- colSums(Szulzewsky_human_sets_C_T>5)<6000
table(drop)
Szulzewsky_human_sets_C_T <- Szulzewsky_human_sets_C_T[,!drop]

rm(drop, keep, htseq_drop_rows)

##
colDataD56 <- data.frame(condition=rep(NA,length(colnames(Szulzewsky_human_sets_C_T))), row.names=colnames(Szulzewsky_human_sets_C_T))
colDataD56$condition[grep("Tumor.",rownames(colDataD56))] <- "Tumor"
colDataD56$condition[grep("Postmortem.",rownames(colDataD56))] <- "Control"
dds <- DESeqDataSetFromMatrix(countData = Szulzewsky_human_sets_C_T, colData = colDataD56,design = ~condition)
dds$condition <- factor(dds$condition,levels=c("Tumor",
                                               "Control"))

dds <- DESeq(dds)


# For comparative and descriptive analysis we use the rlog. blind=FALSE to keep sample group info for calculating normalizations.  
rldMM <- rlog(dds, blind = FALSE)


## PART 4: First descriptive analysis

# Work on sample-to-sample heatmap
sampleDists <- dist(t(assay(rldMM)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)

pdf(file="output/pdf/Szulzewsky_sets/180530_sample-to-sample_heatmap.pdf", height = 24, width = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cellwidth = 16,
         cellheight = 16,
         border_color = NA,
         col=colors)
dev.off()

# Now do the DE analysis 

# load groups of genes as manually currated in the Groups_for_analysis.R file
source("R_scripts/shared/Groups_for_analysis_HUMAN.R")

# Comparison: Tumor versus Control
dds_res_Tumor_Postm <- results(dds, c("condition","Tumor","Control"))
dds_res_Tumor_Postm = dds_res_Tumor_Postm[order(row.names(dds_res_Tumor_Postm)),]
DESeq2::plotMA(dds_res_Tumor_Postm,  alpha = 0.05, main="MA plot Tumor vs. Postmortem", ylim=c(-10,10))

# Let's extract DE and ncount info on selected sets
normalized_counts <- counts(dds, normalized=TRUE)
colnames(normalized_counts)

# Sensome order TP5CP6
sensome_ordered = data.frame(subset(dds_res_Tumor_Postm, rownames(dds_res_Tumor_Postm) %in% HUMAN_microglia_sensome_TP5CP6))
sensome_ordered_order_TP5CP6 = sensome_ordered[match(HUMAN_microglia_sensome_TP5CP6, rownames(sensome_ordered)),]
write.table(sensome_ordered_order_TP5CP6,file="output/tables/Szulzewsky_sets/DE_subsets_C_T/DESeq_Tumor_Postm_HUMAN_microglia_sensome_ordered_order_TP5CP6_180704.txt",col.names=T,quote=F,sep="\t")